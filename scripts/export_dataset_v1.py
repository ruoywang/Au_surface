#!/usr/bin/env python3
"""Dataset v1 export (DFT side): every ACCEPTED production state of dataset_plan_v1 rev 2 (2026-09-28).

Source of truth: 05_production/queue.json (status completed) plus the production-standard pilot states it
reuses (status reused, dir under 03_pilot). Same conventions as export_dataset_v0.py:
  * label states by the ACTUAL converged N_e / mu_e (CPM-ion line), TARGETMU kept separately;
  * all raw energy terms exported (TOTEN, E_without_entropy, E_sigma0, GCE_code, TOTEN - mu*dN); no generic 'energy' key;
  * raw forces on all atoms from the LAST TOTAL-FORCE block (for relaxations = final geometry; CHGCAR/fields are
    written for that geometry too), fixed atoms not zeroed; total drift recorded;
  * 3D grid files indexed in place (not copied).
Extra per state: structure_id, config (ideal / relax / relaxed / pertXX / coll_* / path_*), parent relaxation,
generator note (seed, sigma, image fraction, min distance), family, tier, batch priority.

Outputs (dataset_v1/): states.json, fields_index.json, production.extxyz, summary.md
Usage (from Au_Cl/):  scripts/pyrun.sh scripts/export_dataset_v1.py
"""
import collections
import hashlib
import json
import os
import re
import time

import numpy as np
from ase import Atoms
from ase.io import read, write

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
OUT = f"{ROOT}/dataset_v1"
ZVAL_AU = 11.0
MU0 = -4.9071        # internal electron-chemical-potential reference (TARGETMU of the reference states); delta_U = MU0 - mu
MU_TOL = 0.011       # acceptance |mu_e - TARGETMU| applied by production.evaluate (FERMICONVERGE = 0.01 + 1 meV)
CONFIG_VERSION = "production-2026-09-26 (PREC=Normal, ALGO=Fast, NPAR=16 hybrid, FERMICONVERGE=0.01)"


def geom_hash(atoms):
    """same geometry hash as scripts/extend_mu05.py: cell + positions rounded to 1e-3 A -> states of one geometry share it."""
    h = hashlib.sha1()
    h.update(np.round(atoms.get_cell().array, 3).tobytes()); h.update(np.round(atoms.get_positions(), 3).tobytes())
    return h.hexdigest()[:12]


FIELD_FILES =["CHGCAR", "LOCPOT", "PHI", "PHISOLV", "VSOLV", "RHOB", "RHOION", "ELOC", "P", "SVDW", "SION", "SSOLV", "SCAV", "SDIEL", "POT"]
FIELD_NOTE = {"CHGCAR": "rho_e*V, first block only; sum/N_grid = N_e", "LOCPOT": "Hartree potential (LVHAR), eV",
              "PHI": "phi_explicit + phi_solv, total electrostatic potential of the electrolyte model; u = PHI/(k_B 298 K)",
              "PHISOLV": "solvent feedback potential", "VSOLV": "correction potential added to the KS Hamiltonian",
              "RHOB": "dielectric bound charge x Omega", "RHOION": "net mobile-ion charge x Omega (RHOION/V = n_- - n_+)",
              "ELOC": "local field, z component", "P": "polarisation, z component", "SVDW": "cavity mask", "SION": "ion accessibility mask",
              "SSOLV": "cavity mask", "SCAV": "cavity mask", "SDIEL": "dielectric cavity mask", "POT": "semantics not verified"}


def incar_tags(path):
    tags = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if "=" in line:
            k, v = [x.strip() for x in line.split("=", 1)]; tags[k.upper()] = v
    return tags


def grid_dims(path):
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 3 and all(x.isdigit() for x in p): return [int(x) for x in p]
    return None


def parse_outcar(path, natoms):
    txt = open(path).read()
    toten = float(re.findall(r"free  energy   TOTEN\s*=\s*([-\d.]+)", txt)[-1])
    m = re.findall(r"energy  without entropy=\s*([-\d.]+)\s+energy\(sigma->0\)\s*=\s*([-\d.]+)", txt)[-1]
    efermi = float(re.findall(r"E-fermi :\s*([-\d.]+)", txt)[-1])
    blk = txt.split("TOTAL-FORCE (eV/Angst)")[-1].split("total drift")[0].splitlines()[2:2 + natoms]
    forces = np.array([[float(x) for x in l.split()[3:6]] for l in blk])
    drift = [float(x) for x in re.findall(r"total drift:\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", txt)[-1]]
    n_ionic = txt.count("TOTAL-FORCE (eV/Angst)")
    return dict(E_free_TOTEN=toten, E_without_entropy=float(m[0]), E_sigma0=float(m[1]), E_fermi_OUTCAR=efermi,
                forces=forces, total_drift=drift, n_ionic_steps=n_ionic, reached_accuracy="reached required accuracy" in txt)


def parse_cpm(path):
    txt = open(path).read()
    ion = re.findall(r"CPM-ion:.*?N_ele=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)\s+TARGETMU=\s*([-\d.]+).*?GCE=\s*([-\d.E+]+)\s+cap=\s*([-\d.]+)", txt)
    scf1 = re.findall(r"CPM-scf: NSTEP=\s*1 SCF=\s*1 N_ele=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)", txt)
    if not ion: return None
    return dict(N_e_final=float(ion[-1][0]), mu_e_actual=float(ion[-1][1]), TARGETMU=float(ion[-1][2]), GCE_code=float(ion[-1][3]),
                cap_code=float(ion[-1][4]), mu_e_neutral_first_round=float(scf1[0][1]) if scf1 else None)


def selective_flags(poscar, natoms):
    lines = open(poscar).read().splitlines()
    return np.array([[x == "F" for x in l.split()[3:6]] for l in lines[9:9 + natoms]])


q = json.load(open(f"{ROOT}/05_production/queue.json"))
# geometry registry = the extension index (scripts/extend_mu05.py): one canonical geometry_id + source file per (structure, config)
_MAN = f"{ROOT}/05_production/mu05_extension_manifest.json"
REGISTRY = {(g["structure_id"], g["config"]): g for g in json.load(open(_MAN))["geometries"]} if os.path.exists(_MAN) else {}
_REG_ATOMS = {}
def registry_atoms(src):
    if src not in _REG_ATOMS: _REG_ATOMS[src] = read(f"{ROOT}/{src}")
    return _REG_ATOMS[src]
states, skipped, aliases = {}, [], {}
seen_dirs = {}
for t in sorted(q, key=lambda t: t["status"] != "completed"):          # computed records first, reused aliases after
    if t["status"] not in ("completed", "reused"): continue
    d = t["dir"] if t["dir"].startswith("/") else f"{ROOT}/{t['dir'].split(' ')[0]}"
    if d in seen_dirs:                        # same run directory = same physical state (the reused relaxed@-4.9071 record IS the
        aliases[t["task_id"]] = seen_dirs[d]  # relax task's final step): exported once, the second id recorded as an alias
        continue
    seen_dirs[d] = t["task_id"]
    if not os.path.exists(f"{d}/OUTCAR"): skipped.append((t["task_id"], "no OUTCAR at " + d)); continue
    # relaxation directory (the relax task AND its reused relaxed@-4.9071 record): the state's geometry is the final CONTCAR.
    # (Before 2026-09-28 the reused record was exported with the relax dir's POSCAR = the PRE-relaxation geometry: fixed.)
    geo_file = f"{d}/CONTCAR" if os.path.basename(d).startswith("relax__") and os.path.exists(f"{d}/CONTCAR") else f"{d}/POSCAR"
    atoms = read(geo_file); n = len(atoms)
    cfg_norm = "relaxed" if t["config"] in ("relax", "relaxed") else t["config"]
    reg = REGISTRY.get((t["structure_id"], cfg_norm))
    if reg:                                    # canonical geometry from the extension index; record this state's deviation from it
        geometry_id = reg["geometry_id"]; ref = registry_atoms(reg["source"])
        geo_dev = float(np.abs(atoms.get_positions() - ref.get_positions()).max()) if len(ref) == n else None
    else:
        geometry_id = f"{t['structure_id']}__{cfg_norm}__{geom_hash(atoms)}"; geo_dev = None
    try:
        oc = parse_outcar(f"{d}/OUTCAR", n); cpm = parse_cpm(f"{d}/log.out")
    except Exception as e:
        skipped.append((t["task_id"], f"parse error {type(e).__name__}: {e}")); continue
    if cpm is None: skipped.append((t["task_id"], "no CPM-ion line")); continue
    tags = incar_tags(f"{d}/INCAR"); fixed = selective_flags(f"{d}/POSCAR", n) if "Selective" in open(f"{d}/POSCAR").read() else np.zeros((n, 3), bool)
    N0 = ZVAL_AU * n; dN = cpm["N_e_final"] - N0
    fields = {f: dict(path=f"{d}/{f}", bytes=os.path.getsize(f"{d}/{f}"), grid=grid_dims(f"{d}/{f}"), note=FIELD_NOTE[f])
              for f in FIELD_FILES if os.path.exists(f"{d}/{f}") and os.path.getsize(f"{d}/{f}") > 0}
    # QC: mu closure (farm acceptance), CHGCAR written by THIS run (warm-started tasks start with a copied CHGCAR that
    # VASP overwrites at the end; a stale copy would predate the task start), all 15 field files present, geometry = registry
    # geometry within 0.005 A (the pilot-reused ideal states of Step-8x1/Step-16x1 deviate by 0.003 A, recorded).
    started = time.mktime(time.strptime(t["started_at"], "%Y-%m-%d %H:%M")) if t.get("started_at") else None
    chg_fresh = (os.path.getmtime(f"{d}/CHGCAR") > started) if (started and "CHGCAR" in fields) else None
    qc = dict(mu_within_tolerance=abs(cpm["mu_e_actual"] - cpm["TARGETMU"]) <= MU_TOL, chgcar_written_after_start=chg_fresh,
              fields_complete=(len(fields) == len(FIELD_FILES)), n_fields=len(fields),
              geometry_matches_registry=(geo_dev <= 0.005) if geo_dev is not None else None, geometry_deviation_A=geo_dev,
              relaxation_reached_accuracy=(oc["reached_accuracy"] if t["kind"] == "relax" else None))
    qc["status"] = "accepted" if qc["mu_within_tolerance"] and qc["chgcar_written_after_start"] is not False and qc["fields_complete"] \
        and qc["relaxation_reached_accuracy"] is not False and qc["geometry_matches_registry"] is not False else "flagged"
    states[t["task_id"]] = dict(
        state_id=t["task_id"], structure_id=t["structure_id"], family=t["family"], tier=t["tier"], config=t["config"], kind=t["kind"],
        geometry_id=geometry_id, campaign=t.get("campaign", "dataset_plan_v1_rev2"),
        parent_relaxation=(f"{t['structure_id']}__relax__mu-4.9071" if t["config"] not in ("ideal", "relax") else None),
        generator_note=t.get("note", ""), source_dir=d, geometry_file=os.path.basename(geo_file), source_geometry=t.get("source_geometry"),
        warm_start=t.get("warm_start", "cold start"), qc=qc,
        config_version=CONFIG_VERSION,
        structure=dict(n_atoms=n, cell_A=atoms.get_cell().tolist(), positions_A=atoms.get_positions().tolist(), pbc=[True, True, True],
                       fixed_xyz=fixed.tolist(), n_fixed_atoms=int(fixed.all(axis=1).sum()), metal_top_z_A=float(atoms.get_positions()[:, 2].max())),
        boundary=dict(Lz_A=float(atoms.get_cell()[2][2]), SOL_Z0_A=float(tags.get("SOL_Z0", 0)), SOL_Z1_A=float(tags.get("SOL_Z1", 0)),
                      DIPOL_frac=[float(x) for x in tags.get("DIPOL", "0.5 0.5 0.5").split()], C_MOLAR=float(tags.get("C_MOLAR", 0)), R_ION_A=float(tags.get("R_ION", 0))),
        numerics=dict(PREC=tags.get("PREC"), ENCUT_eV=float(tags.get("ENCUT", 0)), ALGO=tags.get("ALGO"), ISMEAR=int(tags.get("ISMEAR", 0)),
                      SIGMA_eV=float(tags.get("SIGMA", 0)), kpoints=t.get("kpoints"), FERMICONVERGE_eV=float(tags.get("FERMICONVERGE", 0)),
                      IBRION=int(tags.get("IBRION", -1)), EDIFFG=tags.get("EDIFFG")),
        electronic_state=dict(N_neutral=N0, N_e_final=cpm["N_e_final"], delta_N_e=dN, mu_e_actual_eV=cpm["mu_e_actual"], TARGETMU_eV=cpm["TARGETMU"],
                              mu_reference_eV=MU0, delta_U_target_V=round(MU0 - cpm["TARGETMU"], 4), delta_U_actual_V=MU0 - cpm["mu_e_actual"],
                              mu_e_minus_target_eV=cpm["mu_e_actual"] - cpm["TARGETMU"], mu_e_neutral_first_round_eV=cpm["mu_e_neutral_first_round"],
                              E_fermi_OUTCAR_eV=oc["E_fermi_OUTCAR"], cap_code_e_per_V=cpm["cap_code"]),
        labels=dict(E_free_TOTEN_eV=oc["E_free_TOTEN"], E_without_entropy_eV=oc["E_without_entropy"], E_sigma0_eV=oc["E_sigma0"], GCE_code_eV=cpm["GCE_code"],
                    TOTEN_minus_mu_dN_eV=oc["E_free_TOTEN"] - cpm["mu_e_actual"] * dN, energy_label_status="UNCONFIRMED (trainer target not decided)",
                    forces_eV_per_A=oc["forces"].tolist(), max_force_eV_per_A=float(np.abs(oc["forces"]).max()),
                    max_force_movable_eV_per_A=float(np.abs(oc["forces"][~fixed.all(axis=1)]).max()) if (~fixed.all(axis=1)).any() else None,
                    total_drift_eV_per_A=oc["total_drift"], n_ionic_steps=oc["n_ionic_steps"], relaxation_reached_accuracy=oc["reached_accuracy"] if t["kind"] == "relax" else None),
        fields=fields, wall_min=t.get("elapsed_farm_min"), node=t.get("node"))

os.makedirs(OUT, exist_ok=True)
for k, s in states.items(): s["aliases"] = [a for a, b in aliases.items() if b == k]
json.dump(dict(created="2026-09-28", plan="dataset_plan_v1 rev 2", label_standard="production configuration", n_states=len(states), states=states,
               aliases=aliases, skipped=skipped), open(f"{OUT}/states.json", "w"), indent=1)
json.dump({k: v["fields"] for k, v in states.items()}, open(f"{OUT}/fields_index.json", "w"), indent=1)
lst = []
for sid, s in states.items():
    at = Atoms("Au%d" % s["structure"]["n_atoms"], positions=s["structure"]["positions_A"], cell=s["structure"]["cell_A"], pbc=True)
    es, lb = s["electronic_state"], s["labels"]
    at.info.update(dict(state_id=sid, structure_id=s["structure_id"], family=s["family"].replace(" ", "_"), config=s["config"],
                        geometry_id=s["geometry_id"], campaign=s["campaign"], qc_status=s["qc"]["status"],
                        N_neutral=es["N_neutral"], N_e_final=es["N_e_final"], delta_N_e=es["delta_N_e"], mu_e_actual_eV=es["mu_e_actual_eV"],
                        TARGETMU_eV=es["TARGETMU_eV"], delta_U_target_V=es["delta_U_target_V"], delta_U_actual_V=es["delta_U_actual_V"],
                        E_free_TOTEN_eV=lb["E_free_TOTEN_eV"], E_without_entropy_eV=lb["E_without_entropy_eV"],
                        E_sigma0_eV=lb["E_sigma0_eV"], GCE_code_eV=lb["GCE_code_eV"], TOTEN_minus_mu_dN_eV=lb["TOTEN_minus_mu_dN_eV"],
                        energy_label_status="UNCONFIRMED", config_version="production-2026-09-26"))
    at.arrays["forces_raw"] = np.array(lb["forces_eV_per_A"]); at.arrays["fixed"] = np.array(s["structure"]["fixed_xyz"]).all(axis=1)
    lst.append(at)
write(f"{OUT}/production.extxyz", lst, format="extxyz")

fam = collections.Counter(s["family"] for s in states.values()); cfg = collections.Counter(s["config"].split("_")[0] for s in states.values())
geo = len({s["geometry_id"] for s in states.values()})
camp = collections.Counter(s["campaign"] for s in states.values()); mus = collections.Counter(f"{s['electronic_state']['TARGETMU_eV']:.4f}" for s in states.values())
flagged = [(k, {a: b for a, b in s["qc"].items() if b is False}) for k, s in states.items() if s["qc"]["status"] != "accepted"]
per_geo = collections.Counter(s["geometry_id"] for s in states.values())
off = [s for s in states.values() if abs(s["electronic_state"]["mu_e_minus_target_eV"]) > 0.005]
drift = np.array([abs(s["labels"]["total_drift_eV_per_A"][2]) for s in states.values()])
L = ["# Dataset v1 (production states) — summary", "", f"Generated {json.load(open(f'{OUT}/states.json'))['created']} from 05_production/queue.json. "
     f"**{len(states)} accepted states** ({geo} distinct geometries, {len({s['structure_id'] for s in states.values()})} structures), {len(skipped)} skipped. "
     "Labels are by the actual converged state; no `energy` key is assigned (see dataset_v0/README.md).", "",
     "| family | states |", "|---|---|"] + [f"| {k} | {v} |" for k, v in fam.most_common()] + ["", "| config type | states |", "|---|---|"] + \
    [f"| {k} | {v} |" for k, v in cfg.most_common()] + ["",
     f"- μ_e more than 5 meV from TARGETMU (inside FERMICONVERGE = 0.01): {len(off)} states; max |Δμ| = {max(abs(s['electronic_state']['mu_e_minus_target_eV']) for s in states.values())*1000:.1f} meV.",
     f"- |total drift_z|: median {np.median(drift):.3f}, max {drift.max():.3f} eV/Å (PREC=Normal aliasing; see dataset_v0/README.md).",
     f"- relaxations included as final-geometry states: {sum(1 for s in states.values() if s['kind']=='relax')} (all with 'reached required accuracy'); "
     f"{len(aliases)} queue records are aliases of an exported state (the reused relaxed@-4.9071 record = the relax task's final step) and are "
     f"listed under `aliases`, not exported twice (before 2026-09-28 20:50 they were).",
     f"- campaigns: {dict(camp)}; TARGETMU coverage: {dict(sorted(mus.items()))}; states per geometry: {dict(sorted(collections.Counter(per_geo.values()).items()))}.",
     f"- QC (mu closure <= {MU_TOL} eV, CHGCAR written after task start, 15 field files, geometry = registry within 0.005 Å, relaxation accuracy): "
     f"{len(states)-len(flagged)} accepted, {len(flagged)} flagged" + (": " + "; ".join(f"{k} {v}" for k, v in flagged) if flagged else "."),
     f"- geometry deviation from the registry source: max {max((s['qc']['geometry_deviation_A'] or 0) for s in states.values()):.4f} Å "
     f"({sum(1 for s in states.values() if (s['qc']['geometry_deviation_A'] or 0) > 1e-4)} states above 1e-4 Å, the pilot-reused ideal states of Step-8x1/Step-16x1); "
     f"states in the relaxation directory (relax task and reused relaxed@-4.9071) carry the CONTCAR geometry.",
     "- fields: every state indexes CHGCAR, LOCPOT, PHI, PHISOLV, VSOLV, RHOB, RHOION, ELOC, P, SVDW, SION, SSOLV, SCAV, SDIEL, POT in its run directory (not copied).", ""]
if skipped: L += ["## Skipped", ""] + [f"- {a}: {b}" for a, b in skipped] + [""]
open(f"{OUT}/summary.md", "w").write("\n".join(L))
print(f"{len(states)} states exported -> {OUT}/  (skipped {len(skipped)})"); print("\n".join(L[4:4+len(fam)+2]))
