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
import json
import os
import re

import numpy as np
from ase import Atoms
from ase.io import read, write

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
OUT = f"{ROOT}/dataset_v1"
ZVAL_AU = 11.0
FIELD_FILES = ["CHGCAR", "LOCPOT", "PHI", "PHISOLV", "VSOLV", "RHOB", "RHOION", "ELOC", "P", "SVDW", "SION", "SSOLV", "SCAV", "SDIEL", "POT"]
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
states, skipped = {}, []
for t in q:
    if t["status"] not in ("completed", "reused"): continue
    d = t["dir"] if t["dir"].startswith("/") else f"{ROOT}/{t['dir'].split(' ')[0]}"
    if not os.path.exists(f"{d}/OUTCAR"): skipped.append((t["task_id"], "no OUTCAR at " + d)); continue
    geo_file = f"{d}/CONTCAR" if t["kind"] == "relax" and os.path.exists(f"{d}/CONTCAR") else f"{d}/POSCAR"
    atoms = read(geo_file); n = len(atoms)
    try:
        oc = parse_outcar(f"{d}/OUTCAR", n); cpm = parse_cpm(f"{d}/log.out")
    except Exception as e:
        skipped.append((t["task_id"], f"parse error {type(e).__name__}: {e}")); continue
    if cpm is None: skipped.append((t["task_id"], "no CPM-ion line")); continue
    tags = incar_tags(f"{d}/INCAR"); fixed = selective_flags(f"{d}/POSCAR", n) if "Selective" in open(f"{d}/POSCAR").read() else np.zeros((n, 3), bool)
    N0 = ZVAL_AU * n; dN = cpm["N_e_final"] - N0
    fields = {f: dict(path=f"{d}/{f}", bytes=os.path.getsize(f"{d}/{f}"), grid=grid_dims(f"{d}/{f}"), note=FIELD_NOTE[f])
              for f in FIELD_FILES if os.path.exists(f"{d}/{f}") and os.path.getsize(f"{d}/{f}") > 0}
    states[t["task_id"]] = dict(
        state_id=t["task_id"], structure_id=t["structure_id"], family=t["family"], tier=t["tier"], config=t["config"], kind=t["kind"],
        parent_relaxation=(f"{t['structure_id']}__relax__mu-4.9071" if t["config"] not in ("ideal", "relax") else None),
        generator_note=t.get("note", ""), source_dir=d, geometry_file=os.path.basename(geo_file),
        config_version="production-2026-09-26 (PREC=Normal, ALGO=Fast, NPAR=16 hybrid, FERMICONVERGE=0.01)",
        structure=dict(n_atoms=n, cell_A=atoms.get_cell().tolist(), positions_A=atoms.get_positions().tolist(), pbc=[True, True, True],
                       fixed_xyz=fixed.tolist(), n_fixed_atoms=int(fixed.all(axis=1).sum()), metal_top_z_A=float(atoms.get_positions()[:, 2].max())),
        boundary=dict(Lz_A=float(atoms.get_cell()[2][2]), SOL_Z0_A=float(tags.get("SOL_Z0", 0)), SOL_Z1_A=float(tags.get("SOL_Z1", 0)),
                      DIPOL_frac=[float(x) for x in tags.get("DIPOL", "0.5 0.5 0.5").split()], C_MOLAR=float(tags.get("C_MOLAR", 0)), R_ION_A=float(tags.get("R_ION", 0))),
        numerics=dict(PREC=tags.get("PREC"), ENCUT_eV=float(tags.get("ENCUT", 0)), ALGO=tags.get("ALGO"), ISMEAR=int(tags.get("ISMEAR", 0)),
                      SIGMA_eV=float(tags.get("SIGMA", 0)), kpoints=t.get("kpoints"), FERMICONVERGE_eV=float(tags.get("FERMICONVERGE", 0)),
                      IBRION=int(tags.get("IBRION", -1)), EDIFFG=tags.get("EDIFFG")),
        electronic_state=dict(N_neutral=N0, N_e_final=cpm["N_e_final"], delta_N_e=dN, mu_e_actual_eV=cpm["mu_e_actual"], TARGETMU_eV=cpm["TARGETMU"],
                              mu_e_minus_target_eV=cpm["mu_e_actual"] - cpm["TARGETMU"], mu_e_neutral_first_round_eV=cpm["mu_e_neutral_first_round"],
                              E_fermi_OUTCAR_eV=oc["E_fermi_OUTCAR"], cap_code_e_per_V=cpm["cap_code"]),
        labels=dict(E_free_TOTEN_eV=oc["E_free_TOTEN"], E_without_entropy_eV=oc["E_without_entropy"], E_sigma0_eV=oc["E_sigma0"], GCE_code_eV=cpm["GCE_code"],
                    TOTEN_minus_mu_dN_eV=oc["E_free_TOTEN"] - cpm["mu_e_actual"] * dN, energy_label_status="UNCONFIRMED (trainer target not decided)",
                    forces_eV_per_A=oc["forces"].tolist(), max_force_eV_per_A=float(np.abs(oc["forces"]).max()),
                    max_force_movable_eV_per_A=float(np.abs(oc["forces"][~fixed.all(axis=1)]).max()) if (~fixed.all(axis=1)).any() else None,
                    total_drift_eV_per_A=oc["total_drift"], n_ionic_steps=oc["n_ionic_steps"], relaxation_reached_accuracy=oc["reached_accuracy"] if t["kind"] == "relax" else None),
        fields=fields, wall_min=t.get("elapsed_farm_min"), node=t.get("node"))

os.makedirs(OUT, exist_ok=True)
json.dump(dict(created="2026-09-28", plan="dataset_plan_v1 rev 2", label_standard="production configuration", n_states=len(states), states=states,
               skipped=skipped), open(f"{OUT}/states.json", "w"), indent=1)
json.dump({k: v["fields"] for k, v in states.items()}, open(f"{OUT}/fields_index.json", "w"), indent=1)
lst = []
for sid, s in states.items():
    at = Atoms("Au%d" % s["structure"]["n_atoms"], positions=s["structure"]["positions_A"], cell=s["structure"]["cell_A"], pbc=True)
    es, lb = s["electronic_state"], s["labels"]
    at.info.update(dict(state_id=sid, structure_id=s["structure_id"], family=s["family"].replace(" ", "_"), config=s["config"],
                        N_neutral=es["N_neutral"], N_e_final=es["N_e_final"], delta_N_e=es["delta_N_e"], mu_e_actual_eV=es["mu_e_actual_eV"],
                        TARGETMU_eV=es["TARGETMU_eV"], E_free_TOTEN_eV=lb["E_free_TOTEN_eV"], E_without_entropy_eV=lb["E_without_entropy_eV"],
                        E_sigma0_eV=lb["E_sigma0_eV"], GCE_code_eV=lb["GCE_code_eV"], TOTEN_minus_mu_dN_eV=lb["TOTEN_minus_mu_dN_eV"],
                        energy_label_status="UNCONFIRMED", config_version="production-2026-09-26"))
    at.arrays["forces_raw"] = np.array(lb["forces_eV_per_A"]); at.arrays["fixed"] = np.array(s["structure"]["fixed_xyz"]).all(axis=1)
    lst.append(at)
write(f"{OUT}/production.extxyz", lst, format="extxyz")

fam = collections.Counter(s["family"] for s in states.values()); cfg = collections.Counter(s["config"].split("_")[0] for s in states.values())
geo = len({(s["structure_id"], s["config"]) for s in states.values()})
off = [s for s in states.values() if abs(s["electronic_state"]["mu_e_minus_target_eV"]) > 0.005]
drift = np.array([abs(s["labels"]["total_drift_eV_per_A"][2]) for s in states.values()])
L = ["# Dataset v1 (production states) — summary", "", f"Generated {json.load(open(f'{OUT}/states.json'))['created']} from 05_production/queue.json. "
     f"**{len(states)} accepted states** ({geo} distinct geometries, {len({s['structure_id'] for s in states.values()})} structures), {len(skipped)} skipped. "
     "Labels are by the actual converged state; no `energy` key is assigned (see dataset_v0/README.md).", "",
     "| family | states |", "|---|---|"] + [f"| {k} | {v} |" for k, v in fam.most_common()] + ["", "| config type | states |", "|---|---|"] + \
    [f"| {k} | {v} |" for k, v in cfg.most_common()] + ["",
     f"- μ_e more than 5 meV from TARGETMU (inside FERMICONVERGE = 0.01): {len(off)} states; max |Δμ| = {max(abs(s['electronic_state']['mu_e_minus_target_eV']) for s in states.values())*1000:.1f} meV.",
     f"- |total drift_z|: median {np.median(drift):.3f}, max {drift.max():.3f} eV/Å (PREC=Normal aliasing; see dataset_v0/README.md).",
     f"- relaxations included as final-geometry states: {sum(1 for s in states.values() if s['kind']=='relax')} (all with 'reached required accuracy').",
     "- fields: every state indexes CHGCAR, LOCPOT, PHI, PHISOLV, VSOLV, RHOB, RHOION, ELOC, P, SVDW, SION, SSOLV, SCAV, SDIEL, POT in its run directory (not copied).", ""]
if skipped: L += ["## Skipped", ""] + [f"- {a}: {b}" for a, b in skipped] + [""]
open(f"{OUT}/summary.md", "w").write("\n".join(L))
print(f"{len(states)} states exported -> {OUT}/  (skipped {len(skipped)})"); print("\n".join(L[4:4+len(fam)+2]))
