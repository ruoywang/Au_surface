#!/usr/bin/env python3
"""Dataset v0 export (DFT side), 2026-09-27.

Turns the completed CP-DFT single points under 03_pilot/ into a small, self-describing
package that a trainer can read, WITHOUT deciding anything the trainer must decide:

  * groups:  v0_main_production  -- PREC=Normal/ALGO=Fast production configuration (the v0 label standard)
             reference_accurate  -- PREC=Accurate pilot configuration, kept as a separate reference group
             excluded            -- timed-out / cancelled / non-CP tests / speed tests, with the reason
    Same geometry at the same TARGETMU in both groups = configuration duplicate (flagged, NOT extra coverage).
    Same geometry at different actual charge = separate states sharing one geometry_id.
  * per state: structure + fixed-atom mask + boundary/solvent window + numerical config;
    ACTUAL electronic state (N_e, mu_e reached) separate from the TARGET (TARGETMU, neutral N0);
    every raw energy term (TOTEN free energy, energy without entropy, energy(sigma->0), the
    code's GCE and what it numerically is); raw forces on ALL atoms (fixed atoms not zeroed);
    stress and total drift as reported; index of every 3D grid file with grid dims, units and
    normalisation notes -- the grids themselves stay in the run directories (150-400 MB each).
  * NO 'energy' key is assigned: which energy term is the trainer's target (and its
    consistency with the forces) is recorded as UNCONFIRMED (energy_label_status) and is
    selected explicitly at training time via the trainer's key options.  Analysis profiles
    (K_D, dGamma) are NOT included -- they are evaluation quantities, not labels.

Outputs (dataset_v0/):  states.json, fields_index.json, selection.md, README.md,
  v0_main_production.extxyz, reference_accurate.extxyz
Usage (from Au_Cl/):  scripts/pyrun.sh scripts/export_dataset_v0.py
"""
import glob
import json
import os
import re

import numpy as np
from ase import Atoms
from ase.io import read, write

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
OUT = f"{ROOT}/dataset_v0"
ZVAL_AU = 11.0
KBT_CODE = 8.6173857e-5 * 298.0

FIELD_SEMANTICS = {
    "CHGCAR": "electron density rho_e(r)*V on the fine FFT grid (first block only; PAW augmentation block follows). sum/N_grid = N_e electrons. Positive = electrons.",
    "LOCPOT": "LVHAR=.TRUE.: Hartree (electrostatic) potential only, eV, fine grid.",
    "PHI": "phi = phi_explicit + phi_solv (solvation.F:1846): total electrostatic potential of the electrolyte model that drives the ion response; u = PHI/(k_B*298 K) with k_B = 8.6173857e-5 eV/K in the code. No volume factor. Sign: more negative -> anion enrichment.",
    "PHISOLV": "phi_solv (solvation.F:1848): solvent feedback potential alone.",
    "VSOLV": "V_corr (solvation.F:1832): correction potential actually added to the KS Hamiltonian.",
    "RHOB": "dielectric bound charge n_b(r) multiplied by the cell volume Omega (solvation.F:1375/1851): divide by V for e/A^3.",
    "RHOION": "net mobile-ion charge density multiplied by Omega (solvation.F:1384/1852): RHOION/V = n_anion - n_cation in ions/A^3 (validated against the PHI/SION reconstruction, relL2 ~1e-11).",
    "ELOC": "local electric field, z component only (solvation.F:1849).",
    "P": "polarisation, z component only (solvation.F:1850).",
    "SVDW": "cavity/accessibility mask in [0,1] (solvation.F:2048).",
    "SION": "ion accessibility mask in [0,1] (Stern layer excluded), solvation.F:2067.",
    "SSOLV": "cavity/accessibility mask in [0,1] (solvation.F:2085).",
    "SCAV": "cavity/accessibility mask in [0,1] (solvation.F:2104).",
    "SDIEL": "dielectric cavity mask in [0,1] (solvation.F:2124).",
    "POT": "present in the run directories; semantics NOT verified in the source audit (parameter_map.md section C does not cover it).",
}
FIELD_FILES = list(FIELD_SEMANTICS)


def incar_tags(path):
    tags = {}
    for line in open(path):
        line = line.split("#")[0].split("!")[0].strip()
        if "=" in line:
            k, v = [x.strip() for x in line.split("=", 1)]
            tags[k.upper()] = v
    return tags


def grid_dims(path):
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 3 and all(x.isdigit() for x in p):
                return [int(x) for x in p]
    return None


def parse_outcar(path, natoms):
    txt = open(path).read()
    toten = float(re.findall(r"free  energy   TOTEN\s*=\s*([-\d.]+)", txt)[-1])
    m = re.findall(r"energy  without entropy=\s*([-\d.]+)\s+energy\(sigma->0\)\s*=\s*([-\d.]+)", txt)[-1]
    efermi = float(re.findall(r"E-fermi :\s*([-\d.]+)", txt)[-1])
    blk = txt.split("TOTAL-FORCE (eV/Angst)")[-1].split("total drift")[0].splitlines()[2:2 + natoms]
    forces = np.array([[float(x) for x in l.split()[3:6]] for l in blk])
    drift = [float(x) for x in re.findall(r"total drift:\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", txt)[-1]]
    st = re.findall(r"in kB\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", txt)
    stress_kB = [float(x) for x in st[-1]] if st else None
    return dict(E_free_TOTEN=toten, E_without_entropy=float(m[0]), E_sigma0=float(m[1]), E_fermi_OUTCAR=efermi,
                forces=forces, total_drift=drift, stress_kB_as_printed=stress_kB)


def parse_cpm(path):
    txt = open(path).read()
    ion = re.findall(r"CPM-ion:.*?N_ele=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)\s+TARGETMU=\s*([-\d.]+).*?GCE=\s*([-\d.E+]+)\s+cap=\s*([-\d.]+)", txt)
    scf1 = re.findall(r"CPM-scf: NSTEP=\s*1 SCF=\s*1 N_ele=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)", txt)
    if not ion:
        return None
    return dict(N_e_final=float(ion[-1][0]), mu_e_actual=float(ion[-1][1]), TARGETMU=float(ion[-1][2]),
                GCE_code=float(ion[-1][3]), cap_code=float(ion[-1][4]), n_cp_rounds=len(re.findall(r"CPM-scf:", txt)),
                mu_e_neutral_first_round=float(scf1[0][1]) if scf1 else None)


def selective_flags(poscar, natoms):
    lines = open(poscar).read().splitlines()
    assert lines[7].strip().lower().startswith("selective") and lines[8].strip().lower().startswith("cartesian")
    fl = []
    for l in lines[9:9 + natoms]:
        p = l.split()
        fl.append([x == "F" for x in p[3:6]])
    return np.array(fl)


reg = json.load(open(f"{ROOT}/03_pilot/run_registry.json"))
states, excluded, geoms = {}, {}, {}
for name, r in sorted(reg.items()):
    d = f"{ROOT}/{r['dir']}"
    reason = None
    if r.get("last_state") != "COMPLETED":
        reason = f"SLURM state {r.get('last_state')}"
    elif r.get("cpm_ion_lines", 0) == 0:
        reason = "no CP (CPM-ion) result: neutral/convergence test, not a constant-potential state"
    elif r["dir"].startswith("04_speedtest"):
        reason = "speed test"
    if reason:
        excluded[name] = dict(dir=r["dir"], reason=reason); continue
    atoms = read(f"{d}/POSCAR")
    n = len(atoms)
    key = (tuple(np.round(atoms.get_cell().flatten(), 4)), tuple(map(tuple, np.round(atoms.get_positions(), 4))))
    prefix = re.split(r"_(dU|muref)", name)[0].replace("_fastcfg", "")
    gid = geoms.setdefault(key, prefix)
    assert gid == prefix, f"geometry hash of {name} matches a different prefix {gid}"
    tags = incar_tags(f"{d}/INCAR")
    oc = parse_outcar(f"{d}/OUTCAR", n)
    cpm = parse_cpm(f"{d}/log.out")
    fixed = selective_flags(f"{d}/POSCAR", n)
    N0 = ZVAL_AU * n
    dN = cpm["N_e_final"] - N0
    gce_check = cpm["GCE_code"] - (oc["E_sigma0"] - cpm["mu_e_actual"] * dN)
    group = "v0_main_production" if r["config_version"].startswith("production") else "reference_accurate"
    sol_z0, sol_z1 = float(tags["SOL_Z0"]), float(tags["SOL_Z1"])
    fields = {}
    for f in FIELD_FILES:
        p = f"{d}/{f}"
        if os.path.exists(p) and os.path.getsize(p) > 0:
            fields[f] = dict(path=p, bytes=os.path.getsize(p), grid=grid_dims(p), note=FIELD_SEMANTICS[f])
    states[name] = dict(
        state_id=name, geometry_id=gid, group=group, source_dir=d, config_version=r["config_version"],
        structure=dict(n_atoms=n, species=["Au"] * n, cell_A=atoms.get_cell().tolist(), positions_A=atoms.get_positions().tolist(),
                       pbc=[True, True, True], fixed_xyz=fixed.tolist(), n_fixed_atoms=int(fixed.all(axis=1).sum()),
                       metal_top_z_A=float(atoms.get_positions()[:, 2].max()), metal_bottom_z_A=float(atoms.get_positions()[:, 2].min())),
        boundary=dict(Lz_A=float(atoms.get_cell()[2][2]), single_sided_solvent=tags.get("LVAC", ".FALSE.").upper().startswith(".T"),
                      SOL_Z0_A=sol_z0, SOL_Z1_A=sol_z1, SOL_SIGMA_A=float(tags.get("SOL_SIGMA", 0.8)),
                      D_STERN_A_default=2.0, ION_Z0_A_derived=sol_z0 + 2.0, ION_Z1_A_derived=sol_z1 - 2.0,
                      DIPOL_frac=[float(x) for x in tags.get("DIPOL", "0.5 0.5 0.5").split()], IDIPOL=int(tags.get("IDIPOL", 0)),
                      ISOL=int(tags.get("ISOL", 0)), C_MOLAR=float(tags.get("C_MOLAR", 0)), R_ION_A=float(tags.get("R_ION", 0)),
                      T_K_code_default=298.0),
        numerics=dict(PREC=tags.get("PREC"), ENCUT_eV=float(tags.get("ENCUT")), ALGO=tags.get("ALGO"), EDIFF=tags.get("EDIFF"),
                      ISMEAR=int(tags.get("ISMEAR")), SIGMA_eV=float(tags.get("SIGMA")), kpoints=r.get("kpoints"),
                      FERMICONVERGE_eV=float(tags.get("FERMICONVERGE")), NESCHEME=int(tags.get("NESCHEME", 0)),
                      functional="PBE (POTCAR PAW_PBE Au 04Oct2007, ZVAL=11)", parallel=r.get("parallel"), slurm_job_ids=r.get("job_ids")),
        electronic_state=dict(N_neutral=N0, N_e_final=cpm["N_e_final"], delta_N_e=dN, net_electronic_charge_e=-dN,
                              induced_sigma_e_per_A2=-dN / (atoms.get_volume() / atoms.get_cell()[2][2]),
                              mu_e_actual_eV=cpm["mu_e_actual"], TARGETMU_eV=cpm["TARGETMU"], mu_e_minus_target_eV=cpm["mu_e_actual"] - cpm["TARGETMU"],
                              FERMICONVERGE_eV=float(tags.get("FERMICONVERGE")), mu_e_neutral_first_round_eV=cpm["mu_e_neutral_first_round"],
                              E_fermi_OUTCAR_eV=oc["E_fermi_OUTCAR"], cap_code_e_per_V=cpm["cap_code"], n_cp_rounds=cpm["n_cp_rounds"],
                              note="mu_e_actual is the converged value (CPM-ion line); TARGETMU is the request. Label states by mu_e_actual/N_e_final."),
        labels=dict(E_free_TOTEN_eV=oc["E_free_TOTEN"], E_without_entropy_eV=oc["E_without_entropy"], E_sigma0_eV=oc["E_sigma0"],
                    GCE_code_eV=cpm["GCE_code"], GCE_code_identity="GCE_code = E_sigma0 - mu_e_actual*(N_e_final - N_neutral)",
                    GCE_identity_residual_eV=gce_check, TOTEN_minus_mu_dN_eV=oc["E_free_TOTEN"] - cpm["mu_e_actual"] * dN,
                    energy_label_status="UNCONFIRMED: which term maps to the trainer's energy target (and its consistency with the forces, which are derivatives of the free energy TOTEN at fixed N_e) is not decided here; select explicitly at training time.",
                    forces_eV_per_A=oc["forces"].tolist(), forces_note="raw forces on all atoms from the single TOTAL-FORCE block (static run); forces on fixed atoms are NOT zeroed; identical to vasprun.xml",
                    total_drift_eV_per_A=oc["total_drift"], stress_kB_as_printed=oc["stress_kB_as_printed"],
                    stress_note="printed by VASP for the periodic cell including vacuum/solvent; not a physical slab stress",
                    max_force_eV_per_A=float(np.abs(oc["forces"]).max()),
                    max_force_fixed_atoms_eV_per_A=float(np.abs(oc["forces"][fixed.all(axis=1)]).max()) if fixed.all(axis=1).any() else None),
        fields=fields,
    )

# configuration duplicates: same geometry + same TARGETMU in both groups
by_key = {}
for s in states.values():
    by_key.setdefault((s["geometry_id"], round(s["electronic_state"]["TARGETMU_eV"], 4)), []).append(s["state_id"])
for k, ids in by_key.items():
    if len(ids) > 1:
        for i in ids:
            states[i]["configuration_duplicate_of"] = [j for j in ids if j != i]
            states[i]["duplicate_note"] = "same geometry and TARGETMU computed in another configuration: use as a configuration comparison, not as additional structural coverage"

os.makedirs(OUT, exist_ok=True)
json.dump(dict(created="2026-09-27", label_standard_v0="production configuration (PREC=Normal, ALGO=Fast, NPAR=16 hybrid, FERMICONVERGE=0.01)",
               states=states, excluded=excluded), open(f"{OUT}/states.json", "w"), indent=1)
json.dump({sid: s["fields"] for sid, s in states.items()}, open(f"{OUT}/fields_index.json", "w"), indent=1)

# extxyz per group -- keys deliberately explicit; no generic 'energy'/'forces' keys
for group in ("v0_main_production", "reference_accurate"):
    lst = []
    for sid, s in states.items():
        if s["group"] != group:
            continue
        at = Atoms("Au%d" % s["structure"]["n_atoms"], positions=s["structure"]["positions_A"], cell=s["structure"]["cell_A"], pbc=True)
        es, lb = s["electronic_state"], s["labels"]
        at.info.update(dict(state_id=sid, geometry_id=s["geometry_id"], group=group, config_version=s["config_version"],
                            N_neutral=es["N_neutral"], N_e_final=es["N_e_final"], delta_N_e=es["delta_N_e"],
                            mu_e_actual_eV=es["mu_e_actual_eV"], TARGETMU_eV=es["TARGETMU_eV"], mu_e_neutral_eV=es["mu_e_neutral_first_round_eV"],
                            E_free_TOTEN_eV=lb["E_free_TOTEN_eV"], E_without_entropy_eV=lb["E_without_entropy_eV"],
                            E_sigma0_eV=lb["E_sigma0_eV"], GCE_code_eV=lb["GCE_code_eV"], TOTEN_minus_mu_dN_eV=lb["TOTEN_minus_mu_dN_eV"],
                            energy_label_status="UNCONFIRMED", SOL_Z0_A=s["boundary"]["SOL_Z0_A"], SOL_Z1_A=s["boundary"]["SOL_Z1_A"],
                            PREC=s["numerics"]["PREC"], kpoints=s["numerics"]["kpoints"].replace(" ", "x"),
                            configuration_duplicate_of=",".join(s.get("configuration_duplicate_of", [])) or "none"))
        at.arrays["forces_raw"] = np.array(lb["forces_eV_per_A"])
        at.arrays["fixed"] = np.array(s["structure"]["fixed_xyz"]).all(axis=1)
        lst.append(at)
    write(f"{OUT}/{group}.extxyz", lst, format="extxyz")
    print(f"{group}: {len(lst)} states -> {OUT}/{group}.extxyz")

# ---- trainer-format files (FermiMACE fork, env cp-mace-0524) ----
# The fork reads: positions/species; forces via --forces_key; energy via --energy_key; and the
# HARD-CODED info keys  potential = mu_e (target, eV)  and  electron = N_e (model INPUT, must always
# be present -- a missing key silently becomes 0 electrons).  Energy candidates are written under
# their explicit names only; the trainer is told which one to use with --energy_key at run time.
# Group convention found in the existing converters (surface_charge/6-Au/*/vasp2all.py):
# energy = 'energy without entropy', potential = CPM-ion mu_e, electron = CPM-ion N_ele.
os.makedirs(f"{OUT}/mace_input", exist_ok=True)
for group in ("v0_main_production", "reference_accurate"):
    lst = []
    for sid, s in states.items():
        if s["group"] != group:
            continue
        at = Atoms("Au%d" % s["structure"]["n_atoms"], positions=s["structure"]["positions_A"], cell=s["structure"]["cell_A"], pbc=True)
        es, lb = s["electronic_state"], s["labels"]
        at.info.update(dict(config_type="Default", state_id=sid, geometry_id=s["geometry_id"],
                            potential=es["mu_e_actual_eV"], electron=es["N_e_final"],
                            E_free_TOTEN_eV=lb["E_free_TOTEN_eV"], E_without_entropy_eV=lb["E_without_entropy_eV"],
                            E_sigma0_eV=lb["E_sigma0_eV"], GCE_code_eV=lb["GCE_code_eV"], TOTEN_minus_mu_dN_eV=lb["TOTEN_minus_mu_dN_eV"]))
        at.arrays["forces"] = np.array(lb["forces_eV_per_A"])
        lst.append(at)
    write(f"{OUT}/mace_input/{group}.mace.extxyz", lst, format="extxyz")
    print(f"trainer format: {group}: {len(lst)} configs -> {OUT}/mace_input/{group}.mace.extxyz")
open(f"{OUT}/mace_input/KEYS.md", "w").write("""# Keys in *.mace.extxyz (FermiMACE fork, env cp-mace-0524)

- `forces` (arrays, eV/A): pass `--forces_key=forces`. Raw forces on all atoms (fixed atoms not zeroed).
- `potential` (info, eV): actual converged mu_e (CPM-ion line). Literal key required by the fork (data/utils.py:145).
  It is a TARGET (the fork predicts mu_e from structure + N_e).
- `electron` (info): actual converged N_e (CPM-ion line). Literal key required (data/utils.py:146). It is a model INPUT.
- energy candidates (info, eV), choose ONE with `--energy_key=<name>`:
  `E_free_TOTEN_eV` (free energy; the forces are its derivatives at fixed N_e), `E_without_entropy_eV`
  (the group's existing convention in surface_charge/6-Au converters), `E_sigma0_eV`, `GCE_code_eV`
  (= E_sigma0 - mu_e*(N_e - N0)), `TOTEN_minus_mu_dN_eV` (= F - mu_e*(N_e - N0)).
  No generic `energy` key is written on purpose: the mapping to the trainer's energy target is UNCONFIRMED.
- `config_type=Default` (for --config_type_weights), `state_id`, `geometry_id`.
""")

# selection list
lines = ["# Dataset v0 — selection (auto-generated by scripts/export_dataset_v0.py, 2026-09-27)", "",
         "Label standard for v0: the production configuration. The PREC=Accurate pilot runs form a separate reference group; "
         "their energies are not to be merged with the production energies without a verified offset. Same geometry + same TARGETMU in "
         "both groups is a configuration duplicate (comparison only). Same geometry at different actual charge = separate states, one geometry_id.", "",
         "| state | group | geometry | N_e final (neutral) | μ_e actual (target) | TOTEN (eV) | max|F| (eV/Å) | duplicate of |", "|---|---|---|---|---|---|---|---|"]
for sid, s in sorted(states.items(), key=lambda kv: (kv[1]["group"], kv[1]["geometry_id"], kv[1]["electronic_state"]["TARGETMU_eV"])):
    es, lb = s["electronic_state"], s["labels"]
    lines.append(f"| {sid} | {s['group']} | {s['geometry_id']} | {es['N_e_final']:.6f} ({es['N_neutral']:.0f}) | {es['mu_e_actual_eV']:.6f} ({es['TARGETMU_eV']:.4f}) | "
                 f"{lb['E_free_TOTEN_eV']:.5f} | {lb['max_force_eV_per_A']:.3f} | {', '.join(s.get('configuration_duplicate_of', [])) or '—'} |")
lines += ["", "## Excluded", "", "| run | reason |", "|---|---|"] + [f"| {k} | {v['reason']} |" for k, v in sorted(excluded.items())]
gcount = {}
for s in states.values():
    gcount.setdefault(s["group"], set()).add(s["geometry_id"])
lines += ["", "## Coverage", ""] + [f"- {g}: {sum(1 for s in states.values() if s['group']==g)} states, {len(v)} distinct geometries ({', '.join(sorted(v))})" for g, v in gcount.items()]
open(f"{OUT}/selection.md", "w").write("\n".join(lines) + "\n")
print("\n".join(lines[-3:]))
print(f"excluded: {len(excluded)}")
