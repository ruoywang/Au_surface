#!/usr/bin/env python3
"""Low-cost behaviour check of the public Au FLARE potential on the structures this project already has.

No DFT is run. The potential is evaluated, through the in-project LAMMPS build, on:
  1. bulk fcc Au: the lattice constant it prefers, against the 4.158 A used throughout this project;
  2. the DFT-relaxed geometries of a flat terrace, a strip step, both adatoms, an island and a pit: the force
     it reports at the DFT minimum (it should be small but will not be zero: different method), then a
     minimisation with the bottom two layers fixed and the RMS displacement that takes;
  3. the fcc-against-hcp adatom energy difference and the step, island and pit formation energies per
     under-coordinated atom, as sanity signs rather than as numbers to compare with anything here.
Everything is recorded; nothing from it enters a training label.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/check_potential.py
Output: rough_sampling_v1/potential/behaviour_check.json, behaviour_check.md
"""
import json
import os
import subprocess
import sys

import numpy as np
from ase.build import bulk
from ase.io import read

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
sys.path.insert(0, R)
from lmpio import write_data, group_lines, read_data, mic_displacement  # noqa: E402
POT = f"{R}/potential/Au_training/lmp_t0.0001_no_bulk_vac_fix3.header2025.flare"
LMP = f"{R}/env/src/lammps-22Jul2025/build_mpi/lmp"
WORK = f"{R}/potential/check_runs"
CASES = ["T-4x4", "Step-8x2", "A1-fcc", "A1-hcp", "Island-7-compact", "Pit-7-compact"]


def run_lmp(d, script):
    open(f"{d}/in.lammps", "w").write(script)
    p = subprocess.run([LMP, "-in", "in.lammps", "-log", "log.lammps"], cwd=d, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"LAMMPS failed in {d}:\n{p.stdout[-2000:]}\n{p.stderr[-2000:]}")
    return open(f"{d}/log.lammps").read()


def thermo(log, key):
    """Last value of a thermo column from a LAMMPS log."""
    lines = log.splitlines(); val = None
    for i, l in enumerate(lines):
        if l.strip().startswith("Step") and key in l.split():
            cols = l.split(); k = cols.index(key)
            for m in lines[i + 1:]:
                t = m.split()
                if not t or not t[0].lstrip("-").isdigit(): break
                val = float(t[k])
    return val


HEAD = f"""units metal
atom_style atomic
boundary p p {{bz}}
newton on
read_data data.lammps
pair_style flare
pair_coeff * * {POT}
neighbor 2.0 bin
{{groups}}
group mobile subtract all frozen
fix hold frozen setforce 0.0 0.0 0.0
variable fm atom sqrt(fx^2+fy^2+fz^2)
compute fmax mobile reduce max v_fm
thermo_style custom step pe c_fmax
thermo 1
"""


def main():
    os.makedirs(WORK, exist_ok=True)
    res = {}
    # 1. bulk lattice constant
    bulkres = []
    for a0 in np.arange(4.00, 4.33, 0.04):
        d = f"{WORK}/bulk_a{a0:.2f}"; os.makedirs(d, exist_ok=True)
        at = bulk("Au", "fcc", a=a0, cubic=True) * (3, 3, 3)
        write_data(at, f"{d}/data.lammps")
        log = run_lmp(d, HEAD.format(bz="p", groups="group frozen empty") + "run 0\n")
        bulkres.append((float(a0), thermo(log, "PotEng") / len(at)))
    aa = np.array([b[0] for b in bulkres]); ee = np.array([b[1] for b in bulkres])
    c = np.polyfit(aa, ee, 3); r = np.roots(np.polyder(c)); a_min = float(min((x.real for x in r if abs(x.imag) < 1e-9), key=lambda x: np.polyval(c, x)))
    res["bulk"] = dict(scan=bulkres, a0_min_A=a_min, E_per_atom_at_min_eV=float(np.polyval(c, a_min)), project_a0_A=4.158)
    # 2/3. the project's relaxed structures
    gal = json.load(open(f"{ROOT}/analysis/gallery/gallery.json"))
    S = json.load(open(f"{ROOT}/dataset_v1/states.json"))["states"]
    for sid in CASES:
        st = [s for s in S.values() if s["structure_id"] == sid and s["config"] == "relax"]
        src = st[0] if st else [s for s in S.values() if s["structure_id"] == sid and s["config"] == "ideal"][0]
        at = read(f"{src['source_dir']}/{src['geometry_file']}")
        z = at.get_positions()[:, 2]; lay = np.round((z - z.min()) / (4.158 / np.sqrt(3)))
        fixed = lay < 2
        d = f"{WORK}/{sid}"; os.makedirs(d, exist_ok=True)
        write_data(at, f"{d}/data.lammps")
        head = HEAD.format(bz="f", groups=group_lines("frozen", fixed))
        log0 = run_lmp(d, head + "run 0\n")
        e0 = thermo(log0, "PotEng"); fmax0 = thermo(log0, "c_fmax")
        log1 = run_lmp(d, head + "min_style cg\nminimize 1e-8 1e-5 2000 20000\nwrite_data min.data\n")
        e1 = thermo(log1, "PotEng")
        a_in = read_data(f"{d}/data.lammps"); atm = read_data(f"{d}/min.data")
        disp = np.linalg.norm(mic_displacement(a_in.get_cell().array, a_in.get_positions(), atm.get_positions()), axis=1)
        res[sid] = dict(n_atoms=len(at), config=src["config"], E_at_dft_geometry_eV=e0, max_force_at_dft_geometry_eV_per_A=fmax0,
                        E_after_min_eV=e1, dE_min_eV=e1 - e0, rms_disp_mobile_A=float(np.sqrt((disp[~fixed] ** 2).mean())),
                        max_disp_mobile_A=float(disp[~fixed].max()))
    # sanity differences
    def E(k): return res[k]["E_after_min_eV"]
    def N(k): return res[k]["n_atoms"]
    epa = res["bulk"]["E_per_atom_at_min_eV"]
    res["signs"] = dict(
        adatom_hcp_minus_fcc_eV=E("A1-hcp") - E("A1-fcc"),
        adatom_binding_vs_bulk_eV=E("A1-fcc") - E("T-4x4") - epa,
        # Step-8x2 is the 16-site cell of T-4x4 plus an 8-atom half-layer strip: the strip's atoms are added at
        # the bulk energy, the flat slab is NOT rescaled (same area, same number of surfaces)
        step_strip_formation_eV=E("Step-8x2") - E("T-4x4") - 8 * epa,
        step_formation_per_edge_atom_eV=(E("Step-8x2") - E("T-4x4") - 8 * epa) / 16,
        island7_formation_eV=E("Island-7-compact") - E("T-4x4") * (N("Island-7-compact") - 7) / N("T-4x4") - 7 * epa,
        pit7_formation_eV=E("Pit-7-compact") - E("T-4x4") * (N("Pit-7-compact") + 7) / N("T-4x4") + 7 * epa,
        note="formation energies against the flat slab scaled by atom count and the bulk energy per atom; sanity signs only")
    json.dump(res, open(f"{R}/potential/behaviour_check.json", "w"), indent=1)
    L = ["# Behaviour check of the Au FLARE potential on this project's structures", "",
         f"Bulk fcc: energy minimum at a0 = {a_min:.3f} A (project uses 4.158 A); E/atom {res['bulk']['E_per_atom_at_min_eV']:.4f} eV.", "",
         "| structure | atoms | geometry | max force at DFT geometry (eV/A) | energy drop on minimising (eV) | RMS / max displacement of mobile atoms (A) |",
         "|---|---|---|---|---|---|"]
    for sid in CASES:
        r = res[sid]
        L.append(f"| {sid} | {r['n_atoms']} | {r['config']} | {r['max_force_at_dft_geometry_eV_per_A']:.3f} | {r['dE_min_eV']:+.3f} | {r['rms_disp_mobile_A']:.3f} / {r['max_disp_mobile_A']:.3f} |")
    s = res["signs"]
    L += ["", f"- adatom hcp minus fcc: {s['adatom_hcp_minus_fcc_eV']:+.3f} eV",
          f"- adatom binding relative to bulk: {s['adatom_binding_vs_bulk_eV']:+.3f} eV",
          f"- 8-atom step strip (two 8-atom edges) formation: {s['step_strip_formation_eV']:+.3f} eV, {s['step_formation_per_edge_atom_eV']:+.3f} eV per edge atom",
          f"- 7-atom island formation: {s['island7_formation_eV']:+.3f} eV; 7-atom pit formation: {s['pit7_formation_eV']:+.3f} eV", "",
          "These are sanity signs for a candidate-generation tool. None of them is compared with the CP-DFT data and none enters a label."]
    open(f"{R}/potential/behaviour_check.md", "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
