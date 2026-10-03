#!/usr/bin/env python3
"""CP-DFT inputs for the size pairs of the test (NOTHING is submitted here).

A pair = one target centre at one potential in two periodic cells (6x6 and 8x8). Rules:
  * same TARGETMU for both members. For S1/S2 the 8x8 is the kept production reference: the new 6x6 targets its ACTUAL
    converged mu_e (read from its log.out when finished; until then the planned TARGETMU is written and the input is
    marked provisional). For M1-M3 both members are new: the 6x6 is meant to run first and the 8x8 is then matched to
    its actual state (same mechanism, --match-from <6x6 dir>).
  * FERMICONVERGE = 0.01 eV, the one fixed standard for every constant-potential run (user rule 2026-10-03; the six pair
    runs of 2026-10-02 were written with 0.001 before that rule — see README); everything else (ENCUT, k-mesh
    rule, solvent window, smearing, EDIFF, parallel layout) is the production standard imported from
    scripts/production.py through rough_sampling_v1/rough_dft.py.
  * vertical box: the production window needs the highest Au below SOL_Z1 - 15 = 19.603 A with the metal bottom at
    5.0 A. If either member of a pair exceeds it, BOTH members get the same extended setting: Lz and SOL_Z1 raised by
    the same delta (rounded up to 0.5 A, 1.0 A margin), recorded as configuration id "z+<delta>". Nothing is compressed.
  * NELECT start guess as in rough_dft (dataset-wide C and PZC); the CP loop decides the final N_e.

Usage (from Au_Cl/):
  rough_sampling_v1/pyrun_rs.sh size_complex_test_v1/pair_inputs.py --plan size_complex_test_v1/pair_plan.json
pair_plan.json: {"pairs": [{"id": "S1", "target": "...", "U_V": -0.42,
                            "six": {"cell_dir": ".../cells/<id>_6x6"}, "eight": {"cell_dir": ".../cells/<id>_8x8" | "ref_task_dir": ".../dft/<cell>/rough__mu<mu>"}}]}
Output: size_complex_test_v1/dft/<pair>_<size>/{POSCAR,INCAR,KPOINTS,POTCAR,job-run}, pair_queue.json, pair_budget.md
"""
import argparse
import json
import math
import os
import re
import shutil
import sys
import time

import numpy as np
from ase.io import read

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"; T = f"{ROOT}/size_complex_test_v1"
sys.path.insert(0, R); sys.path.insert(0, f"{ROOT}/scripts")
import production as P  # noqa: E402
import rough_dft as RD  # noqa: E402

FERMI_PAIR = 0.01   # USER RULE 2026-10-03: every constant-potential run uses FERMICONVERGE = 0.01, no exceptions (the six
                    # pair runs submitted on 2026-10-02 used 0.001; that is history, never to be repeated). write_member refuses anything else.
Z_LIMIT = P.SOL_Z1 - 15.0      # 19.603 A


def converged_mu(task_dir):
    """actual mu_e of a finished reference (last CPM-ion closure), or None."""
    lo = f"{task_dir}/log.out"
    if not os.path.exists(lo): return None
    ion = re.findall(r"CPM-ion:.*?N_ele=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)", open(lo).read())
    return float(ion[-1][1]) if ion and "General timing" in open(f"{task_dir}/OUTCAR").read() else None


def run_ready(cell_dir, dz=0.0):
    at = read(f"{cell_dir}/cell.extxyz"); pos = at.get_positions(); pos[:, 2] += P.ZMIN - pos[:, 2].min()
    cell = at.get_cell().array.copy(); cell[2] = [0, 0, P.LZ + dz]; at.set_cell(cell, scale_atoms=False); at.set_positions(pos); at.set_pbc((True, True, True))
    movable = pos[:, 2] > 8.0; assert 0 < movable.sum() < len(at)
    return at, movable, float(pos[:, 2].max())


def write_member(d, task, at, movable, mu, U, dz, C_uF, U_pzc, fermi):
    os.makedirs(d, exist_ok=True)
    P.write_poscar(f"{d}/POSCAR", at, movable); kp, n = P.kmesh(at); open(f"{d}/KPOINTS", "w").write(kp)
    a, b = at.cell[0][:2], at.cell[1][:2]; A = abs(a[0] * b[1] - a[1] * b[0]); n_neutral = RD.ZVAL_AU * len(at)
    dN = -C_uF * (U - U_pzc) * RD.E_PER_UC_CM2 * A
    incar = P.INCAR_SP.replace("dataset_plan_v1 rev 2 production single point", "size_complex_test_v1 size-pair single point (production standard K.8, FERMICONVERGE 0.001)")
    if abs(float(fermi) - 0.01) > 1e-12: raise SystemExit(f"FERMICONVERGE must be 0.01 for every constant-potential run (user rule 2026-10-03); got {fermi}")
    incar = incar.format(task=task, mu=f"{mu:.4f}", sol_z0=P.SOL_Z0, sol_z1=round(P.SOL_Z1 + dz, 3))
    assert "FERMICONVERGE = 0.01" in incar, "production INCAR template must carry FERMICONVERGE = 0.01"
    incar += (f"\nNELECT = {n_neutral + dN:.4f}   # start guess only: N_neutral {n_neutral:.0f} + {dN:+.4f} e (C = {C_uF:.2f} uF/cm2, U_pzc = {U_pzc:+.4f} V, A = {A:.1f} A^2, U = {U:+.2f} V)\n")
    if dz > 0: incar += f"# vertical box extended by {dz:.1f} A for this size pair (config z+{dz:.1f}): Lz {P.LZ + dz:.3f}, SOL_Z1 {P.SOL_Z1 + dz:.3f}; SOL_Z0 and the 10 A vacuum above SOL_Z1 unchanged\n"
    open(f"{d}/INCAR", "w").write(incar); shutil.copy(P.POTCAR_SRC, f"{d}/POTCAR")
    wall = RD.walltime_minutes(len(at)); open(f"{d}/job-run", "w").write(P.JOBRUN.format(task=task, walltime=RD.hhmmss(wall), partition=RD.ROUGH_PARTITION))
    return dict(n_atoms=len(at), kpoints=f"{n[0]}x{n[1]}x1", nelect_guess=round(n_neutral + dN, 4), walltime_min=wall)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--plan", default=f"{T}/pair_plan.json"); a = ap.parse_args()
    plan = json.load(open(a.plan)); C_uF, U_pzc = RD.charging_constants(); os.makedirs(f"{T}/dft", exist_ok=True)
    q = []; L = ["# Size-pair inputs (nothing submitted)", "", "| pair | member | source | atoms | k | U (V) | TARGETMU (eV) | mu source | config | FERMICONVERGE | walltime | status |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for pr in plan["pairs"]:
        U = pr["U_V"]; mu_plan = round(RD.MU0_REF - U, 4); mu = mu_plan; mu_src = "planned TARGETMU"
        ref = pr["eight"].get("ref_task_dir")
        if ref:
            m_act = converged_mu(ref)
            if m_act is not None: mu, mu_src = round(m_act, 4), f"actual mu_e of the reference {os.path.basename(os.path.dirname(ref))}"
            else: mu_src = "PROVISIONAL: reference not finished; regenerate with its actual mu_e before submitting"
        # vertical box: the same delta for both members
        tops = []
        for key in ("six", "eight"):
            cd = pr[key].get("cell_dir")
            if cd: tops.append(run_ready(cd)[2])
        excess = max(tops) - Z_LIMIT if tops else -1
        dz = math.ceil((excess + 1.0) * 2) / 2 if excess > 0 else 0.0
        for key, size in (("six", "6x6"), ("eight", "8x8")):
            mem = pr[key]
            if mem.get("ref_task_dir"):
                q.append(dict(pair=pr["id"], member=size, task_dir=mem["ref_task_dir"], source="production reference (FERMICONVERGE 0.01)", status="reference"))
                L.append(f"| {pr['id']} | {size} | reference | - | - | {U:+.2f} | {mu_plan:.4f} | its own | production | 0.01 | - | {'finished' if converged_mu(mem['ref_task_dir']) else 'running/queued'} |"); continue
            at, movable, top = run_ready(mem["cell_dir"], dz)
            task = f"{pr['id']}_{size}__mu{mu:.4f}" + (f"__z+{dz:.1f}" if dz else ""); d = f"{T}/dft/{task}"
            info = write_member(d, task, at, movable, mu, U, dz, C_uF, U_pzc, FERMI_PAIR)
            q.append(dict(pair=pr["id"], member=size, task_id=task, dir=d, cell_dir=mem["cell_dir"], U_V=U, TARGETMU=mu, mu_source=mu_src, config=f"z+{dz:.1f}" if dz else "production",
                          top_A=round(top + (0 if not dz else 0), 3), status="pending (not submitted)" if "PROVISIONAL" not in mu_src else "provisional", **info))
            L.append(f"| {pr['id']} | {size} | new | {info['n_atoms']} | {info['kpoints']} | {U:+.2f} | {mu:.4f} | {mu_src[:40]} | {'z+%.1f' % dz if dz else 'production'} | {FERMI_PAIR} | {RD.hhmmss(info['walltime_min'])} | {q[-1]['status']} |")
    json.dump(dict(created=time.strftime("%Y-%m-%d %H:%M"), fermiconverge_pair=FERMI_PAIR, z_limit_A=Z_LIMIT, tasks=q), open(f"{T}/pair_queue.json", "w"), indent=1)
    new = [t for t in q if "dir" in t]
    L += ["", f"{len(new)} new single points; node-hour estimate (rough_dft cost model, seeded median / cold max): "
          f"{sum(RD.COST['seeded_median'] * RD.scale(t['n_atoms']) for t in new) / 60:.0f} / {sum(RD.COST['cold_max'] * RD.scale(t['n_atoms']) for t in new) / 60:.0f} node-h.",
          "All at FERMICONVERGE 0.01 (the fixed standard for every constant-potential run); a pair's second member targets the first member's ACTUAL converged mu_e."]
    open(f"{T}/pair_budget.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
