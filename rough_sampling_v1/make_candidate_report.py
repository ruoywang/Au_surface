#!/usr/bin/env python3
"""Assemble the candidate report for the rough states from the pipeline's own outputs (no new numbers are computed
here except counts): MD outcome, centre pool, cut/repair outcome, the candidate list with its source distribution,
atom-count / k-mesh distribution, rejection reasons, the first-batch dry run and the staged budget.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/make_candidate_report.py [--out rough200] [--dftdir dft]
Output: rough_sampling_v1/CANDIDATE_REPORT.md
"""
import argparse
import collections
import json
import os
import sys

import numpy as np

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
sys.path.insert(0, R)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="rough200"); ap.add_argument("--dftdir", default="dft"); ap.add_argument("--cells", default="cells")
    a = ap.parse_args()
    man = json.load(open(f"{R}/{a.out}/rough200_manifest.json")); S = man["states"]
    md = json.load(open(f"{R}/md/md_summary.json"))
    cells = [json.loads(l) for l in open(f"{R}/{a.cells}/cells_manifest.jsonl")]
    rej = json.load(open(f"{R}/centres/rejects.json")); fr = json.load(open(f"{R}/centres/frames_used.json"))
    q = json.load(open(f"{R}/{a.dftdir}/queue.json")) if os.path.exists(f"{R}/{a.dftdir}/queue.json") else []
    L = [f"# Rough-state CANDIDATE report ({man['created']})", "",
         f"**Status: {man['status']}.** {man['n_states']} candidate states; nothing submitted; the list is not frozen. "
         "Approval = freeze the list AND accept the budget below; then the first batch of 8 is submitted, recalibrated, then the rest.", "",
         "## 1. MD on the 32 parents (candidate generation only; nothing from it is a label)", "",
         "| parent | wall (min) | frames | T at the 8 sampled 300 K frames (K) | layer counts initial -> final (0..4) | atoms that changed layer |", "|---|---|---|---|---|---|"]
    for r in md:
        L.append(f"| {r['parent']} | {r['wall_min']} | {r['n_frames']} | {min(r['T_at_sampled_frames_K'])}-{max(r['T_at_sampled_frames_K'])} | {r['layers_initial'][:5]} -> {r['layers_final'][:5]} | {r['n_changed_layer']} |")
    T = [t for r in md for t in r["T_at_sampled_frames_K"]]
    L += ["", f"All 32 runs completed the 200 ps protocol (DONE line + dump ending at step 100000). Sampled frames: {fr['frame_times_ps']} ps, "
          f"all in 300 K holds, logged T {min(T)}-{max(T)} K. Detached atoms in any sampled frame: {fr['detached_atoms_per_frame'] or 'none'}.", "",
          "## 2. Centre pool", ""]
    L += open(f"{R}/centres/centres_summary.md").read().splitlines()[2:]
    L += ["", "## 3. Cut and repair", ""] + open(f"{R}/{a.cells}/cells_summary.md").read().splitlines()[2:]
    rep = collections.Counter(c["repair"] for c in cells if c["status"] == "PASS")
    L += ["", f"Repair outcomes over all passing cells: {dict(rep)}.", "", "## 4. The candidate list", ""]
    L += open(f"{R}/{a.out}/rough200_summary.md").read().splitlines()[2:]
    if q:
        L += ["", "## 5. DFT inputs, first batch (dry run) and staged budget", "",
              f"{len(q)} input sets written under `{a.dftdir}/<cell_id>/rough__mu<TARGETMU>/` (production standard K.8; NELECT start guess recorded in each INCAR). "
              "Queue status: " + str(dict(collections.Counter(t["status"] for t in q))) + ".", ""]
        import rough_dft
        first = rough_dft.pick_first(q, 8)
        L += ["First batch that WOULD be submitted (`rough_dft.py submit --first 8`; needs a frozen list and --confirm):", "",
              "| task | class | split | cell | atoms | k-mesh | U (V) | TARGETMU (eV) | NELECT guess | walltime |", "|---|---|---|---|---|---|---|---|---|---|"]
        for t in first:
            L.append(f"| {t['task_id']} | {t['cls']} | {t['split']} | {t['cell']} | {t['n_atoms']} | {t['kpoints']} | {t['U_V']:+.2f} | {t['TARGETMU']:.4f} | {t['nelect_guess']:.2f} | {rough_dft.hhmmss(t['walltime_min'])} |")
        L += [""] + open(f"{R}/{a.dftdir}/budget.md").read().splitlines()[2:]
    figs = sorted(f for f in os.listdir(f"{R}/figures") if f.endswith(".png")) if os.path.isdir(f"{R}/figures") else []
    if figs: L += ["", "## 6. Figures", ""] + [f"- `figures/{f}`" for f in figs]
    L += ["", "## 7. Decisions needed", "",
          "1. Approve or amend the candidate list (`rough200/rough200_cells.csv` has one row per state with its source, CN stratum, U and morphology descriptors).",
          "2. Approve the computing commitment (section 5): node-hours by scenario and ~11 GB of fields per state.",
          "3. Then: `finalize_200.py freeze` -> `rough_dft.py submit --first 8 --confirm` -> recalibrate the cost model on those 8 -> remaining 192."]
    open(f"{R}/CANDIDATE_REPORT.md", "w").write("\n".join(L) + "\n"); print(f"wrote {R}/CANDIDATE_REPORT.md ({len(L)} lines)")


if __name__ == "__main__":
    main()
