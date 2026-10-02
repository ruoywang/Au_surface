#!/usr/bin/env python3
"""6x6 trial cuts: can a 6x6 cell carry the target environment of a centre? Geometry only, no DFT.

Sources
  --frozen (default)   every state of the frozen rough list: same parent, MD frame (by timestep), centre atom, split and
                       potential; target_environment = the centre's CN class + its 6 A neighbourhood; protected = that set
  --centres FILE       jsonl candidates with explicit targets (multi-layer parents): parent_id, step, atom, cls,
                       target_environment, protected_atom_ids, plus parents_dir / md_dir of their own test folder

Rule: CELLS = [(6,6), (8,8)], R_CORE = 6.0 unchanged, re-extracted from the ORIGINAL frame (never cut out of a repaired
8x8). 6x6 is preferred when it passes every existing check (core kept to 1e-3 A, neighbourhood 1:1, no seam contact or
seam-created environment, substrate connected, bottom complete, pit <= 1 layer, per-atom periodic mapping); otherwise the
coded reasons are recorded and 8x8 is tried; no escalation to 10x8. "6x6 is small" is never a reason.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh size_complex_test_v1/size_screen.py [--chunk i/N] [--merge]
                                                                                        [--centres FILE --tag complex]
Output: size_complex_test_v1/size_screen.csv (+ .part*.csv), cells/<cell_id>_6x6/{POSCAR,cell.extxyz}, size_screen_summary.md
"""
import argparse
import collections
import csv
import glob
import json
import os
import sys

import numpy as np
from ase.io import read, write

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
T = f"{ROOT}/size_complex_test_v1"
sys.path.insert(0, R)
import extract_cells as ec  # noqa: E402
import trajio  # noqa: E402

SIZES = [(6, 6), (8, 8)]
FIELDS = ["source", "state_id", "cls", "split", "parent_id", "step", "time_ps", "atom", "cn_stratum", "target_environment", "n_protected",
          "U_V", "TARGETMU_eV", "frozen_cell", "frozen_n_atoms", "preferred_cell", "new_n_atoms", "atoms_saved", "six_status", "six_reasons",
          "eight_status", "eight_reasons", "repair", "seam_affected", "kept_8x8_reference", "dft_status", "new_cell_dir"]


def load_frozen():
    man = json.load(open(f"{R}/rough200/rough200_manifest.json")); q = {t["state_id"]: t for t in json.load(open(f"{R}/dft/queue.json"))}
    keep = set(json.load(open(f"{R}/dft/first_batch.json")))
    rows = []
    for s in man["states"]:
        t = q.get(s["state_id"], {})
        rows.append(dict(source="frozen", state_id=s["state_id"], cls=s["cls"], split=s["split"], parent_id=s["parent_id"], step=s["step"], time_ps=s.get("time_ps"),
                         atom=int(s["cell_id"].split("_a")[1].split("_")[0]), cn_stratum=s["cn_stratum"], U_V=s["U_V"], TARGETMU_eV=s["TARGETMU_eV"],
                         frozen_cell=s["cell"], frozen_n_atoms=s["n_atoms"], target_environment=f"centre ({s['cn_stratum']}) + 6 A neighbourhood",
                         protected=None, kept_8x8_reference=t.get("task_id") in keep, dft_status=t.get("status", ""), parents_dir=f"{R}/parents", md_dir=f"{R}/md"))
    return rows


def load_centres(path):
    rows = []
    for l in open(path):
        c = json.loads(l)
        rows.append(dict(source="complex", state_id=c.get("candidate_id") or f"{c['parent_id']}_f{c['step']:06d}_a{c['atom']:04d}", cls=c["cls"], split=c.get("split", ""),
                         parent_id=c["parent_id"], step=c["step"], time_ps=c.get("time_ps"), atom=int(c["atom"]), cn_stratum=c.get("cn_stratum", ""),
                         U_V=c.get("U_V", ""), TARGETMU_eV=c.get("TARGETMU_eV", ""), frozen_cell="", frozen_n_atoms="", target_environment=c["target_environment"],
                         protected=c.get("protected_atom_ids"), kept_8x8_reference=False, dft_status="", parents_dir=c["parents_dir"], md_dir=c["md_dir"]))
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--centres", default=None); ap.add_argument("--tag", default="frozen")
    ap.add_argument("--chunk", default=None); ap.add_argument("--merge", action="store_true"); ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    os.makedirs(f"{T}/cells", exist_ok=True); ec.CELLDIR = f"{T}/cells"
    rows = load_centres(a.centres) if a.centres else load_frozen()
    if a.limit: rows = rows[:a.limit]
    if a.merge:
        out = [r for p in sorted(glob.glob(f"{T}/size_screen.{a.tag}.part*.csv")) for r in csv.DictReader(open(p))]
        order = {r["state_id"]: k for k, r in enumerate(rows)}; out.sort(key=lambda r: order.get(r["state_id"], 10 ** 9))
        with open(f"{T}/size_screen.{a.tag}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader(); [w.writerow(r) for r in out]
        summary(out, a.tag); return
    groups = collections.OrderedDict()
    for r in rows: groups.setdefault((r["parent_id"], r["step"], r["md_dir"]), []).append(r)
    keys = sorted(groups); part = f"{T}/size_screen.{a.tag}.csv"
    if a.chunk:
        i, N = (int(x) for x in a.chunk.split("/")); keys = keys[i::N]; part = f"{T}/size_screen.{a.tag}.part{i}.csv"
    out = []
    for key in keys:
        pid, step, md_dir = key; trajio.MD_DIR = md_dir
        fr = trajio.frame_by_step(pid, step); p0 = read(f"{groups[key][0]['parents_dir']}/{pid}.extxyz"); fr.set_array("fixed", p0.get_array("fixed"))
        for r in groups[key]:
            base = f"{pid}_f{step:06d}_a{r['atom']:04d}"
            passed, log = ec.extract_one(fr, r["atom"], base, repair=True, sizes=SIZES, first_only=True, protected=r["protected"])
            six = [t for t in log if t["cell"] == "6x6"]; eight = [t for t in log if t["cell"] == "8x8"]
            six_codes = sorted({ec.reason_code(p) for t in six for p in t["problems"]}); eight_codes = sorted({ec.reason_code(p) for t in eight for p in t["problems"]})
            pref = passed[0][0] if passed else "none"
            if passed:
                size, sub, info = passed[0]; cid = f"{base}_{size}"; d = f"{T}/cells/{cid}"; os.makedirs(d, exist_ok=True)
                write(f"{d}/POSCAR", sub, format="vasp", direct=False, sort=False); write(f"{d}/cell.extxyz", sub, format="extxyz")
                json.dump(dict(info, trials=log), open(f"{d}/cut_info.json", "w"), indent=1)
            rec = {k: r.get(k, "") for k in FIELDS}
            rec.update(n_protected=len(r["protected"]) if r["protected"] else 0, preferred_cell=pref,
                       new_n_atoms=len(passed[0][1]) if passed else "", atoms_saved=(int(r["frozen_n_atoms"]) - len(passed[0][1])) if (passed and r["frozen_n_atoms"] != "") else "",
                       six_status="PASS" if pref == "6x6" else "FAIL", six_reasons=";".join(six_codes) if pref != "6x6" else "",
                       eight_status="PASS" if pref == "8x8" else ("not tried" if pref == "6x6" else "FAIL"), eight_reasons=";".join(eight_codes) if pref == "none" else "",
                       repair=passed[0][2]["repair"][:60] if passed else "", seam_affected=passed[0][2]["stats"]["n_seam_affected"] if passed else "",
                       new_cell_dir=f"{T}/cells/{base}_{pref}" if passed else "")
            out.append(rec)
        print(f"  {pid}@{step}: {len(groups[key])} centres done", flush=True)
    with open(part, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader(); [w.writerow(r) for r in out]
    if not a.chunk: summary(out, a.tag)


def summary(out, tag):
    c = collections.Counter(r["preferred_cell"] for r in out); reasons = collections.Counter(x for r in out if r["six_status"] == "FAIL" for x in r["six_reasons"].split(";") if x)
    L = [f"# 6x6 screen ({tag}): {len(out)} centres", "", f"preferred cell: {dict(c)}", "",
         "| 6x6 failure reason (coded) | centres |", "|---|---|"] + [f"| {k} | {v} |" for k, v in reasons.most_common()]
    by = collections.defaultdict(collections.Counter)
    for r in out: by[r["cls"]][r["preferred_cell"]] += 1
    L += ["", "| class | 6x6 | 8x8 | none |", "|---|---|---|---|"] + [f"| {k} | {v['6x6']} | {v['8x8']} | {v['none']} |" for k, v in sorted(by.items())]
    sav = [int(r["atoms_saved"]) for r in out if r["atoms_saved"] not in ("", None) and r["preferred_cell"] == "6x6"]
    if sav: L += ["", f"atoms saved per 6x6 cell vs the frozen cell: median {int(np.median(sav))}, range {min(sav)}-{max(sav)}; "
                  f"new 6x6 sizes {min(int(r['new_n_atoms']) for r in out if r['preferred_cell']=='6x6')}-{max(int(r['new_n_atoms']) for r in out if r['preferred_cell']=='6x6')} atoms"]
    open(f"{T}/size_screen_summary.{tag}.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
