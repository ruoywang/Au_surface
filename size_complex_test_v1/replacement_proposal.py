#!/usr/bin/env python3
"""Proposal (not a decision): which frozen large cells become 6x6, which stay, and which quota slots go to multi-layer
environments. Potentials are never re-drawn: a resized state keeps its U / TARGETMU; a complex candidate that takes a
slot inherits the released state's U bin (drawn inside the bin with the recorded seed), so the 10 x 20 bin balance holds.

Actions
  kept_reference_8x8      the two production references (running/queued); unchanged
  resize_to_6x6           6x6 passes for this centre; the unsubmitted large cell is replaced by the 6x6 of the same
                          centre, frame, split and potential
  keep_8x8                no valid 6x6 (coded reason); the frozen 8x8 stays
  keep_10x8_no_smaller    neither 6x6 nor 8x8 passes for this centre now; the frozen 10x8 stays (or is dropped if the
                          user prefers -- flagged)
  release_for_complex     a resized 6x6 state proposed to give its slot to a multi-layer candidate: chosen as the most
                          redundant single-layer states (smallest combined SOAP+morphology distance to another chosen
                          state of the same class), class-balanced, as many as complex candidates admitted (cap 24)

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh size_complex_test_v1/replacement_proposal.py [--complex size_complex_test_v1/size_screen.complex.csv] [--cap 24]
Output: size_complex_test_v1/replacement_proposal.csv, replacement_proposal.md
"""
import argparse
import collections
import csv
import json
import os
import sys

import numpy as np

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"; T = f"{ROOT}/size_complex_test_v1"
sys.path.insert(0, R)
from finalize_200 import MORPH_KEYS  # noqa: E402

FIELDS = ["state_id", "action", "cls", "split", "cn_stratum", "U_V", "TARGETMU_eV", "old_cell", "old_n_atoms", "new_cell", "new_n_atoms", "reason", "complex_candidate", "complex_target"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--complex", default=f"{T}/size_screen.complex.csv"); ap.add_argument("--cap", type=int, default=24)
    ap.add_argument("--seed", type=int, default=20261002)
    a = ap.parse_args()
    scr = {r["state_id"]: r for r in csv.DictReader(open(f"{T}/size_screen.frozen.csv"))}
    man = json.load(open(f"{R}/rough200/rough200_manifest.json")); S = {s["state_id"]: s for s in man["states"]}
    rows = []
    for sid, r in scr.items():
        s = S[sid]
        base = dict(state_id=sid, cls=r["cls"], split=r["split"], cn_stratum=r["cn_stratum"], U_V=r["U_V"], TARGETMU_eV=r["TARGETMU_eV"], old_cell=r["frozen_cell"],
                    old_n_atoms=r["frozen_n_atoms"], new_cell="", new_n_atoms="", reason="", complex_candidate="", complex_target="")
        if r["kept_8x8_reference"] == "True": base.update(action="kept_reference_8x8", new_cell=r["frozen_cell"], new_n_atoms=r["frozen_n_atoms"], reason="production reference of the size pairs; also has a valid 6x6 (%s atoms)" % r["new_n_atoms"])
        elif r["preferred_cell"] == "6x6": base.update(action="resize_to_6x6", new_cell="6x6", new_n_atoms=r["new_n_atoms"], reason=r["repair"])
        elif r["preferred_cell"] == "8x8": base.update(action="keep_8x8", new_cell="8x8", new_n_atoms=r["frozen_n_atoms"] if r["frozen_cell"] == "8x8" else r["new_n_atoms"], reason="6x6 failed: " + r["six_reasons"])
        else: base.update(action="keep_10x8_no_smaller", new_cell=r["frozen_cell"], new_n_atoms=r["frozen_n_atoms"], reason="6x6 failed: %s; 8x8 failed: %s (frozen 10x8 stays; could be dropped)" % (r["six_reasons"], r["eight_reasons"]))
        rows.append(base)
    # complex candidates admitted (passed the screen), taking slots from the most redundant resized single-layer states
    admitted = []
    if os.path.exists(a.complex):
        admitted = [c for c in csv.DictReader(open(a.complex)) if c["preferred_cell"] in ("6x6", "8x8")][:a.cap]
    if admitted:
        z = np.load(f"{R}/rough200/rough200_soap.npz", allow_pickle=True); cid = list(z["cell_id"]); X = z["X"]; Mv = z["morph"]
        Mz = (Mv - Mv.mean(0)) / np.maximum(Mv.std(0), 1e-9); U = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
        idx = {c: k for k, c in enumerate(cid)}
        resz = [r for r in rows if r["action"] == "resize_to_6x6" and r["split"] == "train"]
        per_class = collections.Counter(c["cls"][1] if c["cls"].startswith("M") else c["cls"] for c in admitted)   # M1->'1'... distribute evenly over A-D
        n_release = len(admitted); release = []
        for cls in "ABCD":
            cand = [r for r in resz if r["cls"] == cls]; ks = [idx[S[r["state_id"]]["cell_id"]] for r in cand]
            if not ks: continue
            Us, Ms = U[ks], Mz[ks][:, :Mv.shape[1]]; d_soap = 1 - Us @ Us.T; d_mor = np.linalg.norm(Ms[:, None] - Ms[None], axis=-1)
            iu = np.triu_indices(len(ks), 1); D = 0.5 * d_soap / max(np.percentile(d_soap[iu], 90), 1e-9) + 0.5 * d_mor / max(np.percentile(d_mor[iu], 90), 1e-9)
            np.fill_diagonal(D, np.inf); red = D.min(1)
            order = np.argsort(red)[: int(np.ceil(n_release / 4))]
            for o in order: release.append((cand[o], float(red[o])))
        release = sorted(release, key=lambda t: t[1])[:n_release]
        rng = np.random.default_rng(a.seed)
        for (r, red), c in zip(release, admitted):
            lo = float(S[r["state_id"]]["U_bin"][0]); Uv = float(np.round(np.clip(np.round(rng.uniform(lo, lo + 0.1), 2), lo + 0.005, lo + 0.095), 2))
            r.update(action="release_for_complex", new_cell=c["preferred_cell"], new_n_atoms=c["new_n_atoms"], complex_candidate=c["state_id"], complex_target=c["target_environment"],
                     reason=f"most redundant single-layer state of its class (nearest-neighbour distance {red:.3f}); slot and U bin [{lo:+.1f},{lo+0.1:+.1f}) go to the multi-layer candidate (U {Uv:+.2f} V)")
    with open(f"{T}/replacement_proposal.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader(); [w.writerow(r) for r in rows]
    c = collections.Counter(r["action"] for r in rows)
    old_atoms = sum(int(r["old_n_atoms"]) for r in rows); new_atoms = sum(int(r["new_n_atoms"]) for r in rows if r["new_n_atoms"])
    L = [f"# Replacement proposal ({len(rows)} states of the frozen list) -- a PROPOSAL for one-time confirmation", "",
         "| action | states |", "|---|---|"] + [f"| {k} | {v} |" for k, v in c.most_common()]
    L += ["", f"Atoms summed over the 200 states: {old_atoms} (frozen) -> {new_atoms} (proposed); with the cost model ~N^1.5 the single-point cost sum scales by "
          f"{sum(int(r['new_n_atoms'])**1.5 for r in rows if r['new_n_atoms']) / sum(int(r['old_n_atoms'])**1.5 for r in rows):.2f}x (an extrapolation; the pairs will measure it).",
          "", "Potentials: every resized state keeps its U; complex candidates inherit the released U bin. Nothing is re-drawn.",
          "" if admitted else "Multi-layer candidates: none admitted yet (complex screen not run)."]
    open(f"{T}/replacement_proposal.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
