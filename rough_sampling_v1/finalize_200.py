#!/usr/bin/env python3
"""Fix the 200 rough states: which reconstructed cells, which split, which potential. Done ONCE.

Inputs : cells/cells_manifest.jsonl (PASS rows from extract_cells.py) and the reconstructed cells themselves.
Outputs: rough200/rough200_manifest.json, rough200_summary.md, rough200_soap.npz

Selection (plan section 5-7):
  * diversity is re-checked on what SURVIVED repair: the SOAP vector (same settings as the centre selection) of
    the centre atom inside the reconstructed cell, not the parent's; farthest-point sampling in that space per
    (class, split) bucket, seeded by the most novel cell against the 52 old geometries;
  * 200 = 160 train / 20 val / 20 test, 50 per class (40 / 5 / 5); every frame, crop and minimised version of a
    parent stays in that parent's split (parent_split.json, written by select_centres.py);
  * 8 size-pair test states: 4 centres cut in BOTH 8x8 and 10x8, one per class where available, both members at
    the SAME potential so the pair isolates the cell-size dependence;
  * one potential U per geometry, 10 bins of 0.1 V over [-0.5, +0.5] V, 20 geometries per bin, assigned by a
    seeded permutation balanced over buckets; U is drawn uniformly inside its bin and rounded to 0.01 V;
    TARGETMU = mu0 - U with mu0 = -4.9071 eV (the project's reference). Nothing is selected by any expected
    response: the geometry is chosen before the potential is drawn, and the potential does not look at the
    geometry.
  * the manifest is written once; re-running refuses to overwrite it unless --force is given, so a restart of
    the DFT stage never re-draws the potentials. The RNG seed is recorded.

The states are labelled later by their converged mu_e, as everywhere in this project; TARGETMU here is the
requested value only.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/finalize_200.py [--cells cells] [--n 200] [--seed 20261002] [--force]
"""
import argparse
import collections
import json
import os
import sys

import numpy as np
from ase.io import read

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
sys.path.insert(0, R)
import select_centres as sc  # noqa: E402

MU0 = -4.9071
U_MIN, U_MAX, N_BINS = -0.5, 0.5, 10
SPLIT_FRAC = {"train": 0.80, "val": 0.10, "test": 0.10}
N_PAIRS = 4


def soap_centre(cell_dir, centre_idx):
    at = read(f"{cell_dir}/cell.extxyz")
    return sc.soap_of(at, [centre_idx])[0]


def farthest_point(U, order_seed, k):
    """indices into U: start at order_seed, then farthest-point by cosine distance."""
    order = [order_seed]; dmin = 1.0 - U @ U[order_seed]
    while len(order) < min(k, len(U)):
        j = int(np.argmax(dmin)); order.append(j); dmin = np.minimum(dmin, 1.0 - U @ U[j])
    return order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", default="cells"); ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261002); ap.add_argument("--force", action="store_true")
    ap.add_argument("--out", default="rough200")
    a = ap.parse_args()
    OUT = f"{R}/{a.out}"; os.makedirs(OUT, exist_ok=True)
    man = f"{OUT}/rough200_manifest.json"
    if os.path.exists(man) and not a.force:
        sys.exit(f"{man} exists: the potentials are fixed. Use --force only to deliberately re-draw (and say so in the report).")
    rows = [json.loads(l) for l in open(f"{R}/{a.cells}/cells_manifest.jsonl")]
    ok = [r for r in rows if r["status"] == "PASS"]
    if not ok: sys.exit("no PASS cells")
    # SOAP of the centre atom in each reconstructed cell + novelty vs the old reference
    sc.OUT = "centres"
    Xref = sc.old_reference(); Uref = Xref / np.maximum(np.linalg.norm(Xref, axis=1, keepdims=True), 1e-12)
    X = np.vstack([soap_centre(f"{R}/{a.cells}/{r['cell_id']}", r["centre_in_cell"]) for r in ok])
    U = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    novelty = 1.0 - (U @ Uref.T).max(axis=1)
    for r, nv in zip(ok, novelty): r["novelty_cell_vs_old"] = float(nv)
    np.savez_compressed(f"{OUT}/rough200_soap.npz", X=X, cell_id=np.array([r["cell_id"] for r in ok]))

    # ---- quotas
    per_class = a.n // 4
    quota = {}
    for cls in "ABCD":
        q = {sp: int(round(per_class * f)) for sp, f in SPLIT_FRAC.items()}
        q["train"] += per_class - sum(q.values())
        for sp, v in q.items(): quota[(cls, sp)] = v
    # ---- size pairs first (test split): centres with both sizes passing
    by_centre = collections.defaultdict(dict)
    for i, r in enumerate(ok): by_centre[r["centre_id"]][r["cell"]] = i
    pair_cands = {cls: [c for c, d in by_centre.items() if len(d) == 2 and ok[d["8x8"]]["split"] == "test" and ok[d["8x8"]]["cls"] == cls] for cls in "ABCD"}
    chosen = []; used_centres = set(); pairs = []
    rng = np.random.default_rng(a.seed)
    classes_for_pairs = [c for c in "ABCD" if pair_cands[c]]
    # one pair per class while possible, then fill from classes with spare candidates
    pick_order = list(classes_for_pairs)
    while len(pairs) < N_PAIRS and pick_order:
        cls = pick_order.pop(0)
        cands = [c for c in pair_cands[cls] if c not in used_centres]
        if not cands: continue
        # most novel pair centre first
        c = max(cands, key=lambda c_: ok[by_centre[c_]["8x8"]]["novelty_cell_vs_old"])
        pairs.append((cls, c)); used_centres.add(c)
        if len(pairs) < N_PAIRS and len(pick_order) == 0:
            pick_order = [k for k in "ABCD" if any(c_ not in used_centres for c_ in pair_cands[k])]
    pair_members = {}
    for k, (cls, c) in enumerate(pairs):
        for size, i in by_centre[c].items():
            chosen.append(i); pair_members[i] = f"pair{k+1}"
        quota[(cls, "test")] -= 2
    # ---- farthest-point fill per bucket (one cell per centre)
    shortfall = {}
    for (cls, sp), q in quota.items():
        if q <= 0: continue
        sub = [i for i, r in enumerate(ok) if r["cls"] == cls and r["split"] == sp and r["centre_id"] not in used_centres and i not in chosen]
        # one cell per centre: prefer 8x8 when both exist
        seen = set(); sub2 = []
        for i in sorted(sub, key=lambda i: (ok[i]["centre_id"], ok[i]["cell"] != "8x8")):
            if ok[i]["centre_id"] in seen: continue
            seen.add(ok[i]["centre_id"]); sub2.append(i)
        if not sub2: shortfall[(cls, sp)] = q; continue
        Us = U[sub2]; seed = int(np.argmax([ok[i]["novelty_cell_vs_old"] for i in sub2]))
        order = farthest_point(Us, seed, q)
        take = [sub2[j] for j in order]
        if len(take) < q: shortfall[(cls, sp)] = q - len(take)
        chosen += take; used_centres.update(ok[i]["centre_id"] for i in take)
    # ---- potentials: 10 bins, balanced over buckets, drawn once
    edges = np.linspace(U_MIN, U_MAX, N_BINS + 1)
    units = []   # a unit = one geometry or one size pair (shares U)
    seen_pairs = set()
    for i in chosen:
        if i in pair_members:
            if pair_members[i] in seen_pairs: continue
            seen_pairs.add(pair_members[i]); units.append([j for j in chosen if pair_members.get(j) == pair_members[i]])
        else: units.append([i])
    # bins assigned round-robin over a shuffled unit order within each bucket, bucket order shuffled too
    bucket_of = lambda u: (ok[u[0]]["cls"], ok[u[0]]["split"])
    buckets = collections.defaultdict(list)
    for u in units: buckets[bucket_of(u)].append(u)
    bin_counts = np.zeros(N_BINS, int); assign = {}
    bkeys = list(buckets); rng.shuffle(bkeys)
    for bk in bkeys:
        us = buckets[bk]; rng.shuffle(us)
        for u in us:
            # least-filled bin, ties broken at random
            cands = np.flatnonzero(bin_counts == bin_counts.min()); b = int(rng.choice(cands)); bin_counts[b] += 1
            Uv = float(np.round(rng.uniform(edges[b], edges[b + 1]), 2))
            Uv = min(max(Uv, edges[b] + 0.005), edges[b + 1] - 0.005) if not (edges[b] <= Uv <= edges[b + 1]) else Uv
            for i in u: assign[i] = (b, Uv)
    # ---- manifest
    states = []
    for i in sorted(chosen, key=lambda i: (ok[i]["split"], ok[i]["cls"], ok[i]["cell_id"])):
        r = ok[i]; b, Uv = assign[i]
        states.append(dict(state_id=f"{r['cell_id']}__U{Uv:+.2f}", cell_id=r["cell_id"], centre_id=r["centre_id"], parent_id=r["parent_id"],
                           cls=r["cls"], split=r["split"], cell=r["cell"], n_atoms=r["n_atoms"], frame=r["frame"], source=r["source"],
                           U_V=Uv, U_bin=[float(edges[b]), float(edges[b + 1])], TARGETMU_eV=round(MU0 - Uv, 4),
                           pair=pair_members.get(i), repair=r["repair"], novelty_cell_vs_old=r["novelty_cell_vs_old"],
                           cell_dir=f"{R}/{a.cells}/{r['cell_id']}"))
    info = dict(created="2026-10-02", seed=a.seed, mu0_eV=MU0, n_states=len(states), bins_V=edges.tolist(), quota={f"{k[0]}/{k[1]}": v for k, v in quota.items()},
                shortfall={f"{k[0]}/{k[1]}": v for k, v in shortfall.items()}, pairs=[dict(id=f"pair{k+1}", cls=c, centre_id=cc) for k, (c, cc) in enumerate(pairs)],
                rule="geometry chosen by farthest-point SOAP of the reconstructed cell's centre per class x split; potential drawn afterwards, once; "
                     "TARGETMU is the request, the state is labelled by its converged mu_e", states=states)
    json.dump(info, open(man, "w"), indent=1)
    # ---- summary
    cnt = collections.Counter((s["cls"], s["split"]) for s in states)
    hist = collections.Counter(tuple(s["U_bin"]) for s in states)
    L = [f"# Rough states fixed: {len(states)} (seed {a.seed})", "", "| class | train | val | test |", "|---|---|---|---|"]
    for cls in "ABCD": L.append(f"| {cls} | {cnt[(cls,'train')]} | {cnt[(cls,'val')]} | {cnt[(cls,'test')]} |")
    L += ["", "Size pairs (test, same U for both members): " + (", ".join(f"{p['id']} {p['cls']} {p['centre_id']}" for p in info["pairs"]) or "none"), "",
          "| U bin (V) | states |", "|---|---|"] + [f"| {lo:+.1f} to {hi:+.1f} | {hist[(lo,hi)]} |" for lo, hi in zip(edges[:-1], edges[1:])]
    L += ["", "Cell sizes: " + str(dict(collections.Counter(s["cell"] for s in states))),
          "Atoms: " + ", ".join(f"{k}: {v}" for k, v in sorted(collections.Counter(s["n_atoms"] for s in states).items())),
          "Repair: " + str(dict(collections.Counter(s["repair"] for s in states)))]
    if shortfall: L += ["", "SHORTFALL (bucket: missing): " + str(info["shortfall"]) + " -- fewer passing cells than the quota; not filled from other buckets."]
    L += ["", f"Novelty of the chosen cells vs the 52 old geometries (1 - max cosine): median {np.median([s['novelty_cell_vs_old'] for s in states]):.4f}."]
    open(f"{OUT}/rough200_summary.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
