#!/usr/bin/env python3
"""Pick candidate CENTRES on the evolved parent surfaces: new local environments, not random crops.

Per the plan (section 5): the descriptor is per SURFACE ATOM, not a slab average, and selection is farthest-
point in that space under a family cover, with three guards: (a) novelty against the old dataset is rewarded
but isolated anomalies are not taken (a centre must have at least one near neighbour among the candidates),
(b) a larger-neighbourhood context (step density, local relief, under-coordinated count within 10 A) is
carried along so the chosen set is not all one kind of site, (c) nothing is selected by any expected response.

Grouping: within each class the 8 parents are split 6 train / 1 validation / 1 test by seed, and every frame,
crop and minimised version of one parent stays in that parent's group. The split is written here, once.

SOAP (dscribe), fixed for selection only, not model hyperparameters:
    species Au, periodic True with pbc (T, T, F), r_cut 6 A, n_max 8, l_max 6, sigma 0.3 A, average off.

Input is the MD trajectory of each parent (traj.lammpstrj) when it exists; otherwise the parent's initial
structure is used as a single frame so the pipeline can be exercised before the MD has run. The output is a
list of 800-1200 centres for the extraction step to try; the final 200 are chosen only after reconstruction,
when diversity is re-checked on what actually survived repair.

Usage (from Au_Cl/):  env PYTHONPATH=rough_sampling_v1/env/pylib:.pyshim <python> rough_sampling_v1/select_centres.py
                      [--n 1000] [--frames-per-parent 8]
Output: rough_sampling_v1/centres/centres.jsonl, centres_summary.md, parent_split.json, old_reference_soap.npz
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
sys.path.insert(0, f"{ROOT}/scripts")
from analysis_spatial import CN_CUT  # noqa: E402   (the same coordination cutoff as the region analysis)

SOAP_KW = dict(species=["Au"], periodic=True, r_cut=6.0, n_max=8, l_max=6, sigma=0.3, average="off")
CONTEXT_R = 10.0        # A; the larger neighbourhood used for the context descriptor
A0 = 4.158
D111 = A0 / np.sqrt(3)
OUT = "centres"         # output folder under rough_sampling_v1 (overridden by --out for test runs)


def coordination(at, rcut=CN_CUT):
    """Neighbour count within rcut, in-plane minimum image, in chunks (4600-atom frames would otherwise need
    a 500 MB temporary per image shift)."""
    P = at.get_positions(); cell = at.get_cell().array; n = len(P); cn = np.zeros(n, int)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = si * cell[0] + sj * cell[1]
            for a in range(0, n, 512):
                d = np.linalg.norm(P[a:a + 512, None, :] - (P[None, :, :] + sh), axis=-1)
                cn[a:a + 512] += ((d < rcut) & (d > 0.1)).sum(1)
    return cn


def surface_atoms(at):
    """Un-buried atoms of the upper surface: fewer than three higher neighbours within 2.35 A laterally, the
    same rule the gallery and the region analysis use."""
    pos = at.get_positions(); cell = at.get_cell().array
    zt = pos[:, 2].max(); cand = np.flatnonzero(pos[:, 2] > zt - 3 * D111)     # only the top three levels can be exposed
    n_above = np.zeros(len(cand), int)
    P = pos[cand]
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = (si * cell[0] + sj * cell[1])[:2]
            for a in range(0, len(P), 512):
                d = np.linalg.norm(P[a:a + 512, None, :2] - (P[None, :, :2] + sh), axis=-1)
                higher = (P[None, :, 2] - P[a:a + 512, None, 2]) > 0.5
                n_above[a:a + 512] += ((d < 2.35) & higher).sum(1)
    return cand[n_above < 3]


def context(at, idx, surf):
    """Per-centre larger-neighbourhood descriptor: relief, under-coordinated count and level count within
    CONTEXT_R in-plane, over surface atoms only."""
    pos = at.get_positions(); cell = at.get_cell().array
    cn = coordination(at)
    P = pos[surf]
    out = []
    for i in idx:
        best = np.full(len(P), np.inf)
        for si in (-1, 0, 1):
            for sj in (-1, 0, 1):
                sh = (si * cell[0] + sj * cell[1])[:2]
                best = np.minimum(best, np.linalg.norm(P[:, :2] + sh - pos[i, :2], axis=1))
        m = best < CONTEXT_R
        z = P[m, 2]
        out.append(dict(cn=int(cn[i]), n_surface_in_10A=int(m.sum()), relief_10A=float(np.ptp(z)) if m.any() else 0.0,
                        n_levels_10A=int(len(np.unique(np.round(z / D111)))),
                        n_under_10A=int((cn[surf[m]] <= 8).sum()), step_density_10A=float((cn[surf[m]] <= 8).mean()) if m.any() else 0.0))
    return out


def soap_of(at, idx):
    from dscribe.descriptors import SOAP
    a = at.copy(); a.set_pbc((True, True, False))
    return SOAP(**SOAP_KW).create(a, centers=list(idx))


# Frames taken from each 200 ps trajectory (one frame per ps; frame k = k ps). Only the 300 K segments are sampled:
# two frames late in the initial 300 K hold and six spread over the final 300 K hold (160-200 ps, after the
# 600 K excursion and the ramp down). The 600 K segment is the exploration that produces the morphology; its
# hot frames are not sampled, so thermal disorder stays at the 300 K level in everything that reaches DFT.
FRAMES_300K = [10, 20, 165, 172, 179, 186, 193, 199]


def frames_for(pid, n_frames, frames=None):
    d = f"{R}/md/{pid}"
    tr = f"{d}/traj.lammpstrj"
    if os.path.exists(tr):
        fr = read(tr, index=":", format="lammps-dump-text")
        # restore species and periodicity
        for f in fr: f.set_chemical_symbols(["Au"] * len(f)); f.set_pbc((True, True, False))
        want = frames if frames is not None else FRAMES_300K
        pick = [k for k in want if k < len(fr)][:n_frames]
        if len(pick) < min(n_frames, len(fr)) // 2:         # short / incomplete trajectory: fall back to an even spread
            pick = np.linspace(0, len(fr) - 1, min(n_frames, len(fr))).astype(int).tolist()
        return [(int(k), fr[k]) for k in pick], "md"
    return [(0, read(f"{R}/parents/{pid}.extxyz"))], "initial"


def old_reference():
    """SOAP vectors of the surface atoms of the 52 selected old geometries: the novelty baseline."""
    p = f"{R}/{OUT}/old_reference_soap.npz"
    if os.path.exists(p): return np.load(p)["X"]
    import csv
    S = json.load(open(f"{ROOT}/dataset_v1/states.json"))["states"]
    by_state = {s["state_id"]: s for s in S.values()}
    gids = sorted({r["geometry_id"] for r in csv.DictReader(open(f"{R}/old_selected_140.csv"))})
    X = []
    for gid in gids:
        s = next(v for v in by_state.values() if v["geometry_id"] == gid)
        at = read(f"{s['source_dir']}/{s['geometry_file']}")
        surf = surface_atoms(at)
        X.append(soap_of(at, surf))
    X = np.vstack(X); os.makedirs(os.path.dirname(p), exist_ok=True); np.savez_compressed(p, X=X); return X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000); ap.add_argument("--frames-per-parent", type=int, default=8)
    ap.add_argument("--out", default="centres", help="output folder under rough_sampling_v1")
    ap.add_argument("--frames", type=int, nargs="*", default=None, help=f"frame indices to sample (default {FRAMES_300K})")
    a = ap.parse_args()
    global OUT; OUT = a.out
    frames_used = {}
    os.makedirs(f"{R}/{OUT}", exist_ok=True)
    M = json.load(open(f"{R}/parents/parents_manifest.json"))
    # split by seed within each class: 6 train / 1 val / 1 test, fixed once
    split = {}
    for cls in "ABCD":
        ps = sorted([m["parent_id"] for m in M if m["cls"] == cls], key=lambda x: (x.split("_s")[0], int(x.split("_s")[1])))
        for k, pid in enumerate(ps): split[pid] = "test" if k == len(ps) - 1 else "val" if k == len(ps) - 2 else "train"
    json.dump(split, open(f"{R}/{OUT}/parent_split.json", "w"), indent=1)

    Xref = old_reference(); nref = np.linalg.norm(Xref, axis=1, keepdims=True); Uref = Xref / np.maximum(nref, 1e-12)
    cands, X = [], []
    src_kind = collections.Counter()
    for m in M:
        pid = m["parent_id"]
        frames, kind = frames_for(pid, a.frames_per_parent, a.frames); src_kind[kind] += 1
        frames_used[pid] = [fi for fi, _ in frames]
        for fi, at in frames:
            surf = surface_atoms(at)
            if len(surf) == 0: continue
            # subsample surface atoms per frame so one big frame cannot flood the pool; the seed is the
            # parent's build seed and the frame index, so a re-run gives the same pool (str hashes would not)
            rng = np.random.default_rng([int(m["seed"]), int(fi), 7919])
            idx = rng.choice(surf, size=min(120, len(surf)), replace=False)
            Xs = soap_of(at, idx); ctx = context(at, idx, surf)
            U = Xs / np.maximum(np.linalg.norm(Xs, axis=1, keepdims=True), 1e-12)
            novelty = 1.0 - (U @ Uref.T).max(axis=1)          # 1 - max cosine similarity to any old surface atom
            for k, i in enumerate(idx):
                cands.append(dict(parent_id=pid, cls=m["cls"], split=split[pid], frame=fi, source=kind, atom=int(i),
                                  novelty_vs_old=float(novelty[k]), **ctx[k]))
                X.append(Xs[k])
    X = np.vstack(X); U = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    # guard (a): no isolated anomalies -- each candidate needs a neighbour within the pool at cosine > 0.98
    sim_nn = np.array([np.sort(U[i] @ U.T)[-2] for i in range(len(U))])
    ok = sim_nn > 0.98
    # guard (c)-adjacent sanity: coordination of a centre must be a surface value
    ok &= np.array([5 <= c["cn"] <= 11 for c in cands])
    pool = np.flatnonzero(ok)
    # farthest-point under a class x split cover: quota 1/4 per class, within class 6:1:1 by parent count
    quota = {}
    for cls in "ABCD":
        for sp, frac in (("train", 6 / 8), ("val", 1 / 8), ("test", 1 / 8)):
            quota[(cls, sp)] = int(round(a.n / 4 * frac))
    chosen = []
    for key, q in quota.items():
        sub = [i for i in pool if (cands[i]["cls"], cands[i]["split"]) == key]
        if not sub: continue
        # seed with the most novel, then farthest-point in SOAP space (cosine distance)
        order = [max(sub, key=lambda i: cands[i]["novelty_vs_old"])]
        dmin = 1.0 - U[sub] @ U[order[0]]
        while len(order) < min(q, len(sub)):
            j = sub[int(np.argmax(dmin))]
            order.append(j); dmin = np.minimum(dmin, 1.0 - U[sub] @ U[j])
        chosen += order
    with open(f"{R}/{OUT}/centres.jsonl", "w") as f:
        for i in chosen: f.write(json.dumps(dict(cands[i], candidate_index=int(i))) + "\n")
    # summary
    by = collections.Counter((cands[i]["cls"], cands[i]["split"]) for i in chosen)
    nov = np.array([cands[i]["novelty_vs_old"] for i in chosen])
    json.dump(dict(frames_used=frames_used, rule="300 K segments only; see FRAMES_300K"), open(f"{R}/{OUT}/frames_used.json", "w"), indent=1)
    L = [f"# Candidate centres", "", f"{len(chosen)} centres from a pool of {len(pool)} (of {len(cands)} surface atoms sampled; "
         f"{int((~ok).sum())} removed as isolated or non-surface). Sources: {dict(src_kind)}; frames per parent: "
         f"{sorted({tuple(v) for v in frames_used.values()})[0] if frames_used else '-'}.", "",
         "| class | train | val | test |", "|---|---|---|---|"]
    for cls in "ABCD": L.append(f"| {cls} | {by[(cls,'train')]} | {by[(cls,'val')]} | {by[(cls,'test')]} |")
    L += ["", f"Novelty vs the 52 old geometries (1 - max cosine): median {np.median(nov):.4f}, "
          f"10th-90th pct {np.percentile(nov,10):.4f}-{np.percentile(nov,90):.4f}.",
          "Context (within 10 A): " + ", ".join(f"{k} median {np.median([cands[i][k] for i in chosen]):.2f}"
                                               for k in ("cn", "relief_10A", "n_under_10A", "n_levels_10A")), "",
          "These are candidates for extraction; the final 200 are fixed only after repair, when diversity is re-checked."]
    open(f"{R}/{OUT}/centres_summary.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
