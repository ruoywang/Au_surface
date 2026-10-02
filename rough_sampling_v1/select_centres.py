#!/usr/bin/env python3
"""Pick candidate CENTRES on the evolved parent surfaces: new local environments, not random crops.

Frames (production mode): the 300 K frames at FRAME_TIMES_PS of every one of the 32 trajectories, located by
LAMMPS TIMESTEP, checked against the protocol stage and the logged temperature (trajio.frames_at). An incomplete
run, a missing frame or an off-stage temperature stops the script and names the parent; there is NO fallback.
`--test-initial` is the only other mode: it uses the initial (pre-MD) parents and labels everything source=initial.

Legality of a centre is decided by geometry, not by coordination number:
  * the atom is exposed (surface_atoms rule), connected to the substrate through the neighbour graph (no
    detached atom or cluster), and has no contact < 2.5 A;
  * CN is recorded and used for STRATIFIED sampling (cn <= 5 adatom/kink, 6-7 edge, 8-9 terrace, >= 10), so
    low-coordinated environments are sampled deliberately instead of being excluded;
  * a centre with no SOAP neighbour at cosine > 0.98 in the pool is flagged rare, not deleted; the geometry
    checks above decide whether it is an error or a rare environment.

Per-frame subsampling keeps every low-CN exposed atom and a random subset of the rest (seeded by the parent's
build seed and the step). Selection is farthest-point in SOAP space within each (class, split, CN stratum)
bucket, seeded by the most novel candidate against the 52 selected old geometries. Novelty rewards new
environments; it never selects by any expected response. The 10 A context (relief, levels, under-coordinated
count, step density) is carried along for the second-layer stratification in finalize_200.py.

SOAP (dscribe), fixed for selection only: species Au, periodic, pbc (T,T,F), r_cut 6, n_max 8, l_max 6, sigma 0.3.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/select_centres.py [--n 1000] [--out centres]
                      [--test-initial]  [--per-frame 60]
Output: <out>/centres.jsonl, centres_summary.md, parent_split.json, frames_used.json, rejects.json, old_reference_soap.npz
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
sys.path.insert(0, f"{ROOT}/scripts"); sys.path.insert(0, R)
from analysis_spatial import CN_CUT  # noqa: E402   (the same coordination cutoff as the region analysis)
import trajio  # noqa: E402

SOAP_KW = dict(species=["Au"], periodic=True, r_cut=6.0, n_max=8, l_max=6, sigma=0.3, average="off")
CONTEXT_R = 10.0        # A; the larger neighbourhood used for the context descriptor
A0 = 4.158
D111 = A0 / np.sqrt(3)
MIN_CONTACT = 2.5
OUT = "centres"
# 300 K frames only: late in the initial hold, and spread over the final hold after the 600 K excursion
FRAME_TIMES_PS = [10.0, 20.0, 165.0, 172.0, 179.0, 186.0, 193.0, 199.0]
STRATA = [("cn<=5", 0, 5), ("cn6-7", 6, 7), ("cn8-9", 8, 9), ("cn>=10", 10, 99)]


def stratum(cn):
    for name, lo, hi in STRATA:
        if lo <= cn <= hi: return name
    return STRATA[-1][0]


def pair_blocks(P, cell, rcut, fn):
    """Apply fn(a0, d_block) to chunked minimum-image distance blocks (9 in-plane images) without the N^2 temporary."""
    n = len(P)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = si * cell[0] + sj * cell[1]
            for a in range(0, n, 512):
                d = np.linalg.norm(P[a:a + 512, None, :] - (P[None, :, :] + sh), axis=-1)
                if si == 0 and sj == 0: d[np.arange(d.shape[0]), np.arange(a, a + d.shape[0])] = np.inf   # the self term only
                fn(a, d)


def neighbour_graph(at, rcut=CN_CUT):
    """cn per atom, min distance per atom, and the adjacency lists within rcut (in-plane minimum image)."""
    P = at.get_positions(); cell = at.get_cell().array; n = len(P)
    cn = np.zeros(n, int); dmin = np.full(n, np.inf); adj = [[] for _ in range(n)]

    def fn(a, d):
        m = d < rcut
        cn[a:a + d.shape[0]] += m.sum(1); dmin[a:a + d.shape[0]] = np.minimum(dmin[a:a + d.shape[0]], d.min(1))
        for i, j in zip(*np.nonzero(m)): adj[a + i].append(int(j))
    pair_blocks(P, cell, rcut, fn)
    return cn, dmin, adj


def connected_to_substrate(adj, roots):
    """Boolean mask: reachable from any root (the fixed bottom atoms) through the neighbour graph."""
    seen = np.zeros(len(adj), bool); stack = list(roots); seen[list(roots)] = True
    while stack:
        i = stack.pop()
        for j in adj[i]:
            if not seen[j]: seen[j] = True; stack.append(j)
    return seen


def coordination(at, rcut=CN_CUT):
    return neighbour_graph(at, rcut)[0]


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


def context(at, idx, surf, cn):
    """Per-centre larger-neighbourhood descriptor: relief, under-coordinated count and level count within
    CONTEXT_R in-plane, over surface atoms only."""
    pos = at.get_positions(); cell = at.get_cell().array
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
        out.append(dict(cn=int(cn[i]), cn_stratum=stratum(int(cn[i])), n_surface_in_10A=int(m.sum()), relief_10A=float(np.ptp(z)) if m.any() else 0.0,
                        n_levels_10A=int(len(np.unique(np.round(z / D111)))),
                        n_under_10A=int((cn[surf[m]] <= 8).sum()), step_density_10A=float((cn[surf[m]] <= 8).mean()) if m.any() else 0.0))
    return out


def soap_of(at, idx):
    from dscribe.descriptors import SOAP
    a = at.copy(); a.set_pbc((True, True, False))
    return SOAP(**SOAP_KW).create(a, centers=list(idx))


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


def frames_for(pid, test_initial):
    if test_initial:
        return [(0, 0.0, float("nan"), read(f"{R}/parents/{pid}.extxyz"))], "initial"
    return trajio.frames_at(pid, FRAME_TIMES_PS), "md"


def allocate(q, sizes):
    """Split quota q over strata proportionally to sqrt(size), at least 1 for a non-empty stratum, never above size."""
    names = [k for k, v in sizes.items() if v > 0]
    if not names: return {}
    w = np.sqrt(np.array([sizes[k] for k in names], float)); raw = q * w / w.sum()
    alloc = {k: max(1, int(np.floor(r))) if q >= len(names) else 0 for k, r in zip(names, raw)}
    alloc = {k: min(v, sizes[k]) for k, v in alloc.items()}
    # distribute the remainder to the strata with the largest fractional parts that still have room
    order = sorted(names, key=lambda k: -(raw[names.index(k)] - np.floor(raw[names.index(k)])))
    while sum(alloc.values()) < min(q, sum(sizes.values())):
        progressed = False
        for k in order:
            if alloc[k] < sizes[k] and sum(alloc.values()) < q: alloc[k] += 1; progressed = True
        if not progressed: break
    return alloc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000); ap.add_argument("--per-frame", type=int, default=60)
    ap.add_argument("--out", default="centres", help="output folder under rough_sampling_v1")
    ap.add_argument("--test-initial", action="store_true", help="TEST MODE: use the initial parents, not MD frames")
    a = ap.parse_args()
    global OUT; OUT = a.out
    os.makedirs(f"{R}/{OUT}", exist_ok=True)
    M = json.load(open(f"{R}/parents/parents_manifest.json"))
    # split by seed within each class: 6 train / 1 val / 1 test, fixed once
    split = {}
    for cls in "ABCD":
        ps = sorted([m["parent_id"] for m in M if m["cls"] == cls], key=lambda x: (x.split("_s")[0], int(x.split("_s")[1])))
        for k, pid in enumerate(ps): split[pid] = "test" if k == len(ps) - 1 else "val" if k == len(ps) - 2 else "train"
    json.dump(split, open(f"{R}/{OUT}/parent_split.json", "w"), indent=1)

    # production mode: every run must be complete BEFORE anything is sampled
    if not a.test_initial:
        bad = [(m["parent_id"], trajio.run_complete(m["parent_id"])[1]) for m in M if not trajio.run_complete(m["parent_id"])[0]]
        if bad:
            for pid, why in bad: print(f"INCOMPLETE {pid}: {why}")
            sys.exit(f"{len(bad)} of {len(M)} parents incomplete; nothing selected (no fallback in production mode)")

    Xref = old_reference(); Uref = Xref / np.maximum(np.linalg.norm(Xref, axis=1, keepdims=True), 1e-12)
    cands, X = [], []
    frames_used = {}; rejects = collections.Counter(); detached_report = {}
    for m in M:
        pid = m["parent_id"]
        frames, kind = frames_for(pid, a.test_initial)
        frames_used[pid] = [dict(step=s, time_ps=t, T_K=T) for s, t, T, _ in frames]
        for step, t_ps, T, at in frames:
            cn, dmin, adj = neighbour_graph(at)
            fixed = np.flatnonzero(at.get_array("fixed").astype(bool)) if "fixed" in at.arrays else np.flatnonzero(at.get_positions()[:, 2] < at.get_positions()[:, 2].min() + 1.5 * D111)
            conn = connected_to_substrate(adj, fixed)
            n_det = int((~conn).sum())
            if n_det: detached_report[f"{pid}@{step}"] = n_det
            surf = surface_atoms(at)
            legal = conn[surf] & (dmin[surf] >= MIN_CONTACT)
            rejects["detached"] += int((~conn[surf]).sum()); rejects["close_contact"] += int((dmin[surf] < MIN_CONTACT).sum())
            surf = surf[legal]
            if len(surf) == 0: continue
            # subsample: keep every low-CN exposed atom, random subset of the rest
            rng = np.random.default_rng([int(m["seed"]), int(step), 7919])
            low = surf[cn[surf] <= 5]; rest = surf[cn[surf] > 5]
            k_rest = max(0, a.per_frame - len(low))
            idx = np.concatenate([low, rng.choice(rest, size=min(k_rest, len(rest)), replace=False) if k_rest and len(rest) else np.array([], int)]).astype(int)
            Xs = soap_of(at, idx); ctx = context(at, idx, surf, cn)
            U = Xs / np.maximum(np.linalg.norm(Xs, axis=1, keepdims=True), 1e-12)
            novelty = 1.0 - (U @ Uref.T).max(axis=1)
            for k, i in enumerate(idx):
                cands.append(dict(parent_id=pid, cls=m["cls"], split=split[pid], step=int(step), time_ps=t_ps, T_K=T, source=kind, atom=int(i),
                                  novelty_vs_old=float(novelty[k]), min_dist_A=float(dmin[i]), **ctx[k]))
                X.append(Xs[k])
    X = np.vstack(X); U = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    # rare-environment flag: no other candidate at cosine > 0.98 (computed in chunks); flagged, never deleted
    sim_nn = np.empty(len(U))
    for a0 in range(0, len(U), 2048):
        S = U[a0:a0 + 2048] @ U.T; S[np.arange(S.shape[0]), np.arange(a0, a0 + S.shape[0])] = -1; sim_nn[a0:a0 + 2048] = S.max(1)
    for c, s in zip(cands, sim_nn): c["rare_flag"] = bool(s <= 0.98); c["nn_cosine_in_pool"] = float(s)
    pool = np.arange(len(cands))
    # quotas: 1/4 per class, 6:1:1 by split, then CN strata proportional to sqrt(size)
    chosen = []
    for cls in "ABCD":
        for sp, frac in (("train", 6 / 8), ("val", 1 / 8), ("test", 1 / 8)):
            q = int(round(a.n / 4 * frac))
            sub = [i for i in pool if (cands[i]["cls"], cands[i]["split"]) == (cls, sp)]
            sizes = {s[0]: sum(cands[i]["cn_stratum"] == s[0] for i in sub) for s in STRATA}
            for st, qs in allocate(q, sizes).items():
                ss = [i for i in sub if cands[i]["cn_stratum"] == st]
                if not ss or qs <= 0: continue
                order = [max(ss, key=lambda i: cands[i]["novelty_vs_old"])]
                dmin_ = 1.0 - U[ss] @ U[order[0]]
                while len(order) < min(qs, len(ss)):
                    j = ss[int(np.argmax(dmin_))]; order.append(j); dmin_ = np.minimum(dmin_, 1.0 - U[ss] @ U[j])
                chosen += order
    with open(f"{R}/{OUT}/centres.jsonl", "w") as f:
        for i in chosen: f.write(json.dumps(dict(cands[i], candidate_index=int(i))) + "\n")
    json.dump(dict(mode="initial (TEST)" if a.test_initial else "md", frame_times_ps=FRAME_TIMES_PS, frames_used=frames_used,
                   detached_atoms_per_frame=detached_report), open(f"{R}/{OUT}/frames_used.json", "w"), indent=1)
    json.dump(dict(rejected_surface_atoms=dict(rejects), n_candidates=len(cands), n_rare_flagged=int(sum(c["rare_flag"] for c in cands)),
                   n_rare_chosen=int(sum(cands[i]["rare_flag"] for i in chosen))), open(f"{R}/{OUT}/rejects.json", "w"), indent=1)
    # summary
    by = collections.Counter((cands[i]["cls"], cands[i]["split"]) for i in chosen)
    bystr = collections.Counter(cands[i]["cn_stratum"] for i in chosen); poolstr = collections.Counter(c["cn_stratum"] for c in cands)
    nov = np.array([cands[i]["novelty_vs_old"] for i in chosen])
    L = [f"# Candidate centres ({'TEST on initial parents' if a.test_initial else 'MD frames'})", "",
         f"{len(chosen)} centres from a pool of {len(cands)} legal exposed atoms sampled on {sum(len(v) for v in frames_used.values())} frames "
         f"({len(frames_used)} parents x {FRAME_TIMES_PS if not a.test_initial else 'initial'}). Surface atoms rejected before sampling: {dict(rejects)}. "
         f"Rare-flagged (no pool neighbour at cosine > 0.98): {sum(c['rare_flag'] for c in cands)} in pool, {sum(cands[i]['rare_flag'] for i in chosen)} chosen (kept: geometry legal).", "",
         "| class | train | val | test |", "|---|---|---|---|"]
    for cls in "ABCD": L.append(f"| {cls} | {by[(cls,'train')]} | {by[(cls,'val')]} | {by[(cls,'test')]} |")
    L += ["", "| CN stratum | in pool | chosen |", "|---|---|---|"] + [f"| {s[0]} | {poolstr[s[0]]} | {bystr[s[0]]} |" for s in STRATA]
    L += ["", f"Novelty vs the 52 old geometries (1 - max cosine): median {np.median(nov):.4f}, 10th-90th pct {np.percentile(nov,10):.4f}-{np.percentile(nov,90):.4f}.",
          "Context (within 10 A), medians over chosen: " + ", ".join(f"{k} {np.median([cands[i][k] for i in chosen]):.2f}" for k in ("cn", "relief_10A", "n_under_10A", "n_levels_10A")),
          f"Detached atoms seen (parent@step: count): {detached_report if detached_report else 'none'}", "",
          "These are candidates for extraction; the final list is fixed only after repair, with the second-layer (whole-cell) stratification."]
    open(f"{R}/{OUT}/centres_summary.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
