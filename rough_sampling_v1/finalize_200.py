#!/usr/bin/env python3
"""Choose the rough states from the reconstructed cells and draw their potentials -- as a CANDIDATE list.

The manifest is written with frozen = false. `freeze` (run only after the user approves the list and the budget)
flips the flag; it never re-draws. rough_dft.py refuses to submit while frozen is false.

Two-layer selection (review of 2026-10-02):
  layer 1  the centre's local environment: SOAP of the centre atom in the RECONSTRUCTED cell (same settings as the
           centre selection) -- this covers new local environments and, by construction, barely changes if the core
           was kept; it therefore cannot by itself tell cells apart or vouch for the periphery;
  layer 2  the whole cell's upper surface: CN histogram of exposed atoms, relief, number of levels, adatom-level
           coverage, true terrace-vacancy fraction, exposed pit-floor fraction, number and size of islands (connected
           adatom-level components) and pits, step density, unregistered (mid-hop) atoms -- standardised and combined
           with layer 1 in the farthest-point distance;
  caps     per parent (<= 1.5 x its fair share of the bucket) and per (parent, frame) (<= 2);
  occupancy compression  lattice-site occupancy of layers 3 and 4, registered on the FROZEN bottom layer and the
           ABC stacking step measured from the two frozen layers (never on 'the first atom of the layer'), canonical
           under in-plane translations and by construction independent of atom order and periodic representation;
           this is verified on every cell (permutation, periodic re-representation, lattice and arbitrary
           translations) before it is used. Cells whose patterns differ in <= 2 sites (same size) are 'the same
           occupancy' -- one is kept. This compresses similar morphologies; it is NOT a claim that two such cells are
           identical geometries (strain and thermal displacements differ), and the dropped ones are listed as
           occupancy duplicates, not geometry duplicates.
  strata   CN strata of the centre get quotas proportional to sqrt(count).

Quotas: 160 train / 20 val / 20 test, 50 per class (40/5/5); 4 size pairs (same core in 8x8 and 10x8, same U) inside
the test split. Potentials: one U per geometry, 10 bins of 0.1 V over [-0.5, +0.5] V, drawn uniformly inside the
bin (rounded to 0.01 V) after the geometry is chosen, balanced over buckets, with a recorded seed; TARGETMU =
mu0 - U, mu0 = -4.9071 eV. With --keep-from <previous manifest>, cells that were already in the previous list keep
their U exactly; replacements take the bins vacated by the cells they replace (no global re-draw).

Usage (from Au_Cl/):
  rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/finalize_200.py [--cells cells] [--out rough200] [--n 200] [--seed 20261002]
                                                                   [--keep-from rough200/rough200_manifest.v1_37a9299.json]
  rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/finalize_200.py freeze [--out rough200]      # after approval only
Output: <out>/rough200_manifest.json (frozen=false), rough200_summary.md, rough200_cells.csv, rough200_soap.npz
"""
import argparse
import collections
import csv
import json
import os
import sys
import time

import numpy as np
from ase.io import read

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
sys.path.insert(0, R); sys.path.insert(0, f"{ROOT}/scripts")
import select_centres as sc  # noqa: E402
from extract_cells import layers_from_z, D111  # noqa: E402
import production as P  # noqa: E402  (kmesh rule only)

MU0 = -4.9071
U_MIN, U_MAX, N_BINS = -0.5, 0.5, 10
SPLIT_FRAC = {"train": 0.80, "val": 0.10, "test": 0.10}
N_PAIRS = 4
CAP_PER_FRAME = 2
CAP_PARENT_FACTOR = 1.5
DEDUP_SITES = 2
REG_TOL = 0.30         # of a primitive vector (~0.9 A): an atom farther from every site of its layer is 'unregistered' (mid-hop)
MORPH_KEYS = ["f_cn<=5", "f_cn6-7", "f_cn8-9", "f_cn>=10", "relief_layers", "n_levels", "adatom_level_ML", "terrace_vacancy_ML",
              "exposed_pit_ML", "n_islands", "max_island_frac", "n_pits", "step_density", "n_unregistered"]


# ------------------------------------------------------------------ registration and occupancy
def registration(at, n1, n2):
    """Fractional coordinates in the primitive basis, layer index, the bottom-layer sublattice offset f0 and the
    stacking step delta (+-1/3 per axis) measured from the two frozen layers. Order-independent (circular means)."""
    cell = at.get_cell().array; a1 = cell[0][:2] / n1; a2 = cell[1][:2] / n2
    f = at.get_positions()[:, :2] @ np.linalg.inv(np.array([a1, a2])); lay = layers_from_z(at)

    def cmean(x):
        x0 = x[0] % 1.0; d = ((x - x0 + 0.5) % 1.0) - 0.5; return (x0 + d.mean(0)) % 1.0
    f0 = cmean(f[lay == 0]); f1 = cmean(f[lay == 1]); delta = ((f1 - f0 + 0.5) % 1.0) - 0.5
    assert np.allclose(np.abs(delta), 1 / 3, atol=0.02), f"stacking step {delta} is not +-1/3"
    return f, lay, f0, delta


def registered_grids(at, n1, n2, tol=REG_TOL):
    """Site-occupancy grids of layers 2, 3, 4 (True = a registered atom on that site), the number of unregistered
    atoms, and the atom count per layer."""
    f, lay, f0, delta = registration(at, n1, n2)
    grids, unreg, counts = {}, 0, {}
    for L in (2, 3, 4):
        g = np.zeros((n1, n2), bool); m = lay == L; counts[L] = int(m.sum())
        if m.any():
            r = f[m] - (f0 + L * delta); k = np.round(r); res = np.abs(r - k).max(1); ok = res <= tol
            unreg += int((~ok).sum()); idx = k[ok].astype(int); g[idx[:, 0] % n1, idx[:, 1] % n2] = True
        grids[L] = g
    return grids, unreg, counts


def all_rolls(grids):
    g3, g4 = grids[3], grids[4]; n1, n2 = g3.shape
    out = np.zeros((n1 * n2, 2 * n1 * n2), bool); k = 0
    for di in range(n1):
        for dj in range(n2):
            out[k] = np.concatenate([np.roll(g3, (di, dj), (0, 1)).ravel(), np.roll(g4, (di, dj), (0, 1)).ravel()]); k += 1
    return out


def canonical(grids):
    return min(np.packbits(r).tobytes() for r in all_rolls(grids))


def min_hamming(rollsA, patternB):
    return int((rollsA != patternB[None, :]).sum(1).min())


def fingerprint_invariance(at, n1, n2, seed=0, trials=3):
    """The canonical occupancy must not change under: atom permutation, moving a third of the atoms to periodic
    equivalents, a whole-cell lattice translation and an arbitrary in-plane translation."""
    base = canonical(registered_grids(at, n1, n2)[0]); rng = np.random.default_rng(seed); n = len(at); cell = at.get_cell().array
    for _ in range(trials):
        t = at[rng.permutation(n)]; Pp = t.get_positions()
        sel = rng.choice(n, size=max(1, n // 3), replace=False); k = rng.integers(-1, 2, size=(len(sel), 2))
        Pp[sel] += k[:, :1] * cell[0] + k[:, 1:2] * cell[1]
        Pp += rng.integers(-2, 3) * cell[0] + rng.integers(-2, 3) * cell[1] + np.r_[rng.uniform(-5, 5, 2), 0.0]
        t.set_positions(Pp)
        if canonical(registered_grids(t, n1, n2)[0]) != base: return False
    return True


# ------------------------------------------------------------------ morphology
def components(members, adj):
    members = set(int(m) for m in members); seen = set(); comps = []
    for s in members:
        if s in seen: continue
        stack = [s]; seen.add(s); comp = [s]
        while stack:
            i = stack.pop()
            for j in adj[i]:
                if j in members and j not in seen: seen.add(j); stack.append(j); comp.append(j)
        comps.append(comp)
    return comps


def morphology(at, n1, n2):
    sites = n1 * n2
    cn, dmin, adj = sc.neighbour_graph(at); surf = sc.surface_atoms(at); lay = layers_from_z(at)
    grids, unreg, counts = registered_grids(at, n1, n2)
    s_cn = cn[surf]; z = at.get_positions()[surf, 2]
    isl = components(np.flatnonzero(lay == 4), adj); pits = components(np.flatnonzero((lay == 2) & np.isin(np.arange(len(at)), surf)), adj)
    return {"f_cn<=5": float((s_cn <= 5).mean()), "f_cn6-7": float(((s_cn >= 6) & (s_cn <= 7)).mean()), "f_cn8-9": float(((s_cn >= 8) & (s_cn <= 9)).mean()),
            "f_cn>=10": float((s_cn >= 10).mean()), "relief_layers": float(np.ptp(z) / D111), "n_levels": int(len(np.unique(np.round(z / D111)))),
            "adatom_level_ML": float(counts[4] / sites),
            "layer3_deficit": float(1 - counts[3] / sites),                       # count-based; negative when mid-hop atoms are counted in layer 3
            "terrace_vacancy_ML": float(1 - grids[3].sum() / sites),              # true unoccupied layer-3 sites
            "exposed_pit_ML": float(((lay == 2) & np.isin(np.arange(len(at)), surf)).sum() / sites), "n_islands": len(isl),
            "max_island_frac": float(max((len(c) for c in isl), default=0) / sites), "n_pits": len(pits), "step_density": float((s_cn <= 8).mean()),
            "n_unregistered": int(unreg)}


def allocate(q, sizes): return sc.allocate(q, sizes)


# ------------------------------------------------------------------ main
def freeze(args):
    man = f"{R}/{args.out}/rough200_manifest.json"; m = json.load(open(man))
    if m.get("frozen"): sys.exit("already frozen")
    m["frozen"] = True; m["frozen_at"] = time.strftime("%Y-%m-%d %H:%M"); m["status"] = "FROZEN (approved list; potentials as drawn, not re-drawn)"
    json.dump(m, open(man, "w"), indent=1); print(f"frozen: {m['n_states']} states, seed {m['seed']}")


def main(args):
    OUT = f"{R}/{args.out}"; os.makedirs(OUT, exist_ok=True); man = f"{OUT}/rough200_manifest.json"
    if os.path.exists(man) and json.load(open(man)).get("frozen") and not args.force:
        sys.exit(f"{man} is FROZEN; re-running would re-select. Use --force only deliberately and say so in the report.")
    prev = json.load(open(f"{R}/{args.keep_from}")) if args.keep_from else None
    rows = [json.loads(l) for l in open(f"{R}/{args.cells}/cells_manifest.jsonl")]
    ok = [r for r in rows if r["status"] == "PASS"]
    if not ok: sys.exit("no PASS cells")
    # ---- descriptors on the reconstructed cells; fingerprint invariance verified on EVERY cell
    sc.OUT = "centres"
    Xref = sc.old_reference(); Uref = Xref / np.maximum(np.linalg.norm(Xref, axis=1, keepdims=True), 1e-12)
    X, morph, rolls, kmesh, inv_fail = [], [], [], [], []
    for r in ok:
        at = read(f"{R}/{args.cells}/{r['cell_id']}/cell.extxyz"); n1, n2 = (int(x) for x in r["cell"].split("x"))
        X.append(sc.soap_of(at, [r["centre_in_cell"]])[0]); morph.append(morphology(at, n1, n2))
        rolls.append(all_rolls(registered_grids(at, n1, n2)[0]))
        if not fingerprint_invariance(at, n1, n2, seed=len(rolls)): inv_fail.append(r["cell_id"])
        kp, n = P.kmesh(at); kmesh.append(f"{n[0]}x{n[1]}x1")
    if inv_fail: sys.exit(f"occupancy fingerprint NOT invariant for {len(inv_fail)} cells, e.g. {inv_fail[:5]}; fix before selecting")
    X = np.vstack(X); U = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    novelty = 1.0 - (U @ Uref.T).max(axis=1)
    Mv = np.array([[m[k] for k in MORPH_KEYS] for m in morph], float)
    Mz = (Mv - Mv.mean(0)) / np.maximum(Mv.std(0), 1e-9)
    for i, r in enumerate(ok): r["novelty_cell_vs_old"] = float(novelty[i]); r["morphology"] = morph[i]; r["kpoints"] = kmesh[i]
    np.savez_compressed(f"{OUT}/rough200_soap.npz", X=X, morph=Mv, cell_id=np.array([r["cell_id"] for r in ok]))
    # ---- occupancy compression (same cell size; patterns within DEDUP_SITES sites under translation)
    dropped = {}
    by_size = collections.defaultdict(list)
    for i, r in enumerate(ok): by_size[r["cell"]].append(i)
    for size, idx in by_size.items():
        idx = sorted(idx, key=lambda i: -novelty[i]); kept = []
        for i in idx:
            pat = rolls[i][0]
            dup = next((j for j in kept if abs(ok[i]["n_atoms"] - ok[j]["n_atoms"]) <= DEDUP_SITES and min_hamming(rolls[j], pat) <= DEDUP_SITES), None)
            if dup is None: kept.append(i)
            else: dropped[ok[i]["cell_id"]] = ok[dup]["cell_id"]
    too_seamy = {r["cell_id"] for r in ok if r.get("stats", {}).get("n_seam_affected", 0) > args.max_seam}
    eligible = [i for i in range(len(ok)) if ok[i]["cell_id"] not in dropped and ok[i]["cell_id"] not in too_seamy]
    # ---- quotas
    per_class = args.n // 4; quota = {}
    for cls in "ABCD":
        q = {sp: int(round(per_class * f)) for sp, f in SPLIT_FRAC.items()}; q["train"] += per_class - sum(q.values())
        for sp, v in q.items(): quota[(cls, sp)] = v
    # ---- size pairs
    by_centre = collections.defaultdict(dict)
    for i in eligible: by_centre[ok[i]["centre_id"]][ok[i]["cell"]] = i
    pair_cands = {cls: [c for c, d in by_centre.items() if len(d) == 2 and ok[d["8x8"]]["split"] == "test" and ok[d["8x8"]]["cls"] == cls] for cls in "ABCD"}
    chosen, used_centres, pairs, pair_members = [], set(), [], {}
    order_cls = [c for c in "ABCD" if pair_cands[c]]
    while len(pairs) < N_PAIRS and order_cls:
        cls = order_cls.pop(0)
        cands = [c for c in pair_cands[cls] if c not in used_centres]
        if not cands: continue
        c = max(cands, key=lambda c_: novelty[by_centre[c_]["8x8"]]); pairs.append((cls, c)); used_centres.add(c)
        for size, i in by_centre[c].items(): chosen.append(i); pair_members[i] = f"pair{len(pairs)}"
        quota[(cls, "test")] -= 2
        if not order_cls and len(pairs) < N_PAIRS: order_cls = [k for k in "ABCD" if any(c_ not in used_centres for c_ in pair_cands[k])]
    # ---- two-layer farthest-point per bucket with caps
    shortfall = {}; cap_report = {}
    for (cls, sp), q in quota.items():
        if q <= 0: continue
        sub = [i for i in eligible if ok[i]["cls"] == cls and ok[i]["split"] == sp and ok[i]["centre_id"] not in used_centres]
        seen = set(); sub2 = []
        for i in sorted(sub, key=lambda i: (ok[i]["centre_id"], ok[i]["cell"] != "8x8")):
            if ok[i]["centre_id"] in seen: continue
            seen.add(ok[i]["centre_id"]); sub2.append(i)
        if not sub2: shortfall[(cls, sp)] = q; continue
        parents = {ok[i]["parent_id"] for i in sub2}
        cap_parent = min(q, int(np.ceil(CAP_PARENT_FACTOR * q / len(parents)))); cap_report[f"{cls}/{sp}"] = dict(parents=len(parents), cap_per_parent=cap_parent)
        Us = U[sub2]; Ms = Mz[sub2]
        d_soap = 1.0 - Us @ Us.T; d_mor = np.linalg.norm(Ms[:, None, :] - Ms[None, :, :], axis=-1)
        iu = np.triu_indices(len(sub2), 1)
        s_soap = max(np.percentile(d_soap[iu], 90), 1e-9) if len(sub2) > 1 else 1.0
        s_mor = max(np.percentile(d_mor[iu], 90), 1e-9) if len(sub2) > 1 else 1.0
        D = 0.5 * d_soap / s_soap + 0.5 * d_mor / s_mor
        sizes = {s[0]: sum(ok[i]["cn_stratum"] == s[0] for i in sub2) for s in sc.STRATA}
        n_parent = collections.Counter(); n_frame = collections.Counter(); taken = []
        for st, qs in allocate(q, sizes).items():
            ss = [k for k, i in enumerate(sub2) if ok[i]["cn_stratum"] == st]
            if not ss or qs <= 0: continue
            dmin = np.full(len(sub2), np.inf); picked = 0
            first = max(ss, key=lambda k: novelty[sub2[k]]); cand_order = [first]
            while picked < qs:
                best = None
                for k in (cand_order if not taken else sorted(ss, key=lambda k: -dmin[k])):
                    i = sub2[k]
                    if k in taken: continue
                    if n_parent[ok[i]["parent_id"]] >= cap_parent or n_frame[(ok[i]["parent_id"], ok[i]["step"])] >= CAP_PER_FRAME: continue
                    best = k; break
                if best is None: break
                taken.append(best); picked += 1; i = sub2[best]
                n_parent[ok[i]["parent_id"]] += 1; n_frame[(ok[i]["parent_id"], ok[i]["step"])] += 1
                dmin = np.minimum(dmin, D[best]); cand_order = []
        got = [sub2[k] for k in taken]
        if len(got) < q: shortfall[(cls, sp)] = q - len(got)
        chosen += got; used_centres.update(ok[i]["centre_id"] for i in got)
    # ---- potentials: kept cells keep their U; new units take the bins vacated by the replaced ones
    rng = np.random.default_rng(args.seed); edges = np.linspace(U_MIN, U_MAX, N_BINS + 1)
    units, seen_pairs = [], set()
    for i in chosen:
        if i in pair_members:
            if pair_members[i] in seen_pairs: continue
            seen_pairs.add(pair_members[i]); units.append([j for j in chosen if pair_members.get(j) == pair_members[i]])
        else: units.append([i])
    assign = {}; bin_counts = np.zeros(N_BINS, int); replaced_out = []; new_in = []
    if prev is not None:
        prevU = {s["cell_id"]: (s["U_V"], s["U_bin"]) for s in prev["states"]}
        chosen_ids = {ok[i]["cell_id"] for i in chosen}
        kept_units, new_units = [], []
        for u in units:
            ids = [ok[i]["cell_id"] for i in u]
            if all(c in prevU for c in ids) and len({prevU[c][0] for c in ids}) == 1: kept_units.append(u)
            else: new_units.append(u)
        for u in kept_units:
            Uv, ub = prevU[ok[u[0]]["cell_id"]]; b = int(np.argmin(np.abs(edges[:-1] - ub[0]))); bin_counts[b] += 1
            for i in u: assign[i] = (b, Uv)
        # vacated bins: previous units no longer chosen (a previous pair counts once)
        prev_units = {}
        for s in prev["states"]:
            key = s["pair"] or s["cell_id"]
            if s["cell_id"] not in chosen_ids: prev_units.setdefault(key, s["U_bin"][0])
        vac = [int(np.argmin(np.abs(edges[:-1] - lo))) for lo in prev_units.values()]
        replaced_out = [s["cell_id"] for s in prev["states"] if s["cell_id"] not in chosen_ids]
        rng.shuffle(vac)
        for u in new_units:
            b = vac.pop() if vac else int(rng.choice(np.flatnonzero(bin_counts == bin_counts.min()))); bin_counts[b] += 1
            Uv = float(np.round(np.clip(np.round(rng.uniform(edges[b], edges[b + 1]), 2), edges[b] + 0.005, edges[b + 1] - 0.005), 2))
            for i in u: assign[i] = (b, Uv); new_in.append(ok[i]["cell_id"])
    else:
        buckets = collections.defaultdict(list)
        for u in units: buckets[(ok[u[0]]["cls"], ok[u[0]]["split"])].append(u)
        bkeys = sorted(buckets); rng.shuffle(bkeys)
        for bk in bkeys:
            us = sorted(buckets[bk], key=lambda u: ok[u[0]]["cell_id"]); rng.shuffle(us)
            for u in us:
                b = int(rng.choice(np.flatnonzero(bin_counts == bin_counts.min()))); bin_counts[b] += 1
                Uv = float(np.round(np.clip(np.round(rng.uniform(edges[b], edges[b + 1]), 2), edges[b] + 0.005, edges[b + 1] - 0.005), 2))
                for i in u: assign[i] = (b, Uv)
    # ---- manifest
    states = []
    for i in sorted(chosen, key=lambda i: (ok[i]["split"], ok[i]["cls"], ok[i]["cell_id"])):
        r = ok[i]; b, Uv = assign[i]
        states.append(dict(state_id=f"{r['cell_id']}__U{Uv:+.2f}", cell_id=r["cell_id"], centre_id=r["centre_id"], parent_id=r["parent_id"], cls=r["cls"],
                           split=r["split"], cell=r["cell"], n_atoms=r["n_atoms"], kpoints=r["kpoints"], step=r["step"], time_ps=r.get("time_ps"), source=r["source"],
                           cn=r["cn"], cn_stratum=r["cn_stratum"], rare_flag=r.get("rare_flag"), U_V=Uv, U_bin=[float(edges[b]), float(edges[b + 1])],
                           TARGETMU_eV=round(MU0 - Uv, 4), pair=pair_members.get(i), repair=r["repair"], novelty_cell_vs_old=r["novelty_cell_vs_old"],
                           n_seam_affected=r.get("stats", {}).get("n_seam_affected"), morphology=r["morphology"], cell_dir=f"{R}/{args.cells}/{r['cell_id']}"))
    info = dict(created=time.strftime("%Y-%m-%d %H:%M"), status="CANDIDATE (not frozen; nothing may be submitted)", frozen=False, seed=args.seed, mu0_eV=MU0,
                n_states=len(states), bins_V=edges.tolist(), quota={f"{k[0]}/{k[1]}": v for k, v in quota.items()},
                shortfall={f"{k[0]}/{k[1]}": v for k, v in shortfall.items()}, caps=dict(per_frame=CAP_PER_FRAME, per_parent=cap_report, max_seam_affected=args.max_seam),
                excluded_too_seamy=sorted(too_seamy), occupancy_duplicates=dropped,
                fingerprint=dict(registration="frozen bottom layer + stacking step from the two frozen layers", tol_sites=REG_TOL, cells_checked=len(ok),
                                 invariance_failures=0, checks="atom permutation, periodic re-representation of a third of the atoms, lattice and arbitrary translation"),
                kept_from=args.keep_from, replaced_out=replaced_out, new_in=new_in,
                pairs=[dict(id=f"pair{k+1}", cls=c, centre_id=cc) for k, (c, cc) in enumerate(pairs)],
                rule="two-layer farthest-point (centre SOAP + whole-cell morphology) per class x split x CN stratum with parent/frame caps and occupancy "
                     "compression; potential drawn afterwards, once (kept cells keep theirs); TARGETMU is the request, the label is the converged mu_e",
                states=states)
    json.dump(info, open(man, "w"), indent=1)
    with open(f"{OUT}/rough200_cells.csv", "w", newline="") as f:
        w = csv.writer(f); keys = MORPH_KEYS + ["layer3_deficit"]
        w.writerow(["state_id", "cls", "split", "cell", "n_atoms", "kpoints", "parent_id", "step", "cn_stratum", "U_V", "TARGETMU_eV", "pair", "repair", "n_seam_affected"] + keys)
        for s in states: w.writerow([s["state_id"], s["cls"], s["split"], s["cell"], s["n_atoms"], s["kpoints"], s["parent_id"], s["step"], s["cn_stratum"], s["U_V"], s["TARGETMU_eV"],
                                     s["pair"] or "", s["repair"], s["n_seam_affected"]] + [round(s["morphology"][k], 3) for k in keys])
    summary(info, states, ok, rows, dropped, OUT, args, prev)


def summary(info, states, ok, rows, dropped, OUT, args, prev):
    cnt = collections.Counter((s["cls"], s["split"]) for s in states); hist = collections.Counter(tuple(s["U_bin"]) for s in states)
    edges = info["bins_V"]
    L = [f"# Rough states: CANDIDATE list of {len(states)} (seed {info['seed']}, frozen = {info['frozen']})", "",
         "Nothing here is submitted or frozen. Approval = `finalize_200.py freeze`, then `rough_dft.py submit --first 8 --confirm`.", "",
         f"Occupancy fingerprint: registered on the frozen bottom layer and the measured stacking step; invariance (permutation, periodic "
         f"re-representation, translations) verified on all {info['fingerprint']['cells_checked']} passing cells, 0 failures.", ""]
    if prev is not None:
        L += ["## Changes against the previous candidate list", "",
              f"Kept from `{info['kept_from']}`: {len(states) - len(info['new_in'])} states with their potentials unchanged. Replaced: {len(info['replaced_out'])} out, "
              f"{len(info['new_in'])} in (new cells took the bins vacated by the replaced ones).", ""]
        if info["replaced_out"]:
            L += ["| out (previous) | in (new) |", "|---|---|"]
            for a_, b_ in zip(info["replaced_out"] + [""] * max(0, len(info["new_in"]) - len(info["replaced_out"])), info["new_in"] + [""] * max(0, len(info["replaced_out"]) - len(info["new_in"]))):
                L.append(f"| {a_} | {b_} |")
            pc_prev = collections.Counter((s["cls"], s["split"]) for s in prev["states"]); st_prev = collections.Counter(s["cn_stratum"] for s in prev["states"])
            st_now = collections.Counter(s["cn_stratum"] for s in states)
            L += ["", "Coverage before -> after: class x split " + ("unchanged" if pc_prev == cnt else f"{dict(pc_prev)} -> {dict(cnt)}") +
                  f"; CN strata {dict(st_prev)} -> {dict(st_now)}; cells per size {dict(collections.Counter(s['cell'] for s in prev['states']))} -> "
                  f"{dict(collections.Counter(s['cell'] for s in states))}."]
        L.append("")
    L += ["## Where they come from", "", "| class | train | val | test |", "|---|---|---|---|"]
    for cls in "ABCD": L.append(f"| {cls} | {cnt[(cls,'train')]} | {cnt[(cls,'val')]} | {cnt[(cls,'test')]} |")
    pc = collections.Counter(s["parent_id"] for s in states); ps = collections.Counter((s["parent_id"], s["step"]) for s in states)
    L += ["", f"Parents used: {len(pc)} of 32; cells per parent min/median/max {min(pc.values())}/{int(np.median(list(pc.values())))}/{max(pc.values())} "
          f"(caps: {json.dumps(info['caps']['per_parent'])}); cells per (parent, frame) max {max(ps.values())} (cap {CAP_PER_FRAME}).",
          "Frames (ps): " + ", ".join(f"{t:.0f} ps: {n}" for t, n in sorted(collections.Counter(s.get('time_ps') or 0 for s in states).items())),
          "", "| CN stratum of the centre | states |", "|---|---|"] + [f"| {st[0]} | {sum(s['cn_stratum'] == st[0] for s in states)} |" for st in sc.STRATA]
    L += ["", "## Size, atoms, k-mesh", "", "Cell sizes: " + str(dict(collections.Counter(s["cell"] for s in states))),
          "Atoms: " + ", ".join(f"{k}: {v}" for k, v in sorted(collections.Counter(s["n_atoms"] for s in states).items())),
          "k-mesh: " + str(dict(collections.Counter(s["kpoints"] for s in states))),
          "Size pairs (test, same U): " + (", ".join(f"{p['id']} {p['cls']} {p['centre_id']}" for p in info["pairs"]) or "none"),
          "Repair: " + str(dict(collections.Counter("no repair" if s["repair"].startswith("no repair") else "seam band minimised" for s in states)))]
    L += ["", "## Morphology coverage of the chosen cells (medians and ranges)", "", "| descriptor | min | median | max |", "|---|---|---|---|"]
    for k in MORPH_KEYS + ["layer3_deficit"]:
        v = np.array([s["morphology"][k] for s in states]); L.append(f"| {k} | {v.min():.3f} | {np.median(v):.3f} | {v.max():.3f} |")
    L += ["", "`terrace_vacancy_ML` counts truly unoccupied layer-3 sites (>= 0); `layer3_deficit` is the count-based 1 - N3/sites and goes negative when an "
          "atom caught between levels is counted in layer 3; `n_unregistered` counts atoms farther than 0.3 of a site spacing from every site of their layer."]
    L += ["", "## Potentials", "", "| U bin (V) | states |", "|---|---|"] + [f"| {lo:+.1f} to {hi:+.1f} | {hist[(lo,hi)]} |" for lo, hi in zip(edges[:-1], edges[1:])]
    fails = collections.Counter(p for r in rows if r["status"] == "FAIL" for t in r["trials"] for p in t["problems"])
    n_fail_centres = sum(r["status"] == "FAIL" for r in rows)
    L += ["", "## Rejections", "", f"Centres with no passing cell: {n_fail_centres} of {n_fail_centres + len({r['centre_id'] for r in ok})}. Reasons over all failed trials:"]
    L += [f"- {k}: {v}" for k, v in fails.most_common(10)]
    L += ["", f"Occupancy duplicates compressed (layer-3/4 site occupancy within {DEDUP_SITES} sites under translation, same size; similar morphology, NOT a claim of "
          f"identical geometry): {len(dropped)}" + ("; e.g. " + "; ".join(f"{a} ~ {b}" for a, b in list(dropped.items())[:5]) if dropped else "")]
    seam = np.array([s["n_seam_affected"] or 0 for s in states]); seam_all = np.array([r.get("stats", {}).get("n_seam_affected", 0) for r in ok])
    L += [f"Seam-affected atoms (CN reduced by >= 2 vs the parent) per chosen cell: median {np.median(seam):.0f}, max {seam.max()} "
          f"(all passing cells: median {np.median(seam_all):.0f}, max {seam_all.max()}; cells above the cap of {args.max_seam} excluded: {len(info['excluded_too_seamy'])})."]
    if info["shortfall"]: L += ["", "SHORTFALL (bucket: missing): " + json.dumps(info["shortfall"]) + " -- fewer eligible cells than the quota under the caps; not filled from other buckets."]
    L += ["", f"Novelty of the chosen cells' centres vs the 52 old geometries (1 - max cosine): median {np.median([s['novelty_cell_vs_old'] for s in states]):.4f}; "
          f"rare-flagged centres chosen: {sum(bool(s['rare_flag']) for s in states)}."]
    open(f"{OUT}/rough200_summary.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", nargs="?", default="select", choices=["select", "freeze"])
    ap.add_argument("--cells", default="cells"); ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261002); ap.add_argument("--force", action="store_true"); ap.add_argument("--out", default="rough200")
    ap.add_argument("--max-seam", type=int, default=16, help="exclude cells with more seam-affected atoms than this (a quarter of a 64-site layer)")
    ap.add_argument("--keep-from", default=None, help="previous manifest (path under rough_sampling_v1): kept cells keep their U; replacements take vacated bins")
    a = ap.parse_args()
    freeze(a) if a.cmd == "freeze" else main(a)
