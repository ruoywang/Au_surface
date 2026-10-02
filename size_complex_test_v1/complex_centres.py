#!/usr/bin/env python3
"""Representative centres on the multi-layer parents, each with a DECLARED target environment and the parent atom ids
that must be protected in the cut (not only the centre's 6 A neighbourhood). About four per parent, from two 300 K
frames, at distinct locations.

Targets (rules evaluated on the MD frame, so what is protected is what actually exists after the exploration):
  M1  upper_edge      exposed atom of a tier >= 5 with CN <= 7 (edge of an upper tier); protected += exposed atoms of
                      the tier below within 5 A laterally (the lower platform and the sidewall foot)
  M1  sidewall_foot   exposed layer-4 atom bonded to a layer-5 atom; protected += upper-edge atoms within 5 A
  M2  middle_terrace  exposed layer-4 atom with layer-5 atoms on one side and exposed layer-3 atoms on the other, both
                      within 7 A; protected += both step edges (layer-5 edge atoms and layer-4 atoms adjacent to layer 3)
                      within 7 A -- the relation 'two steps with a narrow terrace between' must survive the cut
  M2  double_step     exposed layer-5 atom with an exposed layer-3 atom within 4 A laterally (the two risers coincide)
  M3  valley_floor    exposed layer-3 atom with ridge atoms (layer >= 4) within 9 A in two roughly opposite directions
                      (angle > 120 deg); protected += the facing wall atoms of BOTH ridges within 9 A and the floor between
  M3  ridge_top_edge  exposed atom of the highest ridge level with CN <= 7 facing the valley; protected = core only
A rule that finds nothing on a frame is reported as such (e.g. a three-tier island that flattened during the MD is
recorded as what it became; the island top is never frozen to keep an appearance).

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh size_complex_test_v1/complex_centres.py [--per-parent 4] [--frames 165 199]
Output: size_complex_test_v1/complex_centres.jsonl, complex_centres_summary.md
"""
import argparse
import collections
import json
import os
import sys

import numpy as np
from ase.io import read

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"; T = f"{ROOT}/size_complex_test_v1"
sys.path.insert(0, R)
import select_centres as sc  # noqa: E402
import trajio  # noqa: E402
from extract_cells import layers_from_z, R_CORE, neighbours  # noqa: E402

D111 = sc.D111


def lateral_within(P, cell, i, idx, r):
    """indices among idx within r of atom i in-plane (minimum image)."""
    best = np.full(len(idx), np.inf)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = (si * cell[0] + sj * cell[1])[:2]
            best = np.minimum(best, np.linalg.norm(P[idx, :2] + sh - P[i, :2], axis=1))
    return idx[best < r], best


def vectors_to(P, cell, i, idx):
    """minimum-image in-plane vectors from atom i to atoms idx."""
    out = np.zeros((len(idx), 2)); best = np.full(len(idx), np.inf)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = (si * cell[0] + sj * cell[1])[:2]; v = P[idx, :2] + sh - P[i, :2]; d = np.linalg.norm(v, axis=1)
            m = d < best; out[m] = v[m]; best[m] = d[m]
    return out


def targets_for(at, cls):
    """list of (centre, target_name, protected_ids) found on this frame."""
    P = at.get_positions(); cell = at.get_cell().array; lay = layers_from_z(at)
    cn, dmin, adj = sc.neighbour_graph(at); surf = sc.surface_atoms(at); exposed = np.zeros(len(at), bool); exposed[surf] = True
    found = []

    def core(i): return set(int(x) for x in neighbours(at, i, R_CORE)[0]) | {int(i)}
    if cls == "M1":
        for i in surf:
            if lay[i] >= 5 and cn[i] <= 7:
                below = np.flatnonzero(exposed & (lay == lay[i] - 1)); near, _ = lateral_within(P, cell, i, below, 5.0)
                if len(near) >= 3: found.append((int(i), "upper_edge", core(i) | set(int(x) for x in near)))
        for i in surf:
            if lay[i] == 4 and any(lay[j] == 5 for j in adj[i]):
                up = np.flatnonzero(exposed & (lay == 5) & (cn <= 7)); near, _ = lateral_within(P, cell, i, up, 5.0)
                if len(near) >= 2: found.append((int(i), "sidewall_foot", core(i) | set(int(x) for x in near)))
    elif cls == "M2":
        l5 = np.flatnonzero(exposed & (lay == 5)); l3 = np.flatnonzero(exposed & (lay == 3)); l4e = np.flatnonzero((lay == 4) & exposed)
        for i in surf:
            if lay[i] == 4:
                n5, _ = lateral_within(P, cell, i, l5, 7.0); n3, _ = lateral_within(P, cell, i, l3, 7.0)
                if len(n5) >= 3 and len(n3) >= 3:
                    v5 = vectors_to(P, cell, i, n5).mean(0); v3 = vectors_to(P, cell, i, n3).mean(0)
                    if np.dot(v5, v3) < 0:                                     # steps on opposite sides of the terrace atom
                        edge4 = [int(j) for j in l4e if any(lay[k] == 3 for k in adj[j])]; e4, _ = lateral_within(P, cell, i, np.array(edge4, int), 7.0) if edge4 else (np.array([], int), None)
                        found.append((int(i), "middle_terrace", core(i) | set(int(x) for x in n5) | set(int(x) for x in e4)))
        for i in surf:
            if lay[i] == 5:
                n3, _ = lateral_within(P, cell, i, l3, 4.0)
                if len(n3) >= 1: found.append((int(i), "double_step", core(i)))
    elif cls == "M3":
        ridge = np.flatnonzero(exposed & (lay >= 4)); l3 = np.flatnonzero(exposed & (lay == 3)); Lmax = int(lay.max())
        for i in surf:
            if lay[i] == 3:
                near, _ = lateral_within(P, cell, i, ridge, 9.0)
                if len(near) >= 4:
                    v = vectors_to(P, cell, i, near); ang = np.arctan2(v[:, 1], v[:, 0])
                    spread = np.max([(np.cos(a - b) for a in ang for b in ang)]) if False else None
                    # two opposite walls: some pair of ridge atoms seen at > 120 degrees apart
                    cosmin = min(np.cos(a - b) for a in ang for b in ang)
                    if cosmin < -0.5:
                        floor, _ = lateral_within(P, cell, i, l3, 9.0)
                        found.append((int(i), "valley_floor", core(i) | set(int(x) for x in near) | set(int(x) for x in floor)))
        for i in surf:
            if lay[i] == Lmax and cn[i] <= 7 and Lmax >= 5:
                n3, _ = lateral_within(P, cell, i, l3, 9.0)
                if len(n3) >= 2: found.append((int(i), "ridge_top_edge", core(i)))
    return found


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--per-parent", type=int, default=4); ap.add_argument("--frames", type=float, nargs="*", default=[165.0, 199.0])
    a = ap.parse_args()
    trajio.MD_DIR = f"{T}/md"
    M = json.load(open(f"{T}/parents/parents_manifest.json")); out = []; report = []
    for m in M:
        pid = m["parent_id"]; cls = m["cls"]
        try: frames = trajio.frames_at(pid, a.frames)
        except trajio.TrajError as e: report.append(f"{pid}: {e}"); continue
        pool = []
        for step, t_ps, TK, at in frames:
            p0 = read(f"{T}/parents/{pid}.extxyz"); at.set_array("fixed", p0.get_array("fixed"))
            lay = layers_from_z(at); levels = np.bincount(lay).tolist()
            for centre, name, prot in targets_for(at, cls):
                pool.append(dict(parent_id=pid, cls=cls, step=int(step), time_ps=t_ps, T_K=TK, atom=centre, target_environment=name,
                                 protected_atom_ids=sorted(prot), n_protected=len(prot), cn=int(sc.neighbour_graph(at)[0][centre]) if False else None, levels=levels, _P=at.get_positions()[centre, :2].tolist()))
        names = collections.Counter(c["target_environment"] for c in pool)
        report.append(f"{pid} ({cls}): levels after MD {pool[0]['levels'] if pool else '-'}; targets found: {dict(names)}")
        # pick per parent: alternate target types, distinct locations (>= 12 A apart), both frames
        chosen = []
        steps_avail = sorted({c["step"] for c in pool})
        for name in sorted(names, key=lambda n: -names[n]):
            quota = max(1, a.per_parent // len(names) + 1); k = 0
            while len([x for x in chosen if x["target_environment"] == name]) < quota and k < 2 * quota:
                step = steps_avail[k % len(steps_avail)]; k += 1               # alternate frames: locations AND times differ
                cands = sorted([c for c in pool if c["target_environment"] == name and c["step"] == step], key=lambda c: -c["n_protected"])
                for c in cands:
                    if all(np.linalg.norm(np.array(c["_P"]) - np.array(x["_P"])) > 12.0 or c["step"] != x["step"] for x in chosen):
                        chosen.append(c); break
            if len(chosen) >= a.per_parent: break
        for c in chosen[:a.per_parent]:
            c = dict(c); c.pop("_P"); c.pop("cn"); c["candidate_id"] = f"{pid}_f{c['step']:06d}_a{c['atom']:04d}_{c['target_environment']}"
            c["parents_dir"] = f"{T}/parents"; c["md_dir"] = f"{T}/md"; out.append(c)
    with open(f"{T}/complex_centres.jsonl", "w") as f:
        for c in out: f.write(json.dumps(c) + "\n")
    L = ["# Centres on the multi-layer parents", "", f"{len(out)} candidates with declared targets and protected atom sets:", ""] + [f"- {r}" for r in report] + ["",
         "| candidate | target | protected atoms | frame (ps) |", "|---|---|---|---|"] + [f"| {c['candidate_id']} | {c['target_environment']} | {c['n_protected']} | {c['time_ps']:.0f} |" for c in out]
    open(f"{T}/complex_centres_summary.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
