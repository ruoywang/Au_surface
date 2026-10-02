#!/usr/bin/env python3
"""One unified compute list: unique entries keyed by (final geometry, calculation configuration, electronic state),
with the old state ids as aliases. Resolves (a) frozen size pairs whose two members resize to the SAME 6x6 cell and
potential (one entry, two aliases), (b) multi-layer geometries that appear both in the pair test and in the
replacement proposal with different potentials (the TEST potential is adopted, the proposal's value recorded),
(c) the two completed references and the six pair runs counted inside the 200. If the total exceeds 200 the most
redundant resized single-layer entries are dropped (nearest-neighbour distance in the frozen SOAP+morphology space);
nothing is added to reach exactly 200.

Order of execution: pairs (running) -> production 6x6 / 8x8 / multi-layer -> the 9 10x8-only entries last.
Nothing here submits anything.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh size_complex_test_v1/merge_lists.py
Output: size_complex_test_v1/unified_queue.csv, unified_queue.md
"""
import collections
import csv
import json
import os
import sys

import numpy as np

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"; T = f"{ROOT}/size_complex_test_v1"
MU0 = -4.9071
FIELDS = ["entry_id", "status", "order", "cls", "split", "cell", "n_atoms", "geometry_dir", "U_V", "TARGETMU_eV", "fermiconverge", "config", "aliases", "potential_note", "origin"]
COMPLEX_SPLIT = {"PM1d_s11": "train", "PM2a_s11": "train", "PM2b_s11": "val", "PM3a_s11": "train", "PM3b_s11": "val"}


def main():
    prop = list(csv.DictReader(open(f"{T}/replacement_proposal.csv"))); scr = {r["state_id"]: r for r in csv.DictReader(open(f"{T}/size_screen.frozen.csv"))}
    cscr = {r["state_id"]: r for r in csv.DictReader(open(f"{T}/size_screen.complex.csv"))}
    man = json.load(open(f"{R}/rough200/rough200_manifest.json")); S = {s["state_id"]: s for s in man["states"]}
    pq = json.load(open(f"{T}/pair_queue.json")); plan = json.load(open(f"{T}/pair_plan.json"))
    entries = collections.OrderedDict()

    def add(key, **kw):
        if key in entries: entries[key]["aliases"].append(kw.get("alias", "")); return entries[key]
        kw["aliases"] = [kw.pop("alias", "")]; entries[key] = kw; return kw
    # 1. completed references
    for t in json.load(open(f"{R}/dft/queue.json")):
        if t["status"] == "complete":
            s = S[t["state_id"]]
            add(f"{s['cell_id']}|production|{t['TARGETMU']:.4f}", status="complete", order=0, cls=s["cls"], split=s["split"], cell=s["cell"], n_atoms=t["n_atoms"], geometry_dir=s["cell_dir"],
                U_V=s["U_V"], TARGETMU_eV=t["TARGETMU"], fermiconverge=0.01, config="production", alias=t["state_id"], potential_note="as frozen", origin="frozen reference (8x8)")
    # 2. the six pair runs (test potential adopted; S1/S2 at the references' actual mu_e)
    for t in pq["tasks"]:
        if "dir" not in t: continue
        pr = next(p for p in plan["pairs"] if p["id"] == t["pair"]); cd = t["cell_dir"]; cid = os.path.basename(cd)
        U_actual = round(MU0 - t["TARGETMU"], 4)
        cls = pr["id"] if pr["id"].startswith("M") else cid[1]          # S1 -> A, S2 -> B (the reference's class)
        add(f"{cid}|{t['config']}|{t['TARGETMU']:.4f}", status=t["status"], order=1, cls=cls, split="test", cell=t["member"], n_atoms=t["n_atoms"],
            geometry_dir=cd, U_V=U_actual, TARGETMU_eV=t["TARGETMU"], fermiconverge=0.001, config=t["config"], alias=t["task_id"],
            potential_note=f"assigned U {pr['U_V']:+.2f} V; adopted TARGETMU {t['TARGETMU']:.4f} ({t['mu_source'][:60]}) -> actual U {U_actual:+.4f} V", origin=f"size pair {pr['id']}: {pr['target'][:60]}")
    # 3. production entries from the proposal
    for r in prop:
        s = S[r["state_id"]]; a = r["action"]
        if a == "kept_reference_8x8": continue
        if a == "release_for_complex":
            c = cscr[r["complex_candidate"]]; cd = c["new_cell_dir"]; cid = os.path.basename(cd); pid = c["parent_id"]
            note = r["reason"]; Uv = float(note.split("(U ")[1].split(" V")[0]) if "(U " in note else None; mu = round(MU0 - Uv, 4)
            # a geometry already covered by a pair run keeps the TEST potential: merge into that entry as an alias
            pair_keys = [k for k in entries if k.startswith(cid + "|")]
            if pair_keys:
                entries[pair_keys[0]]["aliases"].append(r["state_id"] + " (slot)"); entries[pair_keys[0]]["potential_note"] += f"; proposal would have given U {Uv:+.2f} V -- test potential kept, no second potential computed"; continue
            add(f"{cid}|production|{mu:.4f}", status="pending_production", order=2, cls=c["cls"], split=COMPLEX_SPLIT.get(pid, "train"), cell=c["preferred_cell"], n_atoms=int(c["new_n_atoms"]), geometry_dir=cd,
                U_V=Uv, TARGETMU_eV=mu, fermiconverge=0.01, config="production", alias=r["state_id"] + " (slot)", potential_note=f"inherits the released state's U bin; U {Uv:+.2f} V drawn", origin=f"multi-layer {c['cls']} {c['target_environment']} ({pid})")
            continue
        sc = scr[r["state_id"]]; mu = float(r["TARGETMU_eV"])
        if a == "resize_to_6x6":
            cd = sc["new_cell_dir"]; cid = os.path.basename(cd)
            e = add(f"{cid}|production|{mu:.4f}", status="pending_production", order=2, cls=r["cls"], split=r["split"], cell="6x6", n_atoms=int(sc["new_n_atoms"]), geometry_dir=cd, U_V=float(r["U_V"]), TARGETMU_eV=mu,
                    fermiconverge=0.01, config="production", alias=r["state_id"], potential_note="as frozen", origin="frozen centre resized to 6x6")
            if len(e["aliases"]) > 1: e["potential_note"] = "as frozen; two frozen size-pair members collapse to this one 6x6 entry"
        elif a == "keep_8x8":
            add(f"{s['cell_id']}|production|{mu:.4f}", status="pending_production", order=3, cls=r["cls"], split=r["split"], cell=s["cell"], n_atoms=s["n_atoms"], geometry_dir=s["cell_dir"], U_V=float(r["U_V"]), TARGETMU_eV=mu,
                fermiconverge=0.01, config="production", alias=r["state_id"], potential_note="as frozen", origin="frozen 8x8 kept: " + r["reason"][:70])
        elif a == "keep_10x8_no_smaller":
            add(f"{s['cell_id']}|production|{mu:.4f}", status="pending_last", order=9, cls=r["cls"], split=r["split"], cell=s["cell"], n_atoms=s["n_atoms"], geometry_dir=s["cell_dir"], U_V=float(r["U_V"]), TARGETMU_eV=mu,
                fermiconverge=0.01, config="production", alias=r["state_id"], potential_note="as frozen", origin="frozen 10x8 kept (no smaller valid cell); queued last")
    # 4a. exactly N_COMPLEX multi-layer entries INCLUDING the four pair states: beyond that, complex production entries are
    #     dropped (per class at most 5; within a class the later duplicates of an already covered target type on the same
    #     parent go first) and the single-layer state whose slot they took reverts to its resized 6x6 entry
    N_COMPLEX = 15
    complex_keys = [k for k, e in entries.items() if e["cls"].startswith("M")]
    if len(complex_keys) > N_COMPLEX:
        per_class = collections.defaultdict(list)
        for k in complex_keys: per_class[entries[k]["cls"]].append(k)
        to_drop = []
        for cls, keys in sorted(per_class.items()):
            prod = [k for k in keys if entries[k]["order"] == 2]
            seen = set(); ranked = []
            for k in prod:                       # first occurrence of a (parent, target) is kept preferentially
                tag = entries[k]["origin"].split("multi-layer ")[1] if "multi-layer " in entries[k]["origin"] else k
                ranked.append((tag in seen, k)); seen.add(tag)
            excess = len(keys) - 5
            for dup, k in sorted(ranked, key=lambda x: (not x[0], x[1]))[:max(0, excess)]: to_drop.append(k)
        for k in to_drop[: len(complex_keys) - N_COMPLEX]:
            e = entries.pop(k); slot_alias = next((a.replace(" (slot)", "") for a in e["aliases"] if "(slot)" in a), None)
            if slot_alias and slot_alias in scr and scr[slot_alias]["preferred_cell"] == "6x6":
                sc = scr[slot_alias]; cd = sc["new_cell_dir"]; cid = os.path.basename(cd); s = S[slot_alias]; mu = s["TARGETMU_eV"]
                add(f"{cid}|production|{mu:.4f}", status="pending_production", order=2, cls=s["cls"], split=s["split"], cell="6x6", n_atoms=int(sc["new_n_atoms"]), geometry_dir=cd, U_V=s["U_V"], TARGETMU_eV=mu,
                    fermiconverge=0.01, config="production", alias=slot_alias, potential_note="as frozen (slot returned: multi-layer quota is 15)", origin="frozen centre resized to 6x6")
    # 4b. the 200 cap: drop the most redundant resized single-layer entries if needed (never pad)
    n = len(entries); dropped = []
    if n > 200:
        z = np.load(f"{R}/rough200/rough200_soap.npz", allow_pickle=True); cid = list(z["cell_id"]); X = z["X"]; Mv = z["morph"]; idx = {c: k for k, c in enumerate(cid)}
        U = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12); Mz = (Mv - Mv.mean(0)) / np.maximum(Mv.std(0), 1e-9)
        cand = [(k, e) for k, e in entries.items() if e["origin"].startswith("frozen centre resized") and e["split"] == "train"]
        ks = [idx[S[e["aliases"][0]]["cell_id"]] for _, e in cand]
        Us, Ms = U[ks], Mz[ks]; D = 0.5 * (1 - Us @ Us.T) / 0.1 + 0.5 * np.linalg.norm(Ms[:, None] - Ms[None], axis=-1) / 3.0; np.fill_diagonal(D, np.inf); red = D.min(1)
        for o in np.argsort(red)[: n - 200]: dropped.append(cand[o][0])
        for k in dropped: del entries[k]
    rows = []
    for key, e in entries.items():
        rows.append(dict(entry_id=key, status=e["status"], order=e["order"], cls=e["cls"], split=e["split"], cell=e["cell"], n_atoms=e["n_atoms"], geometry_dir=e["geometry_dir"], U_V=e["U_V"],
                         TARGETMU_eV=e["TARGETMU_eV"], fermiconverge=e["fermiconverge"], config=e["config"], aliases=" | ".join(a for a in e["aliases"] if a), potential_note=e["potential_note"], origin=e["origin"]))
    rows.sort(key=lambda r: (r["order"], r["cls"], r["entry_id"]))
    with open(f"{T}/unified_queue.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader(); [w.writerow(r) for r in rows]
    c = collections.Counter(r["status"] for r in rows); merged = [r for r in rows if " | " in r["aliases"]]
    L = [f"# Unified compute list: {len(rows)} entries (cap 200) -- {dict(c)}", "",
         f"Merged entries (several old ids -> one geometry+potential): {len(merged)}"] + [f"- {r['entry_id']}: {r['aliases']}" for r in merged]
    L += ["", f"Dropped to respect the cap: {len(dropped)}" + ("; " + "; ".join(dropped) if dropped else ""),
          "", "| status | entries | atoms (sum) |", "|---|---|---|"] + [f"| {k} | {v} | {sum(int(r['n_atoms']) for r in rows if r['status'] == k)} |" for k, v in c.items()]
    by = collections.Counter((r["cls"], r["split"]) for r in rows); L += ["", "class x split: " + ", ".join(f"{k[0]}/{k[1]} {v}" for k, v in sorted(by.items()))]
    L += ["", "Pair states (test potential adopted; same-centre pairs share the test split): " + "; ".join(f"{r['entry_id'].split('|')[0]} U {float(r['U_V']):+.4f}" for r in rows if r["order"] == 1)]
    open(f"{T}/unified_queue.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
