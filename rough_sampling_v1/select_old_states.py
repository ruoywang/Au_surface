#!/usr/bin/env python3
"""Select ~140 of the 505 existing DFT states for the rough_sampling_v1 working set, with a reason per row.

Budget (a selection budget, not a quota to be filled):
    16 base-reference geometries   x all available potentials (<= 5)   ~80 states
    24 non-equilibrium geometries  x 2 potentials, preferring +-0.5 V   ~48
    12 supplementary references    x 1 potential                        ~12

Rules, each of which produces a written reason in the output:
  * one actual geometry per base structure, the accepted relaxation preferred, EXCEPT where relaxing destroyed
    the environment the structure exists for (checked, not assumed: R1 must still be hcp-registered, R2 must
    still carry its domain walls, otherwise the ideal geometry is taken and the reason says so);
  * two relaxed geometries that are the same structure under different names are not both kept: every pair of
    base geometries in the same cell with the same atom count is compared by RMSD after translation;
  * non-equilibrium geometries are chosen by farthest-point sampling on a small descriptor of what the
    deformation did (type, family, displacement, closest contact, under-coordinated count), under a family
    cover, not by taking two perturbations of every parent;
  * a state enters only with complete labels. Field completeness is judged per field against siblings on the
    same grid, so a truncated RHOB disqualifies the state from the full-label set without hiding that its
    N_e, mu_e, CHGCAR, PHI and SION are fine; that is recorded in label_availability.jsonl for all 505.

Energy fields keep their original names (E_free_TOTEN_eV, GCE_code_eV, ...). Nothing is renamed "energy".
Nothing is deleted or re-weighted; unselected states stay in the archive untouched.

Usage (from Au_Cl/):  scripts/pyrun.sh rough_sampling_v1/select_old_states.py
Outputs: rough_sampling_v1/old_selected_140.csv, old_excluded.csv, label_availability.jsonl, selection_report.md
"""
import collections
import csv
import json
import os
import sys

import numpy as np
from ase.io import read

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
OUT = f"{ROOT}/rough_sampling_v1"
sys.path.insert(0, f"{ROOT}/scripts")
from analysis_spatial import coordination  # noqa: E402

MU0 = -4.9071
BASE16 = ["T-4x4", "V1", "A1-fcc", "A1-hcp", "Step-8x2", "Step-16x1", "Au211", "Au221", "Kink-edge1",
          "Kink-edge2", "Island-7-compact", "Island-7-elongated", "Pit-7-compact", "Pit-7-trench",
          "R1-hcp-terminated", "R2-stripe-wall"]
SUPP12 = ["V2", "V3", "A3", "Au332", "Au554", "Step-24x1", "C1-island-near-step", "C2-island+pit",
          "Island-7-8x8", "Island-19-8x8", "Pit-7-8x8", "Pit-19-8x8"]
# strictly periodic copies of a base structure: the same local environment in a bigger cell
PERIODIC_COPIES = {"Flat-8x2": "T-4x4", "Flat-16x1": "T-4x4", "Step-8x1": "Step-8x2", "Step-16x2": "Step-16x1"}
N_NONEQ_PERT, N_NONEQ_PATH = 16, 8
RMSD_DUP = 0.15          # A; two "different" base geometries closer than this are the same structure


# ----------------------------------------------------------------------------------------- label availability
def field_completeness(states):
    """Per field, per state: complete if its byte size equals the largest seen on the same grid for the same
    structure. Truncation shows up as a smaller file; absence as no entry."""
    ref = collections.defaultdict(int)
    for s in states:
        for name, f in (s.get("fields") or {}).items():
            key = (s["structure_id"], name, tuple(f.get("grid") or ()))
            ref[key] = max(ref[key], f.get("bytes") or 0)
    rows = []
    for s in states:
        fields = s.get("fields") or {}
        per = {}
        for name, f in fields.items():
            key = (s["structure_id"], name, tuple(f.get("grid") or ()))
            b = f.get("bytes") or 0
            per[name] = "complete" if b == ref[key] else f"truncated ({b}/{ref[key]} bytes)"
        expected = {"CHGCAR", "LOCPOT", "PHI", "SION", "RHOB", "RHOION"}
        for name in expected - set(fields):
            per[name] = "absent"
        e = s["electronic_state"]; lab = s["labels"]
        rows.append(dict(
            state_id=s["state_id"], geometry_id=s["geometry_id"], structure_id=s["structure_id"],
            config=s["config"], campaign=s["campaign"],
            U_V=round(e["mu_reference_eV"] - e["mu_e_actual_eV"], 4), TARGETMU_eV=e["TARGETMU_eV"],
            N_e=e["N_e_final"], mu_e_actual_eV=e["mu_e_actual_eV"],
            mu_within_tolerance=bool((s.get("qc") or {}).get("mu_within_tolerance", True)),
            energy_fields=dict(E_free_TOTEN_eV=lab.get("E_free_TOTEN_eV"), E_sigma0_eV=lab.get("E_sigma0_eV"),
                               GCE_code_eV=lab.get("GCE_code_eV"), status=lab.get("energy_label_status")),
            forces_present=lab.get("forces_eV_per_A") is not None,
            max_force_movable_eV_per_A=lab.get("max_force_movable_eV_per_A"),
            fields=per,
            all_fields_complete=all(v == "complete" for v in per.values()),
            qc_status=(s.get("qc") or {}).get("status"),
        ))
    return rows


# ----------------------------------------------------------------------------------------- geometry helpers
def load_geom(s):
    return read(f"{s['source_dir']}/{s['geometry_file']}")


def rmsd_translated(a, b):
    """RMSD between two same-size structures after removing the mean translation, minimum image in-plane."""
    pa, pb = a.get_positions(), b.get_positions()
    if len(pa) != len(pb): return np.inf
    cell = a.get_cell().array; C = np.array([cell[0][:2], cell[1][:2]]); Ci = np.linalg.inv(C)
    d = pa - pb
    f = d[:, :2] @ Ci; f -= np.round(f); d[:, :2] = f @ C
    d -= d.mean(axis=0)
    return float(np.sqrt((d ** 2).sum(axis=1).mean()))


def registry_kept(sid, geom):
    """Did relaxing R1 keep its hcp termination, and R2 its two walls? Measured, using the gallery's own
    registry classifier so the test is the one the figures use."""
    import render_gallery as rg
    at = rg.flatten_cell(geom)
    cls = rg.registry_class(at, "reconstruction-related")
    top = rg.top_atoms(at)
    c = collections.Counter(int(x) for x in cls[top])
    n = int(top.sum())
    if sid == "R1-hcp-terminated":
        ok = c.get(2, 0) >= 0.8 * n
        return ok, f"hcp-like {c.get(2,0)}/{n} surface atoms after relaxation"
    if sid == "R2-stripe-wall":
        ok = c.get(0, 0) >= 2 and c.get(2, 0) >= 2 and c.get(1, 0) >= 2
        return ok, f"fcc-like {c.get(0,0)} / transition {c.get(1,0)} / hcp-like {c.get(2,0)} of {n}"
    return True, ""


def noneq_descriptor(s, parent_geom, geom):
    """Small, interpretable descriptor of what a non-equilibrium configuration did to its parent."""
    pa, pb = parent_geom.get_positions(), geom.get_positions()
    disp = np.zeros(len(pb))
    if len(pa) == len(pb):
        cell = geom.get_cell().array; C = np.array([cell[0][:2], cell[1][:2]]); Ci = np.linalg.inv(C)
        d = pb - pa; f = d[:, :2] @ Ci; f -= np.round(f); d[:, :2] = f @ C
        disp = np.linalg.norm(d, axis=1)
    cn = coordination(geom)
    n_under = int((cn <= 8).sum())
    # closest Au-Au contact, in-plane periodic
    p = geom.get_positions(); cell = geom.get_cell().array; best = np.inf
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = si * cell[0] + sj * cell[1]
            dd = np.linalg.norm(p[:, None, :] - (p[None, :, :] + sh), axis=-1)
            dd = np.where(dd > 0.1, dd, np.inf); best = min(best, float(dd.min()))
    kind = s["config"].split("_")[0] if s["config"].startswith(("coll", "path")) else s["config"][:4]
    return dict(kind=kind, family=s["family"], max_disp_A=float(disp.max()), mean_disp_A=float(disp.mean()),
                n_moved_gt_0p3=int((disp > 0.3).sum()), min_contact_A=best, n_under_coordinated=n_under,
                n_atoms=len(geom))


def farthest_point(cands, k, keys, families_required):
    """Greedy farthest-point selection on z-scored descriptors, forced to cover each family at least once."""
    if not cands: return []
    X = np.array([[c["desc"][key] for key in keys] for c in cands], float)
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    chosen = []
    # cover families first with the most displaced member of each
    for fam in families_required:
        idx = [i for i, c in enumerate(cands) if c["family"] == fam and i not in chosen]
        if idx: chosen.append(max(idx, key=lambda i: cands[i]["desc"]["max_disp_A"]))
    while len(chosen) < min(k, len(cands)):
        dmin = np.full(len(cands), np.inf)
        for j in chosen: dmin = np.minimum(dmin, np.linalg.norm(X - X[j], axis=1))
        dmin[chosen] = -1
        chosen.append(int(np.argmax(dmin)))
    return chosen


# ----------------------------------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True)
    S = list(json.load(open(f"{ROOT}/dataset_v1/states.json"))["states"].values())
    avail = field_completeness(S)
    with open(f"{OUT}/label_availability.jsonl", "w") as f:
        for r in avail: f.write(json.dumps(r) + "\n")
    complete = {r["state_id"] for r in avail if r["all_fields_complete"] and r["mu_within_tolerance"]}
    by_geom = collections.defaultdict(list)
    for s in S: by_geom[s["geometry_id"]].append(s)
    geoms = {}
    for gid, ss in by_geom.items():
        s0 = ss[0]
        cfg = "relaxed" if s0["config"] in ("relax", "relaxed") else s0["config"]
        geoms[gid] = dict(gid=gid, sid=s0["structure_id"], config=cfg, family=s0["family"],
                          states=sorted(ss, key=lambda s: s["electronic_state"]["TARGETMU_eV"]),
                          relax_ok=any((s.get("labels") or {}).get("relaxation_reached_accuracy") for s in ss))

    selected, excluded = [], []

    def take(g, states, role, reason):
        for s in states:
            e = s["electronic_state"]
            selected.append(dict(state_id=s["state_id"], geometry_id=g["gid"], structure_id=g["sid"],
                                 config=g["config"], family=g["family"], role=role,
                                 U_V=round(e["mu_reference_eV"] - e["mu_e_actual_eV"], 4),
                                 TARGETMU_eV=e["TARGETMU_eV"], reason=reason))

    def drop(g, reason, states=None):
        for s in (states or g["states"]):
            excluded.append(dict(state_id=s["state_id"], geometry_id=g["gid"], structure_id=g["sid"],
                                 config=g["config"], family=g["family"], reason=reason))

    # ---- 1. base references: one geometry per structure, relaxed preferred, environment checked
    base_geoms = {}
    for sid in BASE16:
        cands = [g for g in geoms.values() if g["sid"] == sid and g["config"] in ("ideal", "relaxed", "relax")]
        rel = [g for g in cands if g["config"] in ("relaxed", "relax")]; ide = [g for g in cands if g["config"] == "ideal"]
        pick, why = None, ""
        if rel and rel[0]["relax_ok"]:
            geom = load_geom(rel[0]["states"][0]); ok, note = registry_kept(sid, geom)
            if ok:
                pick, why = rel[0], "accepted relaxation" + (f"; {note}" if note else "")
            else:
                pick, why = ide[0], f"relaxation discarded the target environment ({note}); ideal geometry kept"
                drop(rel[0], f"relaxation lost the environment this structure exists for: {note}")
        elif ide:
            pick, why = ide[0], "no accepted relaxation; ideal geometry"
        if pick is None:
            continue
        base_geoms[sid] = (pick, why)
    # duplicate check among the chosen base geometries
    items = list(base_geoms.items()); dup_notes = {}
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            (sa, (ga, _)), (sb, (gb, _)) = items[i], items[j]
            a = load_geom(ga["states"][0]); b = load_geom(gb["states"][0])
            if len(a) != len(b) or not np.allclose(a.get_cell().array, b.get_cell().array, atol=1e-3): continue
            r = rmsd_translated(a, b)
            if r < RMSD_DUP: dup_notes[sb] = f"same geometry as {sa} after relaxation (RMSD {r:.3f} A)"
    for sid, (g, why) in base_geoms.items():
        if sid in dup_notes:
            drop(g, dup_notes[sid]); continue
        ok = [s for s in g["states"] if s["state_id"] in complete]
        bad = [s for s in g["states"] if s["state_id"] not in complete]
        take(g, ok, "base_reference", f"{why}; all available potentials")
        if bad: drop(g, "incomplete labels (see label_availability.jsonl)", bad)

    # ---- 2. non-equilibrium: farthest-point on a descriptor, under a family cover
    base_lookup = {}
    for s in S:
        if s["config"] in ("relaxed", "relax"): base_lookup[s["structure_id"]] = s
    cands = []
    for g in geoms.values():
        if not g["config"].startswith(("pert", "coll", "path")): continue
        parent = base_lookup.get(g["sid"])
        if parent is None: continue
        geom = load_geom(g["states"][0]); pgeom = load_geom(parent)
        cands.append(dict(g=g, family=g["family"], desc=noneq_descriptor(g["states"][0], pgeom, geom)))
    keys = ["max_disp_A", "mean_disp_A", "n_moved_gt_0p3", "min_contact_A", "n_under_coordinated", "n_atoms"]
    fams = sorted({c["family"] for c in cands})
    pert = [c for c in cands if c["desc"]["kind"] in ("pert", "coll")]
    path = [c for c in cands if c["desc"]["kind"] == "path"]
    for pool, k, tag in ((pert, N_NONEQ_PERT, "perturbation/collective"), (path, N_NONEQ_PATH, "path image")):
        idx = farthest_point(pool, k, keys, [f for f in fams if any(c["family"] == f for c in pool)])
        chosen = {id(pool[i]["g"]) for i in idx}
        for c in pool:
            g = c["g"]; d = c["desc"]
            if id(g) in chosen:
                outer = [s for s in g["states"] if round(s["electronic_state"]["TARGETMU_eV"], 4) in (-5.4071, -4.4071)
                         and s["state_id"] in complete]
                if len(outer) < 2:
                    outer = sorted([s for s in g["states"] if s["state_id"] in complete],
                                   key=lambda s: -abs(s["electronic_state"]["TARGETMU_eV"] - MU0))[:2]
                why = (f"{tag}, farthest-point pick: max displacement {d['max_disp_A']:.2f} A, "
                       f"{d['n_moved_gt_0p3']} atoms moved > 0.3 A, closest contact {d['min_contact_A']:.2f} A, "
                       f"{d['n_under_coordinated']} under-coordinated; two outermost complete potentials")
                take(g, outer, "non_equilibrium", why)
                drop(g, "non-equilibrium geometry kept at two potentials only", [s for s in g["states"] if s not in outer])
            else:
                drop(g, f"{tag} not selected by farthest-point sampling (descriptor too close to a chosen one)")

    # ---- 3. supplementary references: one potential each, alternating sign to balance the set
    sign = 1
    for sid in SUPP12:
        cands = [g for g in geoms.values() if g["sid"] == sid and g["config"] in ("ideal", "relaxed")]
        if not cands: continue
        rel = [g for g in cands if g["config"] == "relaxed" and g["relax_ok"]]
        g = rel[0] if rel else [g for g in cands if g["config"] == "ideal"][0]
        okst = [s for s in g["states"] if s["state_id"] in complete]
        if not okst:
            drop(g, "no state with complete labels"); continue
        want = -5.4071 if sign > 0 else -4.4071
        pick = min(okst, key=lambda s: abs(s["electronic_state"]["TARGETMU_eV"] - want))
        take(g, [pick], "supplementary_reference",
             f"size / spacing / environment reference, one potential (U = {MU0 - pick['electronic_state']['mu_e_actual_eV']:+.2f} V, "
             f"alternating sign across the twelve)")
        drop(g, "supplementary reference kept at one potential only", [s for s in g["states"] if s is not pick])
        sign = -sign
    # periodic copies and everything else not touched above
    taken_geoms = {r["geometry_id"] for r in selected} | {r["geometry_id"] for r in excluded}
    for g in geoms.values():
        if g["gid"] in taken_geoms: continue
        if g["sid"] in PERIODIC_COPIES:
            drop(g, f"strictly periodic copy of {PERIODIC_COPIES[g['sid']]}: same local environment in a larger cell")
        elif g["config"] in ("ideal", "relaxed"):
            drop(g, "ideal/relaxed geometry of a base structure whose other configuration was chosen as the reference")
        else:
            drop(g, "not in the selection budget")

    # ---- write
    with open(f"{OUT}/old_selected_140.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(selected[0].keys())); w.writeheader(); w.writerows(selected)
    with open(f"{OUT}/old_excluded.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(excluded[0].keys())); w.writeheader(); w.writerows(excluded)

    roles = collections.Counter(r["role"] for r in selected)
    gsel = collections.defaultdict(set)
    for r in selected: gsel[r["role"]].add(r["geometry_id"])
    pots = collections.Counter(r["TARGETMU_eV"] for r in selected)
    n_inc = sum(1 for r in avail if not r["all_fields_complete"])
    L = ["# Old-state selection for rough_sampling_v1", "",
         f"From {len(S)} states / {len(geoms)} geometries. Selected **{len(selected)} states** over "
         f"{sum(len(v) for v in gsel.values())} geometries; {len(excluded)} state rows excluded with reasons.", "",
         "| role | geometries | states |", "|---|---|---|"]
    for k in ("base_reference", "non_equilibrium", "supplementary_reference"):
        L.append(f"| {k} | {len(gsel[k])} | {roles[k]} |")
    L += ["", "Potentials in the selected set: " + ", ".join(f"{k}: {v}" for k, v in sorted(pots.items())), "",
          f"States with any incomplete field label: {n_inc} (listed in label_availability.jsonl).", ""]
    if dup_notes: L += ["Base geometries dropped as duplicates:"] + [f"- {k}: {v}" for k, v in dup_notes.items()] + [""]
    L += ["Base-reference choices:"] + [f"- {sid}: {why}" for sid, (g, why) in base_geoms.items()]
    open(f"{OUT}/selection_report.md", "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
