#!/usr/bin/env python3
"""mu05 extension (2026-09-28): add the two OUTER electron chemical potentials to every distinct geometry of the frozen
dataset -- fixed-geometry CP-DFT single points with the source state's production configuration.

    mu_target = mu_0 - dU,  mu_0 = -4.9071 eV:   dU = +0.5 V -> TARGETMU = -5.4071 eV ;  dU = -0.5 V -> TARGETMU = -4.4071 eV
(internal reference, NOT V vs RHE, NOT a per-structure PZC; existing -5.1071/-4.9071/-4.7071 states untouched; no middle points)

Task objects are DISTINCT GEOMETRIES, not queue records:
  * ideal / pert / coll / path states : the POSCAR that was actually computed (never regenerated);
  * relaxed geometry                  : the accepted relaxation's CONTCAR with the relaxation's constraint flags; the relax task and
                                        its reused 'relaxed@-4.9071' record are ONE geometry (extension added once);
  * geometry_id = <structure>__<config>__<sha1 of cell+positions>; dedupe = geometry_id + config_version + TARGETMU (task_id);
    every source record of the same geometry is checked to agree within 1e-3 A; idempotent re-runs add nothing;
  * retired mislabeled structures and unaccepted relaxations never enter (they are not in the queue).
Numerics copied from the source: KPOINTS verbatim, POTCAR, production INCAR (IBRION=-1, NSW=0, window/DIPOL/solvent unchanged; the
pilot inputs of the reused states were verified identical to the production template), cell and positions verbatim (strained
configurations keep their own cell). Warm start: ICHARG=1 with a COPY of the CHGCAR of the completed same-side neighbour
(-5.1071 for -5.4071, -4.7071 for -4.4071; else -4.9071); never ICHARG=11/12; cold start when no neighbour is complete.
Priority 80 (85 for > 200-atom references) -> the farms finish the original plan first. campaign = mu05_extension.
Labels are the ACTUAL converged mu_e / N_e (evaluate() accepts |mu_e - target| <= 0.011 as for all other states).

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/extend_mu05.py --dry-run | (no flag: create tasks) | --report
"""
import hashlib
import json
import os
import re
import shutil
import sys
import time

import numpy as np
from ase.io import read

sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/scripts")
import production as P  # noqa: E402

MU0 = -4.9071
NEW = {"-5.4071": ("+0.5", "-5.1071"), "-4.4071": ("-0.5", "-4.7071")}      # target -> (dU, preferred warm-start neighbour)
OLD = ["-5.1071", "-4.9071", "-4.7071"]
CONFIG_VERSION = "production-K.8-2026-09-26"
CAMPAIGN = "mu05_extension"
SAME_GEOM_TOL = 0.005      # A; source records of one geometry must agree to this (measured: pilot-reused ideal states 0.003 A)
MANIFEST = f"{P.PROD}/mu05_extension_manifest.json"
REPORT = f"{P.PROD}/mu05_extension_manifest.md"
MODE = "dry" if "--dry-run" in sys.argv else "report" if "--report" in sys.argv else "refresh" if "--refresh" in sys.argv else \
       "seed-dry" if "--seed-dry-run" in sys.argv else "create"
# --refresh (run any time, under the queue lock): (a) pending cold-start extension tasks whose same-side neighbour has since
# completed get the warm start (CHGCAR copy + ICHARG=1) -- inputs of PENDING tasks of this campaign only; (b) metadata of
# extension records (source_geometry / note) is re-derived from the index; (c) manifest json + report rewritten.


def geom_hash(atoms):
    h = hashlib.sha1()
    h.update(np.round(atoms.get_cell().array, 3).tobytes()); h.update(np.round(atoms.get_positions(), 3).tobytes())
    return h.hexdigest()[:12]


def read_flags(poscar):
    lines = open(poscar).read().splitlines(); n = int(lines[6].split()[0])
    return np.array([l.split()[3] == "T" for l in lines[9:9 + n]])


def task_dir(t):
    d = t["dir"].split(" ")[0]
    return d if d.startswith("/") else f"{P.ROOT}/{d}"


def source_geometry(t):
    """(Atoms actually computed, movable mask, geometry file, run dir) for a queue record; None if not usable."""
    d = task_dir(t)
    if os.path.basename(d).startswith("relax__"):             # relax task, or its reused relaxed@-4.9071 record
        if not os.path.exists(f"{d}/CONTCAR") or not any(x["task_id"] == f"{t['structure_id']}__relax__mu{P.MU_RELAX}" and
                                                          x["status"] == "completed" for x in Q):
            return None                                       # relaxation not accepted -> no final geometry yet
        return read(f"{d}/CONTCAR"), read_flags(f"{d}/POSCAR"), f"{d}/CONTCAR", d
    if not os.path.exists(f"{d}/POSCAR"): return None
    assert "Selective" in open(f"{d}/POSCAR").read(), d
    return read(f"{d}/POSCAR"), read_flags(f"{d}/POSCAR"), f"{d}/POSCAR", d


def build_index(q):
    """geometry index: (structure_id, config) -> dict(...); config 'relax' folded into 'relaxed'."""
    geoms = {}
    base = [t for t in q if t.get("campaign") != CAMPAIGN]
    base.sort(key=lambda t: (t["dir"].startswith("03_pilot"), t["kind"] != "relax", t["task_id"]))   # canonical source first
    for t in base:
        src = source_geometry(t)
        if src is None: continue
        at, flags, gfile, d = src
        cfg = "relaxed" if t["config"] in ("relax", "relaxed") else t["config"]
        key = (t["structure_id"], cfg)
        if key not in geoms:
            geoms[key] = dict(structure_id=t["structure_id"], family=t["family"], tier=t["tier"], config=cfg, n_atoms=len(at),
                              atoms=at, flags=flags, source_file=gfile, source_task=t["task_id"], kpoints_file=f"{d}/KPOINTS",
                              geometry_id=f"{t['structure_id']}__{cfg}__{geom_hash(at)}", mus={})
        else:
            g = geoms[key]
            assert len(at) == g["n_atoms"] and np.allclose(at.get_cell().array, g["atoms"].get_cell().array, atol=1e-4), \
                f"{t['task_id']}: cell/N differs from {g['source_task']} -> not one geometry, refuse"
            dmax = float(np.linalg.norm(at.get_positions() - g["atoms"].get_positions(), axis=1).max())
            # 0.003 A: the three pilot-reused ideal states (Step-8x1/Step-16x1) differ from the production POSCAR by a rounding
            # of the step-edge top atoms; recorded, treated as the same geometry (the production POSCAR is the source).
            assert dmax <= SAME_GEOM_TOL, f"{t['task_id']}: positions differ from {g['source_task']} by {dmax:.4f} A -> refuse"
            if dmax > 1e-4:
                g.setdefault("geometry_note", []).append(f"{t['task_id']} ({os.path.relpath(gfile, P.ROOT)}) deviates by max {dmax:.4f} A")
            assert (flags == g["flags"]).all(), f"{t['task_id']}: constraint mask differs from {g['source_task']}"
        geoms[key]["mus"][t["mu"]] = dict(task_id=t["task_id"], status=t["status"], dir=d, kind=t["kind"],
                                          elapsed_farm_min=t.get("elapsed_farm_min"), result=t.get("result"))
    for t in q:                                              # extension records register only their potential
        if t.get("campaign") == CAMPAIGN:
            key = (t["structure_id"], t["config"])
            if key in geoms:
                geoms[key]["mus"][t["mu"]] = dict(task_id=t["task_id"], status=t["status"], dir=task_dir(t), kind=t["kind"],
                                                  elapsed_farm_min=t.get("elapsed_farm_min"), result=t.get("result"))
    return geoms


def est_minutes(g, q):
    """measured mean farm time of completed single points of the same structure; else calibrated model c(N) x 0.81."""
    same = [t["elapsed_farm_min"] for t in q if t["status"] == "completed" and t["kind"] == "sp" and t.get("elapsed_farm_min")
            and t["structure_id"] == g["structure_id"] and t.get("campaign") != CAMPAIGN]
    n = g["n_atoms"]
    return (float(np.mean(same)), "measured") if same else (0.55 * n * max(1, n / 72) ** 0.5 * 0.81, "model")


def plan_extension(geoms, q):
    plan, present = [], []
    for key, g in sorted(geoms.items()):
        for mu, (dU, pref) in NEW.items():
            task_id = f"{g['structure_id']}__{g['config']}__mu{mu}"
            if mu in g["mus"] or any(t["task_id"] == task_id for t in q):
                present.append(task_id); continue
            warm = None
            for cand in (pref, "-4.9071"):
                s = g["mus"].get(cand)
                if s and s["status"] in ("completed", "reused") and os.path.exists(f"{s['dir']}/CHGCAR") and \
                        os.path.getsize(f"{s['dir']}/CHGCAR") > 1e6:
                    warm = (cand, f"{s['dir']}/CHGCAR"); break
            e, how = est_minutes(g, q)
            plan.append(dict(task_id=task_id, key=key, g=g, mu=mu, dU=dU, warm=warm, est_min=e, est_how=how))
    return plan, present


def res_ne_mu(s):
    """(N_e, mu_e) actually converged, from the farm's result string; None if the state has none."""
    m = re.search(r"N_e=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)", str(s.get("result") or ""))
    return (float(m.group(1)), float(m.group(2))) if m else None


def nelect_model(geoms, q):
    """Per-geometry linear N_e(mu_e) from its own CONVERGED states, plus an area-scaled fallback slope.

    Validated 2026-09-29 on the 144 +-0.5 V states already computed: predicting N_e at the target from the three
    original potentials is accurate to 0.008 e (median) / 0.018 e (max), while the walk from the neutral count is
    0.41 e (median). Seeding NELECT therefore starts the CP loop ~50x closer than the default neutral start, which is
    where the large cells spend their time (every CP round otherwise restarts from neutral; ICHARG=1 does not help,
    VASP rescales the CHGCAR it reads to the current NELECT). NELECT is only the STARTING count -- the CP loop still
    decides the converged N_e from TARGETMU, and labels remain the actual converged mu_e / N_e."""
    fits, per_area = {}, []
    for key, g in geoms.items():
        pts = [res_ne_mu(s) for s in g["mus"].values() if s["status"] in ("completed", "reused")]
        pts = [p for p in pts if p]
        area = float(np.linalg.norm(np.cross(g["atoms"].get_cell()[0], g["atoms"].get_cell()[1])))
        if len(pts) >= 2:
            slope, icept = np.polyfit([p[1] for p in pts], [p[0] for p in pts], 1)
            fits[key] = dict(slope=float(slope), icept=float(icept), n=len(pts), area=area, basis="own states")
            per_area.append(slope / area)
        elif len(pts) == 1:
            fits[key] = dict(pt=pts[0], area=area, n=1, basis="single state + area-scaled slope")
    med = float(np.median(per_area)) if per_area else None
    for key, f in fits.items():
        if f["n"] == 1 and med is not None:                     # capacitance scales with the cell area
            f["slope"] = med * f["area"]; f["icept"] = f["pt"][0] - f["slope"] * f["pt"][1]
    return fits, med


def seed_nelect(q, geoms, dry=False):
    """Write NELECT (predicted starting electron count) into the INCAR of PENDING extension tasks."""
    fits, med = nelect_model(geoms, q)
    n_set = 0
    for t in q:
        if t.get("campaign") != CAMPAIGN or t["status"] != "pending": continue
        d = task_dir(t)
        if os.path.exists(f"{d}/log.out"): continue                     # a run already touched this directory
        f = fits.get((t["structure_id"], t["config"]))
        if not f or "slope" not in f: continue
        pred = f["slope"] * float(t["mu"]) + f["icept"]
        if abs(pred - t.get("nelect_seed", -1)) < 1e-4: continue        # already seeded with this value
        if dry:
            print(f"  would seed {t['task_id']}: NELECT={pred:.4f} (neutral {11*t['n_atoms']}, {f['basis']}, slope {f['slope']:.3f} e/eV)")
            n_set += 1; continue
        inc = [l for l in open(f"{d}/INCAR").read().splitlines() if not l.startswith("NELECT")]
        i = next((j for j, l in enumerate(inc) if l.startswith("LCEP")), len(inc))
        inc.insert(i, f"NELECT = {pred:.4f}   # predicted start for TARGETMU={t['mu']} from this geometry's converged "
                      f"N_e(mu_e) ({f['basis']}, slope {f['slope']:.4f} e/eV); the CP loop still decides the final N_e")
        open(f"{d}/INCAR", "w").write("\n".join(inc) + "\n")
        t["nelect_seed"] = round(pred, 4); t["nelect_neutral"] = 11 * t["n_atoms"]; t["nelect_basis"] = f["basis"]
        n_set += 1
    return n_set


def res_mu(s):
    """actual mu_e from the farm's result string ('... mu_e=-5.107726'), else None."""
    import re
    m = re.search(r"mu_e=\s*([-\d.]+)", str(s.get("result") or ""))
    return float(m.group(1)) if m else None


def report(q, geoms):
    ext = [t for t in q if t.get("campaign") == CAMPAIGN]
    st = {}
    for t in ext: st[t["status"]] = st.get(t["status"], 0) + 1
    done = [t for t in ext if t["status"] == "completed"]
    act = sum(t.get("elapsed_farm_min") or 0 for t in done) / 60
    est = sum(t.get("est_min") or 0 for t in ext) / 60
    L = [f"# mu05 extension: +-0.5 V end points for every distinct geometry", "",
         f"Generated {time.strftime('%Y-%m-%d %H:%M')} by `scripts/extend_mu05.py --report`. mu_target = mu0 - dU, mu0 = {MU0} eV "
         f"(internal reference). TARGETMU -5.4071 (dU = +0.5 V) and -4.4071 (dU = -0.5 V); existing -5.1071/-4.9071/-4.7071 kept; "
         f"no intermediate points added. Static CP single points (IBRION=-1, NSW=0) on the geometry actually computed "
         f"(relaxed = accepted CONTCAR). Warm start = ICHARG=1 from a copy of the same-side neighbour's CHGCAR.", "",
         f"Distinct geometries: {len(geoms)}. Extension tasks: {len(ext)} -> status {st}. Estimated {est:.0f} node-h "
         f"(measured same-structure single-point means); actual so far {act:.0f} node-h over {len(done)} completed.", "",
         "Status legend: C completed/accepted, R running, P pending, U reused, F failed, - not planned. `w:` warm start source.", "",
         "| structure | config | N | family | -5.1071 | -4.9071 | -4.7071 | -5.4071 (+0.5 V) | -4.4071 (-0.5 V) | est min | actual min |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    code = {"completed": "C", "running": "R", "pending": "P", "reused": "U", "failed": "F"}
    fam_cov = {}
    for key, g in sorted(geoms.items()):
        cells = []
        for mu in OLD:
            s = g["mus"].get(mu); cells.append(code.get(s["status"], "?") if s else "-")
        est_c, act_c = [], []
        for mu in NEW:
            s = g["mus"].get(mu)
            if not s: cells.append("-"); continue
            t = next(x for x in ext if x["task_id"] == s["task_id"])
            w = ("w:" + str(t.get("warm_start_mu"))) if t.get("warm_start", "").startswith("ICHARG") else "cold"
            mu_a = res_mu(s)
            cells.append(code.get(s["status"], "?") + f" ({w})" + (f" mu_e={mu_a:.4f}" if mu_a is not None else ""))
            est_c.append(t.get("est_min") or 0); act_c.append(t.get("elapsed_farm_min") or 0)
            fc = fam_cov.setdefault(g["family"], [0, 0]); fc[0] += 1; fc[1] += s["status"] == "completed"
        L.append(f"| {g['structure_id']} | {g['config']} | {g['n_atoms']} | {g['family']} | " + " | ".join(cells) +
                 f" | {sum(est_c):.0f} | {sum(act_c):.0f} |")
    L += ["", "## Per-family coverage of the new end points (completed / planned)", ""]
    L += [f"- {f}: {c}/{n}" for f, (n, c) in sorted(fam_cov.items())]
    open(REPORT, "w").write("\n".join(L) + "\n")
    print(f"report: {REPORT}  ({len(ext)} extension tasks, {st})")


def write_manifest(q, geoms, plan=None, made=None, est=None, present=None):
    old = json.load(open(MANIFEST)) if os.path.exists(MANIFEST) else {}
    m = dict(created=old.get("created", time.strftime("%Y-%m-%d %H:%M")), updated=time.strftime("%Y-%m-%d %H:%M"), mu0_eV=MU0,
             targets={k: dict(dU_V=v[0], warm_pref=v[1]) for k, v in NEW.items()}, config_version=CONFIG_VERSION, n_geometries=len(geoms),
             n_tasks_created=made if made is not None else old.get("n_tasks_created"),
             n_already_present=len(present) if present is not None else old.get("n_already_present"),
             est_node_h=round(est, 1) if est is not None else old.get("est_node_h"),
             geometries=[dict(geometry_id=g["geometry_id"], structure_id=g["structure_id"], config=g["config"], family=g["family"],
                              n_atoms=g["n_atoms"], source=os.path.relpath(g["source_file"], P.ROOT), geometry_note=g.get("geometry_note", []),
                              states=[dict(mu=mu, task_id=s["task_id"], status=s["status"],
                                           mu_e_actual=res_mu(s),
                                           wall_min=s["elapsed_farm_min"]) for mu, s in sorted(g["mus"].items())])
                         for _, g in sorted(geoms.items())])
    if plan is not None:
        m["tasks"] = [dict(task_id=p["task_id"], geometry_id=p["g"]["geometry_id"], mu_target=p["mu"], dU=p["dU"],
                           warm_start=(p["warm"][0] if p["warm"] else None), est_min=round(p["est_min"], 1), est_how=p["est_how"]) for p in plan]
    else:
        m["tasks"] = old.get("tasks", [])
    json.dump(m, open(MANIFEST, "w"), indent=1)


def apply_warm(t, src_mu, src_chgcar):
    d = task_dir(t)
    shutil.copy(src_chgcar, f"{d}/CHGCAR")                                          # a COPY (VASP rewrites CHGCAR at the end)
    inc = open(f"{d}/INCAR").read()
    assert "IBRION = -1" in inc and "ICHARG" not in inc, d
    open(f"{d}/INCAR", "w").write(inc.replace("IBRION = -1", "ICHARG = 1\nIBRION = -1"))
    t["warm_start"] = f"ICHARG=1 from {src_chgcar}"; t["warm_start_mu"] = src_mu


Q = P.load_queue()
if MODE == "report":
    report(Q, build_index(Q)); sys.exit(0)

if MODE == "seed-dry":
    geoms = build_index(Q); n = seed_nelect(Q, geoms, dry=True)
    print(f"would seed NELECT on {n} pending tasks"); sys.exit(0)

if MODE == "refresh":
    with P.queue_lock():
        q = P.load_queue(); Q = q
        geoms = build_index(q)
        n_warm, n_meta = 0, 0
        for t in q:
            if t.get("campaign") != CAMPAIGN: continue
            g = geoms.get((t["structure_id"], t["config"]))
            if g is None: continue
            if t["source_geometry"] != g["source_file"]:
                t["source_geometry"] = g["source_file"]; t["geometry_id"] = g["geometry_id"]
                t["note"] = f"{CAMPAIGN}: dU={'+0.5' if t['mu']=='-5.4071' else '-0.5'} V of mu0={MU0}; geometry = {os.path.relpath(g['source_file'], P.ROOT)}"
                n_meta += 1
            if t["warm_start"].startswith("ICHARG") and "warm_start_mu" not in t:      # backfill for records created before this field
                src_dir = os.path.dirname(t["warm_start"].split("from ")[1])
                t["warm_start_mu"] = next((mu for mu, s in g["mus"].items() if s["dir"] == src_dir), None); n_meta += 1
            if t["status"] == "pending" and t["warm_start"] == "cold start" and not os.path.exists(f"{task_dir(t)}/log.out"):
                for cand in (NEW[t["mu"]][1], "-4.9071"):
                    s = g["mus"].get(cand)
                    if s and s["status"] in ("completed", "reused") and os.path.exists(f"{s['dir']}/CHGCAR") and os.path.getsize(f"{s['dir']}/CHGCAR") > 1e6:
                        apply_warm(t, cand, f"{s['dir']}/CHGCAR"); n_warm += 1; break
        n_seed = seed_nelect(q, geoms)
        P.save_queue(q)
        if n_warm or n_meta or n_seed:
            P.log(f"[extend_mu05] refresh: {n_warm} pending cold-start tasks given a warm start, {n_meta} records re-pointed "
                  f"to their geometry source, {n_seed} given a predicted NELECT start")
        write_manifest(q, geoms)
        report(q, geoms)
        print(f"refresh: warm starts applied {n_warm}, metadata updated {n_meta}, NELECT seeded {n_seed}")
    sys.exit(0)

rows = P.load_plan()
with P.queue_lock():
    q = P.load_queue(); Q = q
    geoms = build_index(q)
    plan, present = plan_extension(geoms, q)
    est = sum(p["est_min"] for p in plan) / 60
    n_meas = sum(p["est_how"] == "measured" for p in plan)
    print(f"distinct geometries: {len(geoms)}  ({sum(1 for g in geoms.values() if g['config']=='relaxed')} relaxed, "
          f"{sum(1 for g in geoms.values() if g['config']=='ideal')} ideal, {len(geoms)-sum(1 for g in geoms.values() if g['config'] in ('relaxed','ideal'))} pert/coll/path)")
    print(f"extension tasks to add: {len(plan)}   already present: {len(present)}   warm starts available: {sum(1 for p in plan if p['warm'])}")
    print(f"estimated cost: {est:.0f} node-h ({n_meas} of {len(plan)} from measured same-structure times, rest from model c(N)x0.81)")
    by = {}
    for p in plan:
        b = by.setdefault(p["g"]["family"], [0, 0.0]); b[0] += 1; b[1] += p["est_min"] / 60
    for fam, (n, h) in sorted(by.items(), key=lambda kv: -kv[1][1]): print(f"   {fam:28s} {n:4d} tasks  {h:6.1f} node-h")
    big = [p for p in plan if p["g"]["n_atoms"] > 200]
    print(f"   of which > 200-atom references: {len(big)} tasks, {sum(p['est_min'] for p in big)/60:.0f} node-h (priority 85)")
    if MODE == "dry":
        for p in plan[:8]:
            print("  e.g.", p["task_id"], "<-", os.path.relpath(p["g"]["source_file"], P.ROOT), "| warm:", p["warm"][0] if p["warm"] else "cold",
                  f"| est {p['est_min']:.0f} min ({p['est_how']})")
        sys.exit(0)

    made = 0
    for p in plan:
        g = p["g"]; row = rows[g["structure_id"]]
        prio = 85 if g["n_atoms"] > 200 else 80
        t = P.make_task(q, g["structure_id"], row, g["config"], p["mu"], "sp", g["atoms"], g["flags"], priority=prio,
                        note=f"{CAMPAIGN}: dU={p['dU']} V of mu0={MU0}; geometry = {os.path.relpath(g['source_file'], P.ROOT)}")
        if t is None: continue
        d = t["dir"]
        shutil.copy(g["kpoints_file"], f"{d}/KPOINTS")                                   # source k-mesh verbatim
        t.update(campaign=CAMPAIGN, geometry_id=g["geometry_id"], config_type=g["config"], source_geometry=g["source_file"],
                 config_version=CONFIG_VERSION, mu_reference_eV=MU0, delta_U_target_V=float(p["dU"]), mu_target_eV=float(p["mu"]),
                 warm_start="cold start", est_min=round(p["est_min"], 1))
        if p["warm"]: apply_warm(t, p["warm"][0], p["warm"][1])
        made += 1
    q.sort(key=lambda t: (t["priority"], t["n_atoms"], t["task_id"]))
    P.save_queue(q)
    P.log(f"[extend_mu05] {made} extension tasks queued (priority 80/85, campaign {CAMPAIGN}, est {est:.0f} node-h); queue now {len(q)} entries")
    geoms = build_index(q)
    write_manifest(q, geoms, plan=plan, made=made, est=est, present=present)
    print(f"{made} tasks created; manifest: {MANIFEST}")
    report(q, geoms)
