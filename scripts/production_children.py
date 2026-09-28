#!/usr/bin/env python3
"""Children of an accepted relaxation (dataset_plan_v1 rev 2): called by production.py feed as
    production_children.py <structure_id>
right after <structure_id>/relax__mu-4.9071 is accepted. Generates, from the relaxed geometry (CONTCAR):

  relaxed reference   : the relaxation's own final step is recorded as the relaxed state at mu = -4.9071 (its final
                        SCF, forces and grids belong to the final geometry); single points on CONTCAR at the other mu.
  perturbations (2)   : Gaussian displacements of the movable atoms, sigma = 0.05 A (seed 1) and 0.10 A (seed 3);
                        min Au-Au distance >= 2.6 A enforced (next seed otherwise, recorded).
  collective (n_coll) : (a) top plane moved down by 3 % of d111 = 0.072 A; (b) in-plane strain +1 % (cell and xy scaled)
                        -- for strip steps (b) is replaced by an edge-row bend of +-0.15 A across the step, alternating
                        between the ny rows.
  paths (n_path)      : A1-fcc: 1 image, adatom at the bridge point between the relaxed fcc position and the nearest hcp
                        hollow (xy of the nearest layer-2 atom), same z.  Step-8x2: 3 images at 0.25/0.5/0.75 between
                        the relaxed strip (CONTCAR order) and the ideal Step-8x2_edge-vacancy_plus_foot-adatom geometry
                        (same atom order; the moving atom's in-plane path uses the minimum image). Other families: none in Batch A.
  potentials          : all three for N <= 100 (Batch A), the two end points otherwise.
Every child is a static single point; nothing is relaxed again. Tasks are appended to 05_production/queue.json
with priority 25 (after main relaxations/ideal points, before reference points).
"""
import json
import os
import sys

import numpy as np
from ase.io import read

sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/scripts")
import production as P  # noqa: E402

A0 = 4.158
D111 = A0 / np.sqrt(3)
SEEDS = {"pert05": (1, 0.05), "pert10": (3, 0.10)}
LOCK = f"{P.PROD}/queue.lock"


def mic(d, cell):
    """shortest in-plane image of the Cartesian difference d. For a strongly sheared cell the
    component-wise fractional rounding is NOT the shortest vector, so all 9 neighbouring images are tried."""
    C = cell[:2, :2]; f = np.linalg.solve(C.T, d[:2]); f -= np.round(f)
    best = None
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            v = (f + np.array([i, j])) @ C
            if best is None or np.linalg.norm(v) < np.linalg.norm(best): best = v
    return np.array([*best, d[2]])


def min_dist_xy_periodic(at):
    a = at.copy(); a.set_pbc((True, True, False))
    from ase.neighborlist import neighbor_list
    i, j, d = neighbor_list("ijd", a, 3.5)
    return float(d.min()) if len(d) else 9.9


def read_flags(poscar):
    lines = open(poscar).read().splitlines(); n = int(lines[6].split()[0])
    return np.array([l.split()[3] == "T" for l in lines[9:9 + n]])


def main(sid, dry=False, contcar_override=None):
    rows = P.load_plan(); row = rows[sid]
    rdir = f"{P.PROD}/{sid}/relax__mu{P.MU_RELAX}"
    contcar = contcar_override or f"{rdir}/CONTCAR"
    assert os.path.exists(contcar), f"{sid}: no CONTCAR"
    rel = read(contcar); movable = read_flags(f"{rdir}/POSCAR")
    assert len(rel) == row["n_atoms"] and movable.sum() > 0
    pos0 = rel.get_positions(); cell = rel.get_cell().array
    mus = row["mu"] if row["n_atoms"] <= 100 else ["-5.1071", "-4.7071"]
    children = {}                     # config name -> (Atoms, note)

    # perturbations
    if row["n_perturb"] >= 1:
        for name, (seed, sigma) in list(SEEDS.items())[: row["n_perturb"]]:
            best = None
            for s in range(seed, seed + 200):
                rng = np.random.default_rng(s)
                a = rel.copy(); p = pos0.copy()
                p[movable] += rng.normal(0.0, sigma, size=(movable.sum(), 3))
                a.set_positions(p); md = min_dist_xy_periodic(a)
                if best is None or md > best[0]: best = (md, s, a)
                if md >= 2.6:
                    children[name] = (a, f"gaussian sigma={sigma} A on {movable.sum()} movable atoms, seed {s}, min d {md:.2f} A")
                    break
            else:
                # relaxed parents with already-short contacts (kinks, islands) may not admit 2.6 A at sigma = 0.10 within 200 seeds:
                # take the best seed if it clears 2.5 A and RECORD it; otherwise report loudly instead of dropping the configuration
                if best[0] >= 2.5:
                    children[name] = (best[2], f"gaussian sigma={sigma} A on {movable.sum()} movable atoms, seed {best[1]}, min d {best[0]:.2f} A (fallback: best of 200 seeds, 2.6 A not reachable)")
                else:
                    raise RuntimeError(f"{sid} {name}: no seed of 200 gives min distance >= 2.5 A (best {best[0]:.2f})")
    # collective
    if row["n_collective"] >= 1:
        a = rel.copy(); p = pos0.copy(); top = p[:, 2] > p[:, 2].max() - 0.3
        p[top, 2] -= 0.03 * D111; a.set_positions(p)
        children["coll_topspacing-3pct"] = (a, f"top plane ({top.sum()} atoms) moved down by 0.03*d111 = {0.03*D111:.3f} A")
    if row["n_collective"] >= 2:
        if row["family"] == "strip step":
            a = rel.copy(); p = pos0.copy(); strip = p[:, 2] > p[:, 2].max() - 0.3
            frac = a.get_scaled_positions(wrap=True)
            u = np.round(frac[strip, 0], 4); rows_u = np.sort(np.unique(u))
            edge = strip.copy(); edge[strip] = np.isin(u, [rows_u[0], rows_u[-1]])
            # across-step unit vector = a1 perpendicular component
            a1, a2 = cell[0][:2], cell[1][:2]; t = a2 / np.linalg.norm(a2); perp = a1 - np.dot(a1, t) * t; perp /= np.linalg.norm(perp)
            v = np.round(frac[:, 1] * 2) % 2          # ny = 2 rows -> alternate sign
            for i in np.where(edge)[0]:
                p[i, :2] += (0.15 if v[i] == 0 else -0.15) * perp
            a.set_positions(p)
            children["coll_edgebend0.15"] = (a, f"edge rows ({edge.sum()} atoms) displaced +-0.15 A across the step, alternating by row")
        else:
            a = rel.copy(); c = cell.copy(); c[0] *= 1.01; c[1] *= 1.01
            p = pos0.copy(); p[:, :2] *= 1.01
            a.set_cell(c, scale_atoms=False); a.set_positions(p)
            children["coll_strain+1pct"] = (a, "in-plane cell and xy scaled by 1.01")
    # paths
    if row["n_path"] >= 1 and sid == "A1-fcc":
        # fcc adatom on layer 4 sits above a layer-2 atom (xy); the hcp hollow is above a layer-3 atom (the layer two
        # below the adatom) -> target xy = nearest layer-3 atom, |shift| = a/sqrt6 = 1.70 A; bridge image = half-way
        ad = int(np.argmax(pos0[:, 2]))
        l_n2 = np.where(np.abs(pos0[:, 2] - (pos0[ad, 2] - 2 * D111)) < 0.4)[0]
        shift = min((mic(pos0[k] - pos0[ad], cell) for k in l_n2), key=lambda d: np.linalg.norm(d[:2]))
        assert abs(np.linalg.norm(shift[:2]) - A0 / np.sqrt(6)) < 0.15, np.linalg.norm(shift[:2])
        a = rel.copy(); p = pos0.copy(); p[ad, :2] += 0.5 * shift[:2]; a.set_positions(p)
        children["path_bridge"] = (a, f"adatom moved half-way ({0.5*np.linalg.norm(shift[:2]):.3f} A) toward the nearest hcp hollow (bridge site)")
    # --- Batch B path recipes (2 images each: 0.5 = bridge/midpoint, 1.0 = constructed end point, nothing relaxed) ---
    PATH_MIN_D = 2.2          # path images are transition-like geometries on RELAXED parents (linear interpolation, no z-lift); 2.2 A accepted, value recorded
    def path_clearance(mover, target_xy):
        d = np.zeros((len(rel), 3)); d[mover, :2] = mic(np.array([*target_xy, pos0[mover, 2]]) - pos0[mover], cell)[:2]
        a = rel.copy(); a.set_positions(pos0 + 0.5 * d); return min_dist_xy_periodic(a), d
    def path_images(mover, target_xy, name, note):
        md05, d = path_clearance(mover, target_xy)
        assert md05 >= PATH_MIN_D, f"{name}: mid-image min distance {md05:.2f} A < {PATH_MIN_D}"
        for k, f in enumerate((0.5, 1.0), start=1):
            a = rel.copy(); a.set_positions(pos0 + f * d); md = min_dist_xy_periodic(a)
            children[f"{name}{k}"] = (a, f"{note}; image {f:.1f}, moving atom displaced {np.linalg.norm(f*d[mover]):.2f} A; min d {md:.2f} A (path images accepted >= {PATH_MIN_D} A)")
    def inplane_neighbors(k, zsel, lo=2.4, hi=3.4):
        """in-plane neighbours within [lo, hi] A -- relaxed geometries deviate from the ideal 2.94 A by up to ~0.3 A"""
        return [m for m in zsel if m != k and lo <= np.linalg.norm(mic(pos0[m] - pos0[k], cell)[:2]) <= hi]
    if row["n_path"] >= 1 and sid.startswith("Kink-edge"):
        # kink atom = strip atom with only 2 in-plane strip neighbours; it moves one site along the edge (a2 direction)
        strip = np.where(pos0[:, 2] > pos0[:, 2].max() - 0.3)[0]
        cn = {k: len(inplane_neighbors(k, strip)) for k in strip}
        kink = min(cn, key=cn.get); assert cn[kink] == 2, f"kink atom in-plane neighbours = {cn[kink]} (expected 2); counts {sorted(cn.values())}"
        ny = int(round(np.linalg.norm(cell[1]) / (A0 / np.sqrt(2)))); step_vec = cell[1][:2] / ny
        # both directions along the edge are equivalent end points; take the one with the larger mid-image clearance
        tgt = max((pos0[kink, :2] + sgn * step_vec for sgn in (1, -1)), key=lambda c: path_clearance(kink, c)[0])
        path_images(kink, tgt, "path_kinkmove", "kink atom translated one site along the edge (end point = equivalent kink; direction with the larger mid-image clearance)")
    if row["n_path"] >= 1 and sid.startswith("Island-7"):
        # a ring atom (in-plane CN 3) moves to the empty hollow adjacent to it that is farthest from the island centre
        isl = np.where(pos0[:, 2] > pos0[:, 2].max() - 0.3)[0]
        centre = pos0[isl, :2].mean(axis=0)
        cn = {k: len(inplane_neighbors(k, isl)) for k in isl}
        a1v, a2v = cell[0][:2], cell[1][:2]; ny = int(round(np.linalg.norm(a2v) / (A0 / np.sqrt(2)))); nx = int(round(np.linalg.norm(a1v) / (A0 / np.sqrt(2))))
        e1, e2 = a1v / nx, a2v / ny
        occupied = pos0[isl, :2]; best = None
        for mover in sorted(isl, key=lambda k: cn[k]):            # outer atoms (fewest neighbours) first
            if cn[mover] > 3: continue
            for v in (e1, -e1, e2, -e2, e1 - e2, e2 - e1):
                c = pos0[mover, :2] + v
                if all(np.linalg.norm(mic(np.array([*c, 0]) - np.array([*o, 0]), cell)[:2]) > 1.0 for o in occupied):
                    md05, _ = path_clearance(mover, c)
                    if best is None or md05 > best[0]: best = (md05, mover, c)
        assert best is not None, "no island edge atom with a free adjacent hollow"
        path_images(best[1], best[2], "path_detach", "island edge atom moved to an adjacent empty hollow (pair with the largest mid-image clearance)")
    if row["n_path"] >= 1 and sid.startswith("Pit-7"):
        # a rim atom (top-layer atom next to a vacancy) moves into the adjacent vacancy site
        top = np.where(np.abs(pos0[:, 2] - pos0[:, 2].max()) < 0.3)[0]
        a1v, a2v = cell[0][:2], cell[1][:2]; ny = int(round(np.linalg.norm(a2v) / (A0 / np.sqrt(2)))); nx = int(round(np.linalg.norm(a1v) / (A0 / np.sqrt(2))))
        e1, e2 = a1v / nx, a2v / ny
        best = None
        for k in top:
            for v in (e1, -e1, e2, -e2, e1 - e2, e2 - e1):
                c = pos0[k, :2] + v
                if all(np.linalg.norm(mic(np.array([*c, 0]) - np.array([*pos0[m, :2], 0]), cell)[:2]) > 1.0 for m in top):
                    md05, _ = path_clearance(k, c)
                    if best is None or md05 > best[0]: best = (md05, k, c)
        assert best is not None, "no rim atom with an adjacent vacancy found"
        path_images(best[1], best[2], "path_rimin", "pit rim atom moved into an adjacent vacancy site (pair with the largest mid-image clearance)")
    if row["n_path"] >= 1 and sid == "Step-8x2":
        end = read(f"{P.LIB}/Step-8x2_edge-vacancy_plus_foot-adatom.poscar")
        pe = end.get_positions(); pe[:, 2] += P.ZMIN - pe[:, 2].min()
        assert len(end) == len(rel)
        d = np.array([mic(pe[i] - pos0[i], cell) for i in range(len(rel))])
        assert abs(np.linalg.norm(d, axis=1).max() - A0 / np.sqrt(2)) < 0.3, f"moving atom displacement {np.linalg.norm(d, axis=1).max():.2f} A, expected ~2.94"
        for k, f in enumerate((0.25, 0.5, 0.75), start=1):
            a = rel.copy(); a.set_positions(pos0 + f * d)
            md = min_dist_xy_periodic(a)
            assert md >= 2.4, f"path image {k}: min distance {md:.2f} A"
            children[f"path_detach{k}"] = (a, f"linear image {f:.2f} between relaxed Step-8x2 and the edge-vacancy+foot-adatom end point; "
                                              f"moving atom displaced {np.linalg.norm(f*d, axis=1).max():.2f} A; min d {md:.2f} A (path images accepted >= 2.4 A)")

    if dry:
        print(f"[dry-run] {sid}: N={len(rel)} movable={movable.sum()} potentials for children={mus}")
        for cname, (a, note) in children.items():
            dmax = np.linalg.norm(a.get_positions() - pos0, axis=1).max() if a.get_cell().array[0][0] == cell[0][0] else float("nan")
            print(f"   {cname:24s} min d={min_dist_xy_periodic(a):.3f} A  max atom move={dmax:.3f} A  :: {note}")
        print(f"   would queue: relaxed@{[m for m in row['mu'] if m != P.MU_RELAX]} + {len(children)} configs x {len(mus)} mu = "
              f"{(len(row['mu'])-1) + len(children)*len(mus)} single points")
        return

    # queue (flock shared with the farm jobs)
    with P.queue_lock():
        q = P.load_queue(); made = 0
        # relaxed reference: reuse the relaxation at MU_RELAX, single points at the other potentials
        for mu in row["mu"]:
            task = f"{sid}__relaxed__mu{mu}"
            if any(t["task_id"] == task for t in q): continue
            if mu == P.MU_RELAX:
                q.append(dict(task_id=task, structure_id=sid, family=row["family"], tier=row["tier"], config="relaxed", mu=mu, kind="sp",
                              dir=rdir, n_atoms=len(rel), priority=0, status="reused", job_id=None,
                              note="final ionic step of the accepted relaxation (forces/grids of the final geometry)"))
            else:
                a = rel.copy()
                if P.make_task(q, sid, row, "relaxed", mu, "sp", a, movable, priority=25, note="single point on the relaxed geometry (CONTCAR)"): made += 1
        for cname, (a, note) in children.items():
            for mu in mus:
                if P.make_task(q, sid, row, cname, mu, "sp", a, movable, priority=25, note=note): made += 1
        q.sort(key=lambda t: (t["priority"], t["n_atoms"], t["task_id"]))
        P.save_queue(q)
        P.log(f"[children] {sid}: {len(children)} derived configurations ({', '.join(children)}), {made} new single-point tasks queued")


if __name__ == "__main__":
    # production_children.py <structure_id> [--dry-run [<stand-in CONTCAR>]]
    args = sys.argv[1:]
    dry = "--dry-run" in args
    rest = [a for a in args if a != "--dry-run"]
    main(rest[0], dry=dry, contcar_override=(rest[1] if len(rest) > 1 else None))
