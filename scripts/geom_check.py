#!/usr/bin/env python3
"""Generic geometry-check utilities per defect_plan.md's field schema.
Not a black box: each function is narrow and its assumptions are stated in
its docstring, so results can be sanity-checked rather than trusted blindly.
"""
import numpy as np


def ensure_z_vacuum(atoms, vacuum=15.0):
    """fcc111(..., vacuum=None) sets the cell's c-vector to exactly the slab's
    z-extent, with no padding -- the slab is then spuriously periodic in z
    (top atom sits image-adjacent to the bottom atom under MIC). Illustration
    renders don't care, but any MIC-based distance check does. Pads z before
    checking; does not touch atom positions or the in-plane cell."""
    a = atoms.copy()
    cell = a.get_cell()
    pos = a.get_positions()
    zmin, zmax = pos[:, 2].min(), pos[:, 2].max()
    cell[2] = [0, 0, (zmax - zmin) + 2 * vacuum]
    a.set_cell(cell)
    pos[:, 2] += vacuum - zmin
    a.set_positions(pos)
    return a


def basic_stats(atoms):
    atoms = ensure_z_vacuum(atoms)
    cell = atoms.get_cell()
    pos = atoms.get_positions()
    d = atoms.get_all_distances(mic=True)
    n = len(atoms)
    d_off = d[np.triu_indices(n, k=1)]
    return {
        "N_total": n,
        "in_plane_cell_vectors": (tuple(cell[0][:2]), tuple(cell[1][:2])),
        "minimum_Au_Au_distance": float(d_off.min()) if len(d_off) else None,
    }


def added_removed_count(atoms, ref_atoms, tol=0.15):
    """Compare top-of-stack occupancy against a matched flat/undefected
    reference of the SAME in-plane cell. Only meaningful when ref_atoms
    shares atoms' (nx,ny) footprint; caller is responsible for that match."""
    return {
        "N_base": len(ref_atoms),
        "N_total": len(atoms),
        "N_net": len(atoms) - len(ref_atoms),
    }


def periodic_image_isolation(defect_xy, cell2d, merge_threshold=3.5):
    """defect_xy: (k,2) array of the defect's own atom (or vacancy-site) xy
    coordinates. cell2d: (a1_xy, a2_xy). Returns the minimum distance from
    any defect point to any OTHER periodic image of the defect set (i.e. the
    defect-to-its-own-image separation), and a finite/merged classification.
    merge_threshold defaults to ~1.2x the Au-Au nearest-neighbor distance."""
    a1, a2 = np.asarray(cell2d[0]), np.asarray(cell2d[1])
    pts = np.asarray(defect_xy)
    best = np.inf
    for i in range(-1, 2):
        for j in range(-1, 2):
            if i == 0 and j == 0:
                continue
            shift = i * a1 + j * a2
            shifted = pts + shift
            dmat = np.linalg.norm(pts[:, None, :] - shifted[None, :, :], axis=-1)
            best = min(best, dmat.min())
    status = "finite" if best > merge_threshold else "merged"
    return {"defect_image_min_distance": float(best), "periodic_connectivity": status}


def analyze_half_coverage_strip(cell2d):
    """For the specific construction in build_step_structure.py: the strip
    occupies fractional a1 in [0, 0.5), full extent along a2. Step edges run
    parallel to a2 at a1-fraction 0 and 0.5. Returns the true perpendicular
    terrace widths (measured along m_hat = a1's component perpendicular to
    a2, NOT half the raw |a1|), plus L_parallel/L_perpendicular."""
    a1, a2 = np.asarray(cell2d[0], dtype=float), np.asarray(cell2d[1], dtype=float)
    t_hat = a2 / np.linalg.norm(a2)
    a1_perp_component = a1 - np.dot(a1, t_hat) * t_hat
    L_perpendicular = np.linalg.norm(a1_perp_component)
    L_parallel = np.linalg.norm(a2)
    return {
        "step_edge_direction": tuple(t_hat),
        "L_parallel": float(L_parallel),
        "L_perpendicular": float(L_perpendicular),
        "terrace_width_upper": float(L_perpendicular / 2),
        "terrace_width_lower": float(L_perpendicular / 2),
    }


if __name__ == "__main__":
    import sys
    from ase.io import read

    path = sys.argv[1]
    atoms = read(path)
    stats = basic_stats(atoms)
    print(path)
    for k, v in stats.items():
        print(f"  {k}: {v}")
    cell = atoms.get_cell()
    strip = analyze_half_coverage_strip((cell[0][:2], cell[1][:2]))
    print("  -- if this is a half-coverage strip along a2 --")
    for k, v in strip.items():
        print(f"  {k}: {v}")
