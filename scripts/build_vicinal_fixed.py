#!/usr/bin/env python3
"""Corrected crystallographic vicinal slabs Au(221), Au(332), Au(554) -- 2026-09-27.

Why: the previous builds fed conventional-looking Miller indices to ase.build.surface() with the
fcc PRIMITIVE cell (bulk('Au','fcc') default). Indices are then interpreted in the primitive
reciprocal basis: primitive (221) = cubic (113), primitive (332) = cubic (112), primitive (554)
= cubic (223). Verified on the old files from their atomic positions (surface normal recovered
from nearest-neighbour vectors; plane spacings 1.254 / 0.849 / 0.504 A). Only Au211 (fcc211) was
right. The old files are kept under their TRUE indices with a _retired suffix.

Fix: bulk('Au','fcc', a, cubic=True) + surface(); number of atomic planes chosen for ~9 A metal
thickness (Au211 in the library: 12 planes x 0.849 A = 9.3 A); the step direction (period
a/sqrt2 = 2.94 A) repeated x2 for the MAIN face Au(221) (allows along-edge-independent
perturbations) and x1 for the reference faces Au(332)/Au(554) (ideal only).
Validation printed for every slab: recovered cubic normal, plane spacing vs a/(2*sqrt(h2+k2+l2)),
minimum Au-Au distance (periodic), thickness, in-plane cell.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_vicinal_fixed.py
"""
import itertools
import json
import os
import shutil

import numpy as np
from ase.build import bulk, surface
from ase.io import read, write
from ase.neighborlist import neighbor_list

A0 = 4.158
OUT = "03_pilot/all_defect_structures"
TARGET_THICKNESS = 9.0          # A, metal top to bottom (Au211: 9.34)
# (h,k,l), repeats along the 2.94 A step direction, role
TARGETS = [((2, 2, 1), 2, "main"), ((3, 3, 2), 1, "reference"), ((5, 5, 4), 1, "reference")]


def pad_z(atoms, vacuum=15.0):
    pos = atoms.get_positions(); cell = atoms.get_cell()
    zmin, zmax = pos[:, 2].min(), pos[:, 2].max()
    cell[2] = [0, 0, (zmax - zmin) + 2 * vacuum]
    atoms.set_cell(cell); pos[:, 2] += vacuum - zmin; atoms.set_positions(pos)
    return atoms


def cubic_frame(atoms):
    pos = atoms.get_positions(); cell = atoms.get_cell()
    zmid = (pos[:, 2].min() + pos[:, 2].max()) / 2
    core = [i for i in range(len(atoms)) if abs(pos[i, 2] - zmid) < 2.5][:6]
    vecs = []
    for i in core:
        for j in range(len(atoms)):
            if i == j: continue
            for n1 in (-1, 0, 1):
                for n2 in (-1, 0, 1):
                    d = pos[j] + n1 * cell[0] + n2 * cell[1] - pos[i]
                    if abs(np.linalg.norm(d) - A0 / np.sqrt(2)) < 0.05: vecs.append(d)
    axes = []
    for u, v in itertools.combinations(np.array(vecs), 2):
        if abs(np.dot(u, v)) < 1e-3:
            s = u + v
            if abs(np.linalg.norm(s) - A0) < 0.05:
                s = s / np.linalg.norm(s)
                if not any(abs(abs(np.dot(s, a)) - 1) < 1e-3 for a in axes): axes.append(s)
    for trio in itertools.combinations(range(len(axes)), 3):
        M = np.array(axes)[list(trio)]
        if np.allclose(np.abs(M @ M.T), np.eye(3), atol=1e-3): return M
    return None


def recovered_hkl(atoms):
    M = cubic_frame(atoms)
    n = np.cross(atoms.cell[0], atoms.cell[1]); n /= np.linalg.norm(n)
    h = M @ n
    h = h / np.abs(h[np.abs(h) > 1e-6]).min()
    return tuple(sorted(np.abs(np.round(h, 3)), reverse=True))


def describe(atoms, hkl):
    z = np.sort(np.unique(np.round(atoms.get_positions()[:, 2], 3)))
    d_expect = A0 / (2 * np.sqrt(sum(x * x for x in hkl)))     # all three targets have mixed parity
    i, j, dd = neighbor_list("ijd", atoms, 3.2)
    lens = [np.linalg.norm(atoms.cell[0]), np.linalg.norm(atoms.cell[1])]
    return dict(n_atoms=len(atoms), det_cell_A3=float(np.linalg.det(atoms.get_cell().array)),
                planes=len(z), plane_spacing_A=float(np.diff(z).mean()), plane_spacing_expected_A=float(d_expect),
                thickness_A=float(z[-1] - z[0]), min_AuAu_A=float(dd.min()), in_plane_A=[float(x) for x in lens],
                recovered_cubic_normal=[float(x) for x in recovered_hkl(atoms)],
                atoms_per_plane=float(len(atoms) / len(z)))


def primitive_indices(hkl):
    """Cubic (h,k,l) -> indices in ASE's fcc primitive basis a1=(0,1,1)a/2, a2=(1,0,1)a/2, a3=(1,1,0)a/2:
    h'_i = G.a_i / (2*pi) -> ((k+l)/2, (h+l)/2, (h+k)/2), scaled to integers. Check: (2,2,1)->(3,3,4),
    (3,3,2)->(5,5,6), (5,5,4)->(9,9,10); back-transformed with b_i they give (2,2,1),(3,3,2),(5,5,4)."""
    h, k, l = hkl
    v = np.array([(k + l), (h + l), (h + k)], dtype=float)
    g = np.gcd.reduce(v.astype(int))
    return tuple(int(x) for x in v / g)


# surface() with the PRIMITIVE cell and correctly transformed indices returns the primitive surface
# cell (one atom per plane, contains the a/sqrt2 step vector); with the conventional cell ASE returns a
# 9.3 x 9.3 A supercell for these faces, which is why the conventional route is not used here.
bulk_prim = bulk("Au", "fcc", a=A0)
record = {}
# retire the mislabeled files under their true indices (once)
for old, true in (("Au221", "Au113"), ("Au332", "Au211b"), ("Au554", "Au223")):
    src = f"{OUT}/{old}.poscar"; dst = f"{OUT}/{true}_retired.poscar"
    if os.path.exists(src) and not os.path.exists(dst):
        at = read(src); rec = recovered_hkl(at)
        shutil.move(src, dst)
        record[f"{true}_retired"] = dict(from_file=f"{old}.poscar", recovered_cubic_normal=[float(x) for x in rec], n_atoms=len(at),
                                          note="built with primitive-cell Miller indices; not the face its old name claimed; retired, not in the plan")
        print(f"retired {old}.poscar -> {true}_retired.poscar  (true normal {rec}, {len(at)} atoms)")

for hkl, nrep, role in TARGETS:
    d = A0 / (2 * np.sqrt(sum(x * x for x in hkl)))
    n_planes = int(round(TARGET_THICKNESS / d)) + 1
    # ASE 'layers' = repeats of the surface unit cell along the normal; find the value giving n_planes atomic planes
    best = None
    hkl_p = primitive_indices(hkl)
    for L in range(2, 60):
        s = surface(bulk_prim, hkl_p, layers=L, vacuum=None, periodic=True)
        z = np.unique(np.round(s.get_positions()[:, 2], 3))
        if len(z) >= n_planes:
            best = (L, s); break
    L, slab = best
    # ASE returns SOME basis of the 2-D surface lattice, not the reduced one -> Gauss (Lagrange) reduction
    # of the in-plane vectors so that the shortest one is the a/sqrt2 step vector
    a, b = slab.cell[0].copy(), slab.cell[1].copy()
    for _ in range(50):
        if np.linalg.norm(b) < np.linalg.norm(a): a, b = b, a
        m = int(round(np.dot(a, b) / np.dot(a, a)))
        if m == 0: break
        b = b - m * a
    cell = slab.get_cell(); cell[0], cell[1] = a, b
    slab.set_cell(cell, scale_atoms=False); slab.wrap(eps=1e-8)
    # repeat along the step direction (the in-plane vector of length a/sqrt2)
    lens = [np.linalg.norm(slab.cell[0]), np.linalg.norm(slab.cell[1])]
    istep = int(np.argmin([abs(l - A0 / np.sqrt(2)) for l in lens]))
    assert abs(lens[istep] - A0 / np.sqrt(2)) < 1e-3, f"no 2.94 A step-direction vector in the surface cell: {lens}"
    rep = [1, 1, 1]; rep[istep] = nrep
    slab = slab.repeat(tuple(rep))
    slab = pad_z(slab)
    # right-handed cell (2026-09-27 review): Gauss reduction may swap/flip in-plane vectors and leave det(A) < 0.
    # Keep Cartesian positions and +z, invert one in-plane vector, re-wrap in-plane. Same periodic structure.
    cell = slab.get_cell().array.copy()
    if np.linalg.det(cell) < 0:
        cell[1] *= -1
        slab.set_cell(cell, scale_atoms=False)
        slab.wrap(pbc=(True, True, False), eps=1e-8)
    assert np.linalg.det(slab.get_cell().array) > 0, "left-handed cell after fix"
    info = describe(slab, hkl)
    # direction check: recovered normal (sorted |components|, ratio-normalised) parallel to the target (h,k,l)
    rn = np.array(info["recovered_cubic_normal"]); tn = np.array(sorted(hkl, reverse=True), dtype=float)
    parallel = abs(np.dot(rn, tn) / np.linalg.norm(rn) / np.linalg.norm(tn) - 1.0) < 1e-4
    ok = (parallel and abs(info["plane_spacing_A"] - info["plane_spacing_expected_A"]) < 0.01 and info["min_AuAu_A"] > 2.9)
    name = "Au%d%d%d" % hkl
    info.update(role=role, step_repeats=nrep, ase_layers=L, primitive_indices_used=list(hkl_p), validation="PASS" if ok else "FAIL")
    record[name] = info
    print(f"{name}: {info['n_atoms']} atoms, {info['planes']} planes x {info['plane_spacing_A']:.3f} A (expected {info['plane_spacing_expected_A']:.3f}), "
          f"thickness {info['thickness_A']:.2f} A, in-plane {info['in_plane_A'][0]:.2f} x {info['in_plane_A'][1]:.2f} A, min d {info['min_AuAu_A']:.3f}, "
          f"recovered normal {info['recovered_cubic_normal']}, {info['validation']}")
    if ok:
        write(f"{OUT}/{name}.poscar", slab, format="vasp", direct=False, sort=True)

os.makedirs("03_pilot/report_assets/batch1", exist_ok=True)
json.dump(record, open("03_pilot/report_assets/batch1/vicinal_rebuild_validation.json", "w"), indent=1)
print("wrote 03_pilot/report_assets/batch1/vicinal_rebuild_validation.json")
