"""Minimal LAMMPS data-file I/O for the rough-sampling scripts, shared by md_driver, check_potential and
extract_cells so there is one copy of the conventions.

Conventions
  * ONE atom type. pair_style flare indexes its species by the LAMMPS atom type (pair_flare.cpp: beta_matrices
    [itype-1], cutoff_matrix(itype-1, jtype-1)), so a second "fixed" type would address a species the Au-only
    potential does not have. Frozen atoms are therefore selected by ID ranges (group_lines), never by type.
  * atom ids are the ASE index + 1 and the order is preserved, so every LAMMPS output maps back 1:1.
  * the cell is written in LAMMPS' lower-triangular form; positions are rotated into that frame exactly
    (r' = frac(r) @ cell'), which is the identity for the parents (already lower-triangular).
  * read_data parses "Atoms # atomic" (id type x y z [ix iy iz]) and sorts by id; positions come back as LAMMPS
    left them (wrapped into the periodic box), so displacements must be taken with mic_displacement.
"""
import numpy as np
from ase import Atoms

AU_MASS = 196.967


def lammps_cell(cell):
    """Lower-triangular cell with the same lengths and angles as `cell` (rows a, b, c)."""
    A, B, C = np.asarray(cell, float)
    ax = np.linalg.norm(A)
    bx = B @ A / ax; by = np.sqrt(B @ B - bx * bx)
    cx = C @ A / ax; cy = (B @ C - bx * cx) / by; cz = np.sqrt(C @ C - cx * cx - cy * cy)
    return np.array([[ax, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]])


def to_lammps_frame(at):
    """Positions of `at` expressed in the lower-triangular frame (exact rotation / reflection)."""
    cell = at.get_cell().array; lc = lammps_cell(cell)
    frac = at.get_positions() @ np.linalg.inv(cell)
    return frac @ lc, lc


def write_data(at, path, comment="Au, one atom type; frozen atoms are selected by id in the input script"):
    pos, lc = to_lammps_frame(at)
    n = len(at)
    skew = abs(lc[1, 0]) > 1e-8 or abs(lc[2, 0]) > 1e-8 or abs(lc[2, 1]) > 1e-8
    L = [f"# {comment}", "", f"{n} atoms", "1 atom types", "",
         f"0.0 {lc[0,0]:.10f} xlo xhi", f"0.0 {lc[1,1]:.10f} ylo yhi", f"0.0 {lc[2,2]:.10f} zlo zhi"]
    if skew: L.append(f"{lc[1,0]:.10f} {lc[2,0]:.10f} {lc[2,1]:.10f} xy xz yz")
    L += ["", "Masses", "", f"1 {AU_MASS}", "", "Atoms # atomic", ""]
    L += [f"{i+1} 1 {x:.10f} {y:.10f} {z:.10f}" for i, (x, y, z) in enumerate(pos)]
    open(path, "w").write("\n".join(L) + "\n")
    return lc


def group_lines(name, mask, per_line=40):
    """`group <name> id a:b c d:e ...` lines selecting the atoms where mask is True (LAMMPS ids = index+1).
    Repeated `group` commands with the same name take the union. An all-False mask gives an empty group."""
    ids = np.flatnonzero(np.asarray(mask, bool)) + 1
    if len(ids) == 0: return f"group {name} empty"
    runs = np.split(ids, np.flatnonzero(np.diff(ids) > 1) + 1)
    toks = [f"{r[0]}:{r[-1]}" if len(r) > 1 else f"{r[0]}" for r in runs]
    return "\n".join(f"group {name} id " + " ".join(toks[k:k + per_line]) for k in range(0, len(toks), per_line))


def read_data(path):
    """Atoms (all Au) from a LAMMPS data file written by write_data or by LAMMPS write_data, sorted by id."""
    lines = open(path).read().splitlines()
    n = None; xhi = yhi = zhi = None; xy = xz = yz = 0.0; rows = []
    i = 0
    while i < len(lines):
        s = lines[i].split()
        if len(s) >= 2 and s[1] == "atoms" and s[0].isdigit(): n = int(s[0])
        elif len(s) >= 4 and s[2] == "xlo": xhi = float(s[1]) - float(s[0])
        elif len(s) >= 4 and s[2] == "ylo": yhi = float(s[1]) - float(s[0])
        elif len(s) >= 4 and s[2] == "zlo": zhi = float(s[1]) - float(s[0])
        elif len(s) >= 6 and s[3] == "xy": xy, xz, yz = float(s[0]), float(s[1]), float(s[2])
        elif s and s[0] == "Atoms":
            i += 1
            while i < len(lines) and len(rows) < n:
                t = lines[i].split()
                if t and t[0].isdigit(): rows.append((int(t[0]), float(t[2]), float(t[3]), float(t[4])))
                i += 1
            continue
        i += 1
    rows.sort(); pos = np.array([[r[1], r[2], r[3]] for r in rows])
    cell = np.array([[xhi, 0, 0], [xy, yhi, 0], [xz, yz, zhi]])
    return Atoms("Au" * len(pos), positions=pos, cell=cell, pbc=(True, True, False))


def mic_displacement(cell, r0, r1):
    """r1 - r0 with in-plane minimum image (LAMMPS wraps atoms into the periodic box)."""
    cell = np.asarray(cell); C = cell[:2, :2]; Ci = np.linalg.inv(C)
    d = r1 - r0
    f = d[:, :2] @ Ci; f -= np.round(f); d = d.copy(); d[:, :2] = f @ C
    return d
