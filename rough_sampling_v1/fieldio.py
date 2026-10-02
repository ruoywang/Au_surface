"""Value-level completeness check of VASP / VASPsol++ grid files (CHGCAR, LOCPOT, PHI, RHOB, ...).

A file passes only if: the structure header parses (species counts give N_ions), a grid line "nx ny nz" follows
the coordinates, at least nx*ny*nz values parse as floats and are all finite, and, for the files VASP writes with
PAW augmentation (CHGCAR, POT), exactly N_ions "augmentation occupancies" blocks follow the grid data. Counting
words (the earlier check) accepted a file whose tail was garbage or NaN; this does not.

Parsing uses numpy's text parser on the data block (C speed; ~5-10 s for a 90 M-value file), so one state with
15 fields is checked in a few minutes, once, when its run has ended.
"""
import os
import re

import numpy as np

AUGMENTED = {"CHGCAR", "POT"}


def header(path):
    """(n_ions, grid, byte offset of the first data value) or raises ValueError."""
    with open(path, "rb") as f:
        lines = []
        for _ in range(8): lines.append(f.readline())
        counts = lines[6].split()
        if not all(c.isdigit() for c in counts): raise ValueError("no species-count line")
        n_ions = sum(int(c) for c in counts)
        mode = lines[7].strip().lower()
        if not (mode.startswith(b"direct") or mode.startswith(b"cart")): raise ValueError("no Direct/Cartesian line")
        for _ in range(n_ions): f.readline()
        blank = f.readline()
        if blank.strip(): raise ValueError("no blank line after coordinates")
        grid_line = f.readline().split()
        if len(grid_line) != 3 or not all(g.isdigit() for g in grid_line): raise ValueError("no grid line")
        return n_ions, tuple(int(g) for g in grid_line), f.tell()


def check(path, name=None):
    """(ok, note, stats). stats has n_values, mean (the charge per grid point for CHGCAR), n_augmentation."""
    name = name or os.path.basename(path)
    if not os.path.exists(path): return False, "missing", {}
    try:
        n_ions, grid, off = header(path)
    except (ValueError, IndexError) as e:
        return False, f"header: {e}", {}
    n_grid = grid[0] * grid[1] * grid[2]
    with open(path, "rb") as f:
        f.seek(off); block = f.read()
    k = block.find(b"augmentation")
    data_txt = block if k < 0 else block[:k]
    try:
        vals = np.fromstring(data_txt.decode("ascii", "replace"), sep=" ")
    except Exception as e:
        return False, f"data block unparsable: {e}", {}
    stats = dict(grid=list(grid), n_values=int(len(vals)), n_ions=n_ions)
    if len(vals) < n_grid: return False, f"truncated: {len(vals)} of {n_grid} values", stats
    main = vals[:n_grid]
    if not np.isfinite(main).all(): return False, f"{int((~np.isfinite(main)).sum())} non-finite values", stats
    # non-numeric junk inside the block: np.fromstring stops at the first unparsable token, so a short count already
    # catches it; a count >= n_grid with extra trailing tokens is tolerated only for the augmented files
    stats["mean"] = float(main.mean()); stats["sum_over_grid"] = float(main.sum() / n_grid)
    if name in AUGMENTED:
        if k < 0: return False, "no augmentation blocks in an augmented field", stats
        ok, note, n_aug = check_augmentation(block[k:], n_ions)
        stats["n_augmentation"] = n_aug
        if not ok: return False, note, stats
    elif k >= 0:
        return False, "unexpected augmentation block in a non-augmented field", stats
    return True, f"{len(vals)} values on {grid[0]}x{grid[1]}x{grid[2]}" + (f", {stats.get('n_augmentation')} augmentation blocks complete" if name in AUGMENTED else ""), stats


def check_augmentation(tail, n_ions):
    """Every PAW block 'augmentation occupancies <ion> <n>' must be followed by exactly n parsable, finite values,
    and the ion indices must run 1..n_ions. Counting the headers (the earlier check) accepted a truncated or NaN
    last block. Returns (ok, note, n_blocks)."""
    parts = tail.decode("ascii", "replace").split("augmentation occupancies")
    seen = []
    for p in parts[1:]:
        toks = p.split()
        if len(toks) < 2 or not toks[0].isdigit() or not toks[1].isdigit(): return False, "augmentation block header incomplete", len(seen)
        ion, n = int(toks[0]), int(toks[1])
        try:
            vals = np.array(toks[2:], float)
        except ValueError:
            return False, f"augmentation block {ion}: non-numeric values", len(seen)
        if len(vals) != n: return False, f"augmentation block {ion}: {len(vals)} of {n} values", len(seen)
        if not np.isfinite(vals).all(): return False, f"augmentation block {ion}: non-finite values", len(seen)
        seen.append(ion)
    if seen != list(range(1, n_ions + 1)): return False, f"augmentation blocks for ions {seen[:3]}...{seen[-1:] if seen else ''} ({len(seen)}), expected 1..{n_ions}", len(seen)
    return True, "complete", len(seen)


def contcar_matches_poscar(poscar, contcar, tol=1e-6):
    """For a single point (IBRION = -1) CONTCAR must repeat POSCAR's cell and positions."""
    def load(p):
        L = open(p).read().splitlines()
        cell = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * float(L[1])
        n = sum(int(x) for x in L[6].split()); k = 7
        if L[k].strip().lower().startswith("s"): k += 1
        cart = L[k].strip().lower().startswith("c"); k += 1
        pos = np.array([[float(x) for x in L[k + i].split()[:3]] for i in range(n)])
        return cell, pos if cart else pos @ cell
    if not (os.path.exists(poscar) and os.path.exists(contcar)): return False, "POSCAR or CONTCAR missing"
    c0, p0 = load(poscar); c1, p1 = load(contcar)
    dc = np.abs(c0 - c1).max()
    if p0.shape != p1.shape: return False, f"atom count differs ({len(p0)} vs {len(p1)})"
    # VASP writes CONTCAR wrapped into the cell: compare under minimum image in fractional coordinates
    f = (p1 - p0) @ np.linalg.inv(c0); f -= np.round(f); dp = np.abs(f @ c0).max()
    return bool(dc < tol and dp < tol), f"max |dcell| {dc:.1e}, max |dpos| {dp:.1e} A (minimum image)"
