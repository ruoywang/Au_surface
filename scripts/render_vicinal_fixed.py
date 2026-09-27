#!/usr/bin/env python3
"""Re-render the gallery figures of the corrected Au(221)/Au(332)/Au(554) slabs (2026-09-27), replacing
the figures of the retired mislabeled builds. Reuses render_batch1.do_vicinal unchanged; its side-view
panel tiles the slab along x, so each slab is first rotated in-plane so that its terrace-width vector
(the long in-plane vector) points along +x and the 2.94 A step vector has a positive y component.
Rotated copies are written next to the figures (render-only; library and production POSCARs untouched).
Usage (from Au_Cl/):  scripts/pyrun.sh scripts/render_vicinal_fixed.py
"""
import os
import sys

import numpy as np
from ase.io import read, write

sys.path.insert(0, "scripts")
import render_batch1 as R  # noqa: E402  (module-level: loads MANIFEST, defines helpers; main guarded)

OUT = "03_pilot/report_assets/batch1"
CASES = [  # tag, family text, microfacet, n_periods for the cross-section tiling
    ("Au221", "Au(221) = 4(111)x(111)", "n(111)x(111), with Au(332)/Au(554)", "{111} -> B-type", 3),
    ("Au332", "Au(332) = 6(111)x(111)", "n(111)x(111), with Au(221)/Au(554)", "{111} -> B-type", 2),
    ("Au554", "Au(554) = 10(111)x(111)", "n(111)x(111), with Au(221)/Au(332)", "{111} -> B-type", 2),
]
for tag, facet, family, micro, nper in CASES:
    a = read(f"{R.SRC}/{tag}.poscar")
    cell = a.get_cell().array.copy()
    iw = int(np.argmax([np.linalg.norm(cell[0]), np.linalg.norm(cell[1])]))   # width vector = long one
    w = cell[iw][:2]; ang = -np.arctan2(w[1], w[0])
    Rm = np.array([[np.cos(ang), -np.sin(ang), 0], [np.sin(ang), np.cos(ang), 0], [0, 0, 1]])
    pos = a.get_positions() @ Rm.T; cell2 = cell @ Rm.T
    if cell2[1 - iw][1] < 0:                      # mirror y so the step vector has +y (keeps det sign via z? no: flip y and z-order)
        pos[:, 1] *= -1; cell2[:, 1] *= -1
        cell2[2] = [0, 0, cell2[2][2]]
    if np.linalg.det(cell2) < 0:
        cell2[1 - iw] *= -1
    b = a.copy(); b.set_cell(cell2, scale_atoms=False); b.set_positions(pos); b.wrap(pbc=(True, True, False), eps=1e-8)
    tmp = f"_render_{tag}.poscar"
    write(f"{R.SRC}/{tmp}", b, format="vasp", direct=False, sort=True)
    try:
        R.do_vicinal(tag, tmp, facet, family, micro, iw, nper)
        print(f"{tag}: rendered from {tmp} (width axis index {iw}, |width| = {np.linalg.norm(cell2[iw]):.2f} A, {nper} periods) -> {OUT}/{tag}.png, {tag}_whitebg.png")
    finally:
        os.remove(f"{R.SRC}/{tmp}")
