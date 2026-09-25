#!/usr/bin/env python3
"""Build the Step-8x4 / Step-16x4 / Step-24x4 terrace-width series per
defect_plan.md Correction 2: 4-layer base, half the top face gets one extra
full Au layer, forming a raised strip bounded by two step edges.

Bug fixed here vs the first version: the added layer's in-plane registry
must NOT come from a separately-built 5-layer reference slab. ASE's
fcc111() starts the ABC stacking sequence at a different phase depending on
the requested layer count, so a 5-layer build's top layer does not
reliably sit at the fcc-continuation site relative to a 4-layer build's top
layer -- in the first version it silently landed registry-identical to the
existing top layer (atoms 2.4 A apart, i.e. stacked directly on top of each
other, not in a hollow site). Fix: derive the new layer's registry from
INSIDE the same 4-layer build (same convention already used correctly by
add_adatom() in build_all_defect_types.py: the fcc-continuation site shares
the registry of the layer TWO below the top, i.e. index nlayers-3)."""
import os
import numpy as np
from ase.build import fcc111
from ase.io import write
from ase import Atom

A0 = 4.158
NLAYERS = 4
OUT = "03_pilot/all_defect_structures"
os.makedirs(OUT, exist_ok=True)

SIZES = {"Step-8x4": (8, 4), "Step-16x4": (16, 4), "Step-24x4": (24, 4)}


def base_slab(nx, ny, nlayers=NLAYERS):
    slab = fcc111('Au', size=(nx, ny, nlayers), a=A0, vacuum=None, orthogonal=False, periodic=True)
    pos = slab.get_positions()
    pos[:, 2] += (5.0 - pos[:, 2].min())
    slab.set_positions(pos)
    return slab


for tag, (nx, ny) in SIZES.items():
    base = base_slab(nx, ny)
    posb = base.get_positions()
    frac = base.get_scaled_positions()
    zb = np.sort(np.unique(np.round(posb[:, 2], 2)))
    d111 = zb[-1] - zb[-2]

    # fcc-continuation registry = same xy sites as the layer TWO below the
    # top (index nlayers-3), read from THIS SAME build -- no cross-build ref
    registry_z = zb[NLAYERS - 3]
    registry_idx = np.where(np.isclose(posb[:, 2], registry_z, atol=0.05))[0]

    keep = [i for i in registry_idx if frac[i, 0] < 0.5]
    for i in keep:
        x, y, _ = posb[i]
        base.append(Atom('Au', position=[x, y, zb[-1] + d111]))

    write(f"{OUT}/{tag}.poscar", base, format="vasp", direct=False, sort=True)
    print(tag, len(base), "atoms —", len(keep), "added-layer atoms; expected N_base=",
          4 * nx * ny, "N_strip=", int(0.5 * nx * ny))
