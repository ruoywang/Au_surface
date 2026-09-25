#!/usr/bin/env python3
"""Step-8x1 / Step-8x2 calculation-ready inputs for the numerical
consistency test (same physical structure, cell repeated x2 along the
step direction -- does the result change?).

Window convention copied EXACTLY from build_pilot_structures.py (the
already-validated T/V1/A1-fcc batch): Z_PAD_BOTTOM=5.0, LIQUID_EXTENT=20.0,
VACUUM_TOTAL=10.0, D_STERN=2.0, SOL_Z0 inset by 1 layer, bottom 2 layers
fixed, DIPOL=[0.5,0.5,0.5]. Step-8x1 has an added top (strip) layer, same
as A1-fcc's adatom -- so the "z_top_max = z_top_flat + d111" tallest-
feature convention transfers directly.

Step-8x2 is built by tiling the FINISHED Step-8x1 (window already set)
x2 along a2 -- not independently re-centered -- so both share identical
window boundaries and only the along-step repeat count differs.
"""
import json
import os
import numpy as np
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.io import write
from ase import Atom

A0 = 4.158
NLAYERS = 4
BOTTOM_FIXED_LAYERS = 2
Z_PAD_BOTTOM = 5.0
LIQUID_EXTENT = 20.0
VACUUM_TOTAL = 10.0
D_STERN = 2.0
SOL_Z0_INSET_LAYERS = 1
NX = 8
OUTBASE = "03_pilot/structures"


def base_slab(nx, ny, nlayers=NLAYERS):
    slab = fcc111('Au', size=(nx, ny, nlayers), a=A0, vacuum=None, orthogonal=False, periodic=True)
    pos = slab.get_positions()
    pos[:, 2] += (Z_PAD_BOTTOM - pos[:, 2].min())
    slab.set_positions(pos)
    return slab


def build_step_8x1():
    base = base_slab(NX, 1)
    posb = base.get_positions()
    frac = base.get_scaled_positions()
    layer_z = np.sort(np.unique(np.round(posb[:, 2], 3)))
    assert len(layer_z) == NLAYERS
    d111 = layer_z[-1] - layer_z[-2]
    registry_z = layer_z[NLAYERS - 3]
    registry_idx = np.where(np.isclose(posb[:, 2], registry_z, atol=0.05))[0]
    keep = [i for i in registry_idx if frac[i, 0] < 0.5]
    for i in keep:
        x, y, _ = posb[i]
        base.append(Atom('Au', position=[x, y, layer_z[-1] + d111]))

    z_top_flat = layer_z[-1]
    z_top_max = z_top_flat + d111  # tallest feature = the strip, same convention as A1-fcc's adatom
    sol_z0 = layer_z[NLAYERS - 1 - SOL_Z0_INSET_LAYERS]
    sol_z1 = z_top_max + LIQUID_EXTENT
    Lz = Z_PAD_BOTTOM + (z_top_flat - layer_z[0]) + d111 + LIQUID_EXTENT + VACUUM_TOTAL

    cell = base.get_cell()
    cell[2, 2] = Lz
    base.set_cell(cell, scale_atoms=False)
    base.set_pbc(True)

    zs = base.get_positions()[:, 2]
    fixed_mask = np.array([np.any(np.isclose(z, layer_z[:BOTTOM_FIXED_LAYERS], atol=0.05)) for z in zs])
    base.set_constraint(FixAtoms(mask=fixed_mask.tolist()))

    geo = dict(layer_z=layer_z.tolist(), d111=float(d111), z_top_flat=float(z_top_flat),
               z_top_max=float(z_top_max), sol_z0=float(sol_z0), sol_z1=float(sol_z1), Lz=float(Lz),
               n_strip_atoms=len(keep))
    return base, geo


def write_pair(atoms1, geo):
    for tag, mult in [("Step-8x1", 1), ("Step-8x2", 2)]:
        d = f"{OUTBASE}/{tag}"
        os.makedirs(d, exist_ok=True)
        atoms = atoms1 * (1, mult, 1) if mult > 1 else atoms1.copy()
        # rebuild the constraint fresh from z-position rather than trust tiling of the constraint object
        zs = atoms.get_positions()[:, 2]
        layer_z = np.array(geo["layer_z"])
        fixed_mask = np.array([np.any(np.isclose(z, layer_z[:BOTTOM_FIXED_LAYERS], atol=0.05)) for z in zs])
        atoms.set_constraint(FixAtoms(mask=fixed_mask.tolist()))
        write(f"{d}/POSCAR", atoms, format='vasp', direct=False, sort=True)

        ion_z0 = geo["sol_z0"] + D_STERN
        ion_z1 = geo["sol_z1"] - D_STERN
        manifest = dict(
            parent_id=tag, generator="scripts/build_step_dft_pair.py",
            purpose="numerical consistency check: same structure, cell repeated x{} along step direction".format(mult),
            a0_A=A0, supercell_nx_ny_nlayers=[NX, mult, NLAYERS], n_atoms=len(atoms),
            bottom_fixed_layers=BOTTOM_FIXED_LAYERS, boundary_mode="single_sided_CEP_DIP",
            cell_Lz_A=geo["Lz"], SOL_Z0_A=geo["sol_z0"], SOL_Z1_A=geo["sol_z1"],
            ION_Z0_A_derived=float(ion_z0), ION_Z1_A_derived=float(ion_z1),
            D_STERN_A=D_STERN, DIPOL_frac=[0.5, 0.5, 0.5],
            note="window copied from build_pilot_structures.py (T/V1/A1-fcc batch), unmodified. "
                 "Step-8x2 = Step-8x1 tiled x2 along a2 AFTER the window was set, not independently centered.",
        )
        with open(f"{d}/manifest.json", "w") as f:
            json.dump(manifest, f, indent=2)
        print(tag, len(atoms), "atoms ->", d)


if __name__ == "__main__":
    atoms1, geo = build_step_8x1()
    print("window: SOL_Z0=%.3f SOL_Z1=%.3f Lz=%.3f" % (geo["sol_z0"], geo["sol_z1"], geo["Lz"]))
    write_pair(atoms1, geo)
