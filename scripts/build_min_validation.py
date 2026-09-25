#!/usr/bin/env python3
"""Build the flat-T, single-sided minimal validation case (no defect).

Purpose: empirically check the interface conclusions in 00_audit/model_conventions.md
(NESCHEME=3 static convergence, ISOL=2 PHI semantics, single-sided window profile,
charge closure) on a small, cheap case, before generating the full 8-defect batch.
This is NOT a physically calibrated production point: a0 and TARGETMU below are
explicitly flagged placeholders, not validated values.
"""
import json
import numpy as np
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.io import write

# ---- inputs (see README.md in the output directory for justification) ----
A0_PBE_PLACEHOLDER = 4.180   # Angstrom. Literature-typical PBE/PAW Au fcc lattice constant.
                              # NOT locally validated: v2 doc Sec.3.1 requires a local 5-7 volume
                              # EOS fit before production use. Placeholder only for this plumbing test.
NX, NY, NLAYERS = 4, 4, 6
BOTTOM_FIXED_LAYERS = 2
Z_PAD_BOTTOM = 5.0           # Angstrom, vacuum gap below the fixed bottom Au layer (user request:
                              # shift the metal up ~5-6 A so DIPOL=0.5,0.5,0.5 sits in a sane spot
                              # rather than right at the periodic boundary)
# LIQUID_EXTENT/VACUUM_TOTAL: v2 doc's own "40/20 A" (line 94) has NO stated derivation (confirmed
# by grep -- not from an explicit-water reference, not from a Debye-length estimate, just an
# unlabeled starting guess). Re-derived here from the 1:1 electrolyte Debye length at 1 M,
# lambda_D[nm] = 0.304/sqrt(C[M]) => lambda_D(1M) = 3.0 A. A few lambda_D (~15 A) should reach the
# bulk plateau; + ~5 A genuine plateau margin (v2 Sec.5 acceptance criterion) + far-window buffer.
# STILL a starting point, not a validated value -- must be checked against the Sec.5 plateau tests.
LIQUID_EXTENT = 20.0         # Angstrom, from topmost Au atom to the far solvent-vacuum window (SOL_Z1)
VACUUM_TOTAL = 10.0          # Angstrom, connected vacuum thickness above the liquid window
D_STERN = 2.0                # Angstrom, must match the code's D_STERN default (solvation.F:1722)
SOL_Z0_INSET_LAYERS = 1      # place SOL_Z0 at the layer just below the top (metal-interior anchor)
SEED = 20260921

OUTDIR = "01_min_validation/_shared/T_flat"


def build():
    slab = fcc111('Au', size=(NX, NY, NLAYERS), a=A0_PBE_PLACEHOLDER,
                   vacuum=None, orthogonal=False, periodic=True)
    pos = slab.get_positions()
    pos[:, 2] += (Z_PAD_BOTTOM - pos[:, 2].min())
    slab.set_positions(pos)

    zs = slab.get_positions()[:, 2]
    layer_z = np.sort(np.unique(np.round(zs, 3)))
    assert len(layer_z) == NLAYERS, f"expected {NLAYERS} distinct layers, got {len(layer_z)}"

    layer_index = np.array([int(np.argmin(np.abs(layer_z - z))) for z in zs])  # 0 = bottom
    fixed_mask = layer_index < BOTTOM_FIXED_LAYERS
    slab.set_constraint(FixAtoms(mask=fixed_mask.tolist()))

    z_top_au = zs.max()
    slab_thickness = zs.max() - zs.min()
    Lz = Z_PAD_BOTTOM + slab_thickness + LIQUID_EXTENT + VACUUM_TOTAL

    cell = slab.get_cell()
    cell[2, 2] = Lz
    slab.set_cell(cell, scale_atoms=False)
    slab.set_pbc(True)

    # single-sided window (solvation.F:1503-1541): window == 1 for SOL_Z0 < z < SOL_Z1
    sol_z0 = layer_z[NLAYERS - 1 - SOL_Z0_INSET_LAYERS]   # one layer below the top Au layer
    sol_z1 = z_top_au + LIQUID_EXTENT
    ion_z0 = sol_z0 + D_STERN
    ion_z1 = sol_z1 - D_STERN

    # DIPOL: user requested the plain default (0.5, 0.5, 0.5) rather than a computed centroid --
    # per audit finding, DIPOL is only the multipole reference center (not the jump-plane location),
    # so this simple, conventional choice is fine as long as the metal isn't sitting right at the
    # periodic boundary (hence Z_PAD_BOTTOM=5.0 above).
    dipol_frac = 0.5

    import os
    os.makedirs(OUTDIR, exist_ok=True)
    write(f"{OUTDIR}/POSCAR", slab, format='vasp', direct=False, sort=True)

    manifest = dict(
        parent_id="T_flat_single_sided_minvalidation",
        generator="scripts/build_min_validation.py",
        seed=SEED,
        a0_A=A0_PBE_PLACEHOLDER,
        a0_status="PLACEHOLDER_NOT_LOCALLY_VALIDATED",
        supercell=[NX, NY, NLAYERS],
        n_atoms=len(slab),
        bottom_fixed_layers=BOTTOM_FIXED_LAYERS,
        defect_centers=[],
        defect_type="T (flat, unreconstructed, no defect)",
        boundary_mode="single_sided_CEP_DIP",
        z_top_au_A=float(z_top_au),
        slab_thickness_A=float(slab_thickness),
        liquid_extent_A=LIQUID_EXTENT,
        vacuum_total_A=VACUUM_TOTAL,
        cell_Lz_A=float(Lz),
        SOL_Z0_A=float(sol_z0),
        SOL_Z1_A=float(sol_z1),
        ION_Z0_A_derived=float(ion_z0),
        ION_Z1_A_derived=float(ion_z1),
        D_STERN_A=D_STERN,
        DIPOL_frac=[0.5, 0.5, float(dipol_frac)],
    )
    with open(f"{OUTDIR}/manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest


if __name__ == "__main__":
    m = build()
    print(json.dumps(m, indent=2))
