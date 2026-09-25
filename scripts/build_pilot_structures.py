#!/usr/bin/env python3
"""Build T / V1 / A1-fcc (v2 doc Sec.3.2) with the locally-fitted a0=4.158 A
(02_bulk_eos/RESULT.md), sharing ONE common single-sided window across all
three (v2 Sec.2.2: same defects batch must share z-layout / window, sized to
cover the tallest feature -- here A1-fcc's added atom)."""
import json
import os
import numpy as np
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.io import write

A0 = 4.158  # locally fitted, see 02_bulk_eos/RESULT.md
NX, NY, NLAYERS = 4, 4, 4  # project decision 2026-09-22: 4 layers for speed (see 00_audit/parameter_map.md H);
                            # overrides v2 doc Sec.3.2's "4 layers = low-cost pretest only" caution
BOTTOM_FIXED_LAYERS = 2
Z_PAD_BOTTOM = 5.0
LIQUID_EXTENT = 20.0
VACUUM_TOTAL = 10.0
D_STERN = 2.0
SOL_Z0_INSET_LAYERS = 1
SEED = 20260921

OUTBASE = "03_pilot/structures"


def base_slab():
    slab = fcc111('Au', size=(NX, NY, NLAYERS), a=A0, vacuum=None, orthogonal=False, periodic=True)
    pos = slab.get_positions()
    pos[:, 2] += (Z_PAD_BOTTOM - pos[:, 2].min())
    slab.set_positions(pos)
    return slab


def common_geometry(slab):
    zs = slab.get_positions()[:, 2]
    layer_z = np.sort(np.unique(np.round(zs, 3)))
    assert len(layer_z) == NLAYERS
    z_top_flat = layer_z[-1]
    # A1-fcc raises the top by one (111) layer spacing -> tallest feature in this batch
    d111 = layer_z[-1] - layer_z[-2]
    z_top_max = z_top_flat + d111
    sol_z0 = layer_z[NLAYERS - 1 - SOL_Z0_INSET_LAYERS]
    sol_z1 = z_top_max + LIQUID_EXTENT
    Lz = Z_PAD_BOTTOM + (z_top_flat - layer_z[0]) + d111 + LIQUID_EXTENT + VACUUM_TOTAL
    return dict(layer_z=layer_z, z_top_flat=z_top_flat, d111=d111, z_top_max=z_top_max,
                sol_z0=sol_z0, sol_z1=sol_z1, Lz=Lz)


def finalize_and_write(tag, atoms, geo, extra_manifest=None):
    d = f"{OUTBASE}/{tag}"
    os.makedirs(d, exist_ok=True)
    cell = atoms.get_cell()
    cell[2, 2] = geo["Lz"]
    atoms.set_cell(cell, scale_atoms=False)
    atoms.set_pbc(True)

    zs = atoms.get_positions()[:, 2]
    layer_z = geo["layer_z"]
    fixed_mask = np.array([np.any(np.isclose(z, layer_z[:BOTTOM_FIXED_LAYERS], atol=0.05)) for z in zs])
    atoms.set_constraint(FixAtoms(mask=fixed_mask.tolist()))

    write(f"{d}/POSCAR", atoms, format='vasp', direct=False, sort=True)

    ion_z0 = geo["sol_z0"] + D_STERN
    ion_z1 = geo["sol_z1"] - D_STERN
    manifest = dict(
        parent_id=tag, generator="scripts/build_pilot_structures.py", seed=SEED,
        a0_A=A0, a0_status="LOCALLY_FITTED_7PT_EOS_see_02_bulk_eos_RESULT.md",
        supercell=[NX, NY, NLAYERS], n_atoms=len(atoms),
        bottom_fixed_layers=BOTTOM_FIXED_LAYERS,
        boundary_mode="single_sided_CEP_DIP",
        cell_Lz_A=float(geo["Lz"]),
        SOL_Z0_A=float(geo["sol_z0"]), SOL_Z1_A=float(geo["sol_z1"]),
        ION_Z0_A_derived=float(ion_z0), ION_Z1_A_derived=float(ion_z1),
        D_STERN_A=D_STERN, DIPOL_frac=[0.5, 0.5, 0.5],
        note="SOL_Z0/SOL_Z1 shared across T/V1/A1-fcc batch (v2 Sec.2.2), sized to A1-fcc's raised atom.",
    )
    if extra_manifest:
        manifest.update(extra_manifest)
    with open(f"{d}/manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(tag, len(atoms), "atoms ->", d)
    return manifest


def build_T(geo):
    slab = base_slab()
    return finalize_and_write("T", slab, geo, dict(defect_type="T (flat, unreconstructed, no defect)", defect_centers=[]))


def build_V1(geo):
    slab = base_slab()
    zs = slab.get_positions()[:, 2]
    top_idx = np.where(np.isclose(zs, geo["layer_z"][-1], atol=0.05))[0]
    pos = slab.get_positions()
    xy = pos[top_idx, :2]
    center = xy.mean(axis=0)
    d2 = ((xy - center) ** 2).sum(axis=1)
    remove_i = top_idx[np.argmin(d2)]
    removed_pos = pos[remove_i].tolist()
    del slab[remove_i]
    return finalize_and_write("V1", slab, geo,
                              dict(defect_type="V1 (single top-layer vacancy)",
                                   defect_centers=[removed_pos]))


def build_A1_fcc(geo):
    slab = base_slab()
    pos = slab.get_positions()
    frac = slab.get_scaled_positions()
    top_idx = np.where(np.isclose(pos[:, 2], geo["layer_z"][-1], atol=0.05))[0]
    # fcc continuation site: same in-plane stacking registry as the layer TWO below the top
    below2_idx = np.where(np.isclose(pos[:, 2], geo["layer_z"][-3], atol=0.05))[0]
    ref_xy_frac = frac[below2_idx[0], :2]
    new_frac = np.array([ref_xy_frac[0], ref_xy_frac[1], 0.0])
    new_cart = slab.get_cell().T @ [new_frac[0], new_frac[1], 0.0]
    new_z = geo["z_top_flat"] + geo["d111"]
    from ase import Atom
    slab.append(Atom('Au', position=[new_cart[0], new_cart[1], new_z]))
    return finalize_and_write("A1_fcc", slab, geo,
                              dict(defect_type="A1-fcc (single adatom, fcc continuation site)",
                                   defect_centers=[[float(new_cart[0]), float(new_cart[1]), float(new_z)]]))


if __name__ == "__main__":
    geo = common_geometry(base_slab())
    print("shared geometry:", {k: v for k, v in geo.items() if k != "layer_z"})
    build_T(geo)
    build_V1(geo)
    build_A1_fcc(geo)
