#!/usr/bin/env python3
"""Systematic v2 Sec.5-style QC on all 9 pilot points: bulk PHI plateau flatness,
SDIEL/SION plateau, RHOB+RHOION charge closure vs. the CP-loop's own Delta q_e."""
import numpy as np

Lz = 44.603
SOL_Z0, SOL_Z1 = 9.801, 34.603
NEUTRAL = {"T": 704, "V1": 693, "A1_fcc": 715}
POINTS = ["T_dUm02", "T_dUp00", "T_dUp02", "V1_dUm02", "V1_dUp00", "V1_dUp02",
          "A1_fcc_dUm02", "A1_fcc_dUp00", "A1_fcc_dUp02"]
FINAL_NELE = {
    "T_dUm02": 703.835517, "T_dUp00": 704.000000, "T_dUp02": 704.166912,
    "V1_dUm02": 692.843135, "V1_dUp00": 693.010903, "V1_dUp02": 693.180361,
    "A1_fcc_dUm02": 714.761478, "A1_fcc_dUp00": 714.938718, "A1_fcc_dUp02": 715.117065,
}


def read_field(path):
    with open(path) as f:
        lines = f.readlines()
    dim_idx = None
    for i, l in enumerate(lines):
        parts = l.split()
        if len(parts) == 3 and all(p.isdigit() for p in parts):
            dim_idx = i
            ngx, ngy, ngz = map(int, parts)
            break
    data_lines = lines[dim_idx + 1:]
    flat = np.fromstring("".join(data_lines), sep=" ")
    n = ngx * ngy * ngz
    flat = flat[:n]
    grid = flat.reshape((ngz, ngy, ngx))
    return grid, ngx, ngy, ngz


for tag in POINTS:
    d = tag
    tag_struct = tag.split("_dU")[0]
    try:
        phi, ngx, ngy, ngz = read_field(f"{d}/PHI")
        sdiel, *_ = read_field(f"{d}/SDIEL")
        sion, *_ = read_field(f"{d}/SION")
        rhob, *_ = read_field(f"{d}/RHOB")
        rhoion, *_ = read_field(f"{d}/RHOION")
    except FileNotFoundError as e:
        print(f"{tag}: MISSING FIELD FILE {e}")
        continue

    z = (np.arange(ngz) / ngz) * Lz
    phi_z = phi.mean(axis=(1, 2))
    sdiel_z = sdiel.mean(axis=(1, 2))
    sion_z = sion.mean(axis=(1, 2))

    mask = (z > SOL_Z0 + 8) & (z < SOL_Z1 - 6)
    mask_ion = (z > SOL_Z0 + 12) & (z < SOL_Z1 - 10)

    phi_bulk_mean = phi_z[mask].mean()
    phi_bulk_std = phi_z[mask].std()
    sdiel_plateau = sdiel_z[mask].mean()
    sion_plateau = sion_z[mask_ion].mean()
    sion_std = sion_z[mask_ion].std()

    cell_a = 11.7606
    area = 119.781711
    dV = (area * Lz) / (ngx * ngy * ngz)
    closure = (rhob.sum() + rhoion.sum()) * dV

    dn_cp = FINAL_NELE[tag] - NEUTRAL[tag_struct]

    # expected: solvent/ion response charge should cancel the metal's excess/deficit,
    # i.e. integral(RHOB+RHOION) approx -dn_cp for overall periodic-cell neutrality
    residual = closure - (-dn_cp)
    print(f"{tag:16s} phi_bulk={phi_bulk_mean:+.2e}+-{phi_bulk_std:.1e}  "
          f"SDIEL={sdiel_plateau:.6f}  SION={sion_plateau:.6f}(+-{sion_std:.1e})  "
          f"int(RHOB+RHOION)={closure:+.4f}e  -dn_cp(expected)={-dn_cp:+.4f}e  residual={residual:+.2e}e")
