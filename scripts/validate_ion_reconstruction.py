#!/usr/bin/env python3
"""Validate the ion-species reconstruction against a CHARGED reference
(T_dUm02, TARGETMU=-5.1071, q_e=-0.164483 -- a real, non-near-zero signal),
not the near-zero Step-8x1_muref point. Reconstruction is done on the FULL
3D grid (Boltzmann relation is nonlinear -- must not average phi/SION first
and apply the formula to the averages). Tests multiple psi zero-point
choices and both sign conventions explicitly, rather than assuming either.
"""
import numpy as np
from ase.io import read

KBT = 8.617333262e-5 * 298.15
N_BULK = 1.0 * 6.02214076e-4
R_ION = 4.0
D_ION = 2 ** (5 / 6) * R_ION
N_MAX = 1.0 / D_ION ** 3
Lz = 44.603


def read_field(path):
    with open(path) as f:
        lines = f.readlines()
    for i, l in enumerate(lines):
        parts = l.split()
        if len(parts) == 3 and all(p.isdigit() for p in parts):
            ngx, ngy, ngz = map(int, parts); dim_idx = i; break
    flat = np.fromstring("".join(lines[dim_idx + 1:]), sep=" ")[:ngx * ngy * ngz]
    return flat.reshape((ngz, ngy, ngx)), ngx, ngy, ngz


def species_3d(psi3d, sion3d):
    u = np.clip(psi3d / KBT, -50, 50)
    Dfac = 1 + (2 * N_BULK / N_MAX) * (np.cosh(u) - 1)
    n_A = sion3d * N_BULK * np.exp(-u) / Dfac  # "exp(-u)" species per v2 doc's n_plus slot
    n_B = sion3d * N_BULK * np.exp(u) / Dfac   # "exp(+u)" species per v2 doc's n_minus slot
    return n_A, n_B


D = "T_dUm02"
V = read(f"{D}/POSCAR").get_volume()
phi, nx, ny, nz = read_field(f"{D}/PHI")
sion, *_ = read_field(f"{D}/SION")
rhoion_raw, *_ = read_field(f"{D}/RHOION")
rhoion = rhoion_raw / V  # e/A^3, full 3D

z = (np.arange(nz) / nz) * Lz
sion_z = sion.mean(axis=(1, 2))
phi_z = phi.mean(axis=(1, 2))
bulk_mask = sion_z > 0.9999
interface_mask = (sion_z > 0.01) & (sion_z < 0.99)
print(f"bulk z: {z[bulk_mask].min():.2f}-{z[bulk_mask].max():.2f}  "
      f"interface z: {z[interface_mask].min():.2f}-{z[interface_mask].max():.2f}")

zero_points = {
    "raw (no shift)": 0.0,
    "bulk-window mean": phi_z[bulk_mask].mean(),
    "deep-vacuum mean (z>40)": phi_z[z > 40].mean(),
}


def rel_l2(a, b, mask=None):
    if mask is not None:
        a = a[mask]; b = b[mask]
    num = np.sqrt(np.sum((a - b) ** 2))
    den = np.sqrt(np.sum(b ** 2))
    return num / den if den > 0 else float("nan")


print(f"\nreference: RHOION full-cell integral = {rhoion_raw.sum()/(nx*ny*nz):.6f} e "
      f"(dn_cp for this point = -0.164483 e)\n")

best = None
for zp_name, phi0 in zero_points.items():
    psi3d = phi - phi0  # broadcast: phi is (ngz,ngy,ngx), phi0 is scalar
    n_A, n_B = species_3d(psi3d, sion)
    for sign_name, recon in [("n_A - n_B", n_A - n_B), ("n_B - n_A", n_B - n_A)]:
        recon_z = recon.mean(axis=(1, 2))
        err_full = rel_l2(recon_z, rhoion_z := rhoion.mean(axis=(1, 2)))
        err_bulk = rel_l2(recon_z, rhoion_z, bulk_mask)
        err_iface = rel_l2(recon_z, rhoion_z, interface_mask)
        recon_integral = recon.sum() / (nx * ny * nz)
        print(f"zero={zp_name:26s} sign={sign_name:10s} "
              f"integral={recon_integral:+.4f}e (ref -0.1645e)  "
              f"relL2: full={err_full:.3f} bulk={err_bulk:.3f} interface={err_iface:.3f}")
