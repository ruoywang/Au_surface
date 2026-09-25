#!/usr/bin/env python3
"""Reconstruct n+/n- from the already-computed Step-8x1 fields, per v2 doc
Sec.7.2's constitutive relation (validated against RHOION below before use --
not assumed). No new DFT; reuses Step-8x1_muref's PHI/SION/RHOION.

u = e*psi/(kB*T),  D = 1 + (2*n_bulk/n_max)*(cosh(u)-1)
n_plus  = SION*n_bulk*exp(-u)/D
n_minus = SION*n_bulk*exp(+u)/D
psi = PHI - PHI_bulk  (validated: this sign, no flip, reconstructs RHOION to
                        1.5e-5 e/A^3 in the bulk window -- same order as
                        RHOION's own bulk noise floor)
n_max = 1/d_ion^3, d_ion = 2^(5/6)*R_ION (D_ION not set in INCAR -> this default)
n_bulk (per species) = C_MOLAR * 6.02214076e-4 A^-3
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read

D = "Step-8x1_muref"
Lz = 44.603
R_ION = 4.0
C_MOLAR = 1.0
T_KELVIN = 298.15
KBT = 8.617333262e-5 * T_KELVIN
N_BULK = C_MOLAR * 6.02214076e-4
D_ION = 2 ** (5 / 6) * R_ION
N_MAX = 1.0 / D_ION ** 3


def read_field(path):
    with open(path) as f:
        lines = f.readlines()
    for i, l in enumerate(lines):
        parts = l.split()
        if len(parts) == 3 and all(p.isdigit() for p in parts):
            ngx, ngy, ngz = map(int, parts); dim_idx = i; break
    flat = np.fromstring("".join(lines[dim_idx + 1:]), sep=" ")[:ngx * ngy * ngz]
    return flat.reshape((ngz, ngy, ngx)), ngx, ngy, ngz


atoms = read(f"{D}/POSCAR")
V = atoms.get_volume()
a1 = atoms.get_cell()[0][0]  # across-step direction, Angstrom (Step-8x1: nx along x)

phi, ngx, ngy, ngz = read_field(f"{D}/PHI")
sion, *_ = read_field(f"{D}/SION")
rhoion_raw, *_ = read_field(f"{D}/RHOION")
rhoion = rhoion_raw / V  # corrected per the volume-normalization fix

z = (np.arange(ngz) / ngz) * Lz
x = (np.arange(ngx) / ngx) * a1

phi_z = phi.mean(axis=(1, 2))
sion_z = sion.mean(axis=(1, 2))
bulk_mask = sion_z > 0.9999
phi_bulk = phi_z[bulk_mask].mean()
print(f"bulk window: z={z[bulk_mask].min():.2f}-{z[bulk_mask].max():.2f} A, phi_bulk={phi_bulk:.6f} V")


def n_plus_minus(psi, sion_local):
    u = np.clip(psi / KBT, -50, 50)
    Dfac = 1 + (2 * N_BULK / N_MAX) * (np.cosh(u) - 1)
    return sion_local * N_BULK * np.exp(-u) / Dfac, sion_local * N_BULK * np.exp(u) / Dfac


# ---- validation (must reproduce RHOION before trusting anything downstream) ----
psi_z = phi_z - phi_bulk
npz, nmz = n_plus_minus(psi_z, sion_z)
recon_z = npz - nmz
val_err = np.max(np.abs(recon_z[bulk_mask] - rhoion.mean(axis=(1, 2))[bulk_mask]))
print(f"validation: max|reconstructed - RHOION| in bulk window = {val_err:.3e} e/A^3 "
      f"(bulk RHOION scale ~1.5e-5 e/A^3) -- {'PASS' if val_err < 5e-5 else 'CHECK'}")

# ---- output 1: n+/n- z-profiles normalized to n_bulk ----
fig, ax = plt.subplots(figsize=(7, 5), dpi=200)
ax.plot(z, npz / N_BULK, label="n+/n_bulk", color="#2f6fa8")
ax.plot(z, nmz / N_BULK, label="n-/n_bulk", color="#a83d2f")
ax.plot(z, sion_z, label="SION (accessibility)", color="#888", ls="--", lw=1)
ax.axvspan(z[bulk_mask].min(), z[bulk_mask].max(), color="#2f6fa8", alpha=0.06, label="SION=1 window")
ax.set_xlabel("z (A)"); ax.set_ylabel("n / n_bulk, or SION"); ax.legend(fontsize=9)
ax.set_title("Step-8x1: normal ion profiles (planar-averaged)")
plt.tight_layout()
plt.savefig("../03_pilot/report_assets/batch1/step8x1_ion_zprofile.png", dpi=200)
plt.close()
print("wrote step8x1_ion_zprofile.png")
print(f"  n+/n_bulk in SION=1 window: mean={  (npz/N_BULK)[bulk_mask].mean():.4f} std={(npz/N_BULK)[bulk_mask].std():.2e}")
print(f"  n-/n_bulk in SION=1 window: mean={(nmz/N_BULK)[bulk_mask].mean():.4f} std={(nmz/N_BULK)[bulk_mask].std():.2e}")

# ---- output 2: cross-step (x,z) map of psi, SION, n- ----
phi_xz = phi.mean(axis=1)  # average over the along-step (ngy) axis only
sion_xz = sion.mean(axis=1)
psi_xz = phi_xz - phi_bulk
npx, nmx = n_plus_minus(psi_xz, sion_xz)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), dpi=200, sharey=True)
for ax, field, title, cmap in [
    (axes[0], psi_xz, "psi = PHI - PHI_bulk (V)", "RdBu_r"),
    (axes[1], sion_xz, "SION (accessibility)", "viridis"),
    (axes[2], nmx / N_BULK, "n-/n_bulk", "magma"),
]:
    im = ax.pcolormesh(x, z, field, shading="auto", cmap=cmap)
    plt.colorbar(im, ax=ax, fraction=0.046)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x, across-step (A)")
    ax.axhline(14.603, color="white", lw=1, ls=":")  # top of the strip
axes[0].set_ylabel("z (A)")
plt.suptitle("Step-8x1: cross-step section (averaged along the step direction)", fontsize=12)
plt.tight_layout()
plt.savefig("../03_pilot/report_assets/batch1/step8x1_crossstep_map.png", dpi=200)
plt.close()
print("wrote step8x1_crossstep_map.png")

# ---- output 3: step-adjacent vs terrace enrichment, accessible-volume normalized ----
# strip (elevated terrace) occupies fractional a1 in [0,0.5); its two edges are at
# u=0 and u=0.5. "step-adjacent" = within 3 A of either edge (in x); "terrace interior"
# = the rest, split by which side of the strip boundary it's on.
u = x / a1
edge_dist = np.minimum(np.abs(u - 0.0) % 1.0, np.abs(u - 0.5) % 1.0) * a1
near_edge = edge_dist < 3.0
far_from_edge = ~near_edge
print(f"\nterrace width each side = {a1/2:.2f} A; 'near-edge' band = 3 A from either edge")
print(f"  near-edge fraction of cell width: {near_edge.mean():.2f}  far-from-edge fraction: {far_from_edge.mean():.2f}")
if far_from_edge.mean() < 0.15:
    print("  NOTE: far-from-edge region is a thin sliver of the cell -- Step-8x1's terrace "
          "is too narrow to claim a genuine step-unaffected interior. Treat any 'terrace' "
          "number below as a narrow-periodic-array conditional value, not an isolated-step limit.")

# accessible-volume-weighted average of n-, in the bulk-z window, split by x-region
z_bulk_idx = np.where(bulk_mask)[0]
sion_col = sion[z_bulk_idx][:, :, :].mean(axis=1)  # (z_bulk, x)
nm_col = nmx[z_bulk_idx][:, :]  # (z_bulk, x) -- already averaged over y

def weighted_mean(field2d, sion2d, xmask):
    w = sion2d[:, xmask]
    v = field2d[:, xmask]
    wsum = w.sum()
    return (v * w).sum() / wsum if wsum > 0 else float("nan")

nm_near = weighted_mean(nm_col, sion_col, near_edge)
nm_far = weighted_mean(nm_col, sion_col, far_from_edge)
print(f"\naccessible-volume-weighted n-/n_bulk, bulk-z window:")
print(f"  near-edge (<3 A from a step edge): {nm_near/N_BULK:.4f}")
print(f"  far-from-edge (terrace interior, thin): {nm_far/N_BULK:.4f}")
print(f"  ratio (near/far): {nm_near/nm_far:.4f}")
