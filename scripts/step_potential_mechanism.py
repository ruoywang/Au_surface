#!/usr/bin/env python3
"""Mechanism check for Step-8x1_muref vs Step-8x1_dUp02: is the weaker
near-edge response a real, robust local effect, or an artifact of (a) fixed
near-edge bandwidth, (b) lumping edge1+edge2 together, or (c) sampling by
absolute height rather than distance from the (corrugated) accessible
boundary? No new DFT -- reuses the already-computed, already-validated
fields (T=298.0 fix applied)."""
import numpy as np
from ase.io import read

KBT = 8.6173857e-5 * 298.0
N_BULK = 1.0 * 6.02214076e-4
R_ION = 4.0
N_MAX = 1.0 / (2 ** (5 / 6) * R_ION) ** 3
Lz = 44.603
Z_METAL_TOP = 14.603
NEAR_SURFACE_DEPTH = 15.0


def read_field(path):
    with open(path) as f:
        lines = f.readlines()
    for i, l in enumerate(lines):
        parts = l.split()
        if len(parts) == 3 and all(p.isdigit() for p in parts):
            ngx, ngy, ngz = map(int, parts); dim_idx = i; break
    flat = np.fromstring("".join(lines[dim_idx + 1:]), sep=" ")[:ngx * ngy * ngz]
    return flat.reshape((ngz, ngy, ngx)), ngx, ngy, ngz


def load(D):
    atoms = read(f"{D}/POSCAR")
    V = atoms.get_volume()
    a1 = atoms.get_cell()[0][0]
    phi, ngx, ngy, ngz = read_field(f"{D}/PHI")
    sion, *_ = read_field(f"{D}/SION")
    u = np.clip(phi / KBT, -50, 50)
    Dfac = 1 + (2 * N_BULK / N_MAX) * (np.cosh(u) - 1)
    n_anion = sion * N_BULK * np.exp(-u) / Dfac
    z = (np.arange(ngz) / ngz) * Lz
    x = (np.arange(ngx) / ngx) * a1
    return dict(phi=phi, sion=sion, n_anion=n_anion, z=z, x=x, a1=a1, ngx=ngx, ngy=ngy, ngz=ngz)


A = load("Step-8x1_muref")
B = load("Step-8x1_dUp02")
assert A["ngx"] == B["ngx"] and A["ngz"] == B["ngz"], "grids must match for a same-position comparison"

sion_z_A = A["sion"].mean(axis=(1, 2))
near_surface = (sion_z_A > 0.01) & (A["z"] > Z_METAL_TOP) & (A["z"] < Z_METAL_TOP + NEAR_SURFACE_DEPTH)


def edge_dist(x, a1, ue):
    u = x / a1
    return np.abs(((u - ue + 0.5) % 1.0) - 0.5) * a1  # a1 already scaled correctly for THIS cell


def KD(state, xmask, zmask=near_surface):
    n = state["n_anion"][zmask][:, :, xmask].sum()
    s = state["sion"][zmask][:, :, xmask].sum()
    return n / (N_BULK * s) if s > 0 else float("nan")


print("=== 1. Bandwidth sensitivity (near-edge band = 2, 3, 4 A; both edges lumped) ===")
for band in [2.0, 3.0, 4.0]:
    d = edge_dist(A["x"], A["a1"], 0.0)
    d = np.minimum(d, edge_dist(A["x"], A["a1"], 0.5))
    near = d < band; far = ~near
    KA_n, KA_f = KD(A, near), KD(A, far)
    KB_n, KB_f = KD(B, near), KD(B, far)
    dKn, dKf = KB_n - KA_n, KB_f - KA_f
    print(f"band={band:.0f}A  near_frac={near.mean():.2f}  "
          f"K_A(near/far)={KA_n:.4f}/{KA_f:.4f}  K_B(near/far)={KB_n:.4f}/{KB_f:.4f}  "
          f"dK_near={dKn:.4f} dK_far={dKf:.4f}  ratio={dKn/dKf:.3f}")

print("\n=== 2. Edge1 vs edge2 separately (band=3A, not lumped) ===")
for ue, name in [(0.0, "edge1 (u=0)"), (0.5, "edge2 (u=0.5)")]:
    d = edge_dist(A["x"], A["a1"], ue)
    near = d < 3.0; far = ~near
    KA_n, KA_f = KD(A, near), KD(A, far)
    KB_n, KB_f = KD(B, near), KD(B, far)
    dKn, dKf = KB_n - KA_n, KB_f - KA_f
    print(f"{name}: K_A(near/far)={KA_n:.4f}/{KA_f:.4f}  K_B(near/far)={KB_n:.4f}/{KB_f:.4f}  "
          f"dK_near={dKn:.4f} dK_far={dKf:.4f}  ratio={dKn/dKf:.3f}")

print("\n=== 3. Accessibility-relative sampling (per-column local SION=0.5 boundary, "
      "fixed 10A window ABOVE that local boundary, not a fixed absolute z) ===")
sion_col_A = A["sion"].mean(axis=1)  # (ngz, ngx), averaged over the short along-step axis
z_bound = np.full(A["ngx"], np.nan)
for i in range(A["ngx"]):
    col = sion_col_A[:, i]
    above = np.where(col > 0.5)[0]
    if len(above):
        z_bound[i] = A["z"][above[0]]  # first z where accessibility exceeds 0.5, per column
print(f"local boundary z range across columns: {np.nanmin(z_bound):.2f} - {np.nanmax(z_bound):.2f} A "
      f"(corrugation from the step)")

WIN = 10.0
n_anion_acc = {}
for label, state in [("A", A), ("B", B)]:
    vals = []
    for i in range(state["ngx"]):
        zlo, zhi = z_bound[i], z_bound[i] + WIN
        zmask_col = (state["z"] > zlo) & (state["z"] < zhi)
        n_col = state["n_anion"][zmask_col, :, i].sum()
        s_col = state["sion"][zmask_col, :, i].sum()
        vals.append(n_col / (N_BULK * s_col) if s_col > 0 else np.nan)
    n_anion_acc[label] = np.array(vals)

d3 = edge_dist(A["x"], A["a1"], 0.0)
d3 = np.minimum(d3, edge_dist(A["x"], A["a1"], 0.5))
near3, far3 = d3 < 3.0, d3 >= 3.0
KA_n = np.nanmean(n_anion_acc["A"][near3]); KA_f = np.nanmean(n_anion_acc["A"][far3])
KB_n = np.nanmean(n_anion_acc["B"][near3]); KB_f = np.nanmean(n_anion_acc["B"][far3])
dKn, dKf = KB_n - KA_n, KB_f - KA_f
print(f"K_A(near/far)={KA_n:.4f}/{KA_f:.4f}  K_B(near/far)={KB_n:.4f}/{KB_f:.4f}  "
      f"dK_near={dKn:.4f} dK_far={dKf:.4f}  ratio={dKn/dKf:.3f}")

print("\n=== 4. Potential difference map, same positions (psi_B - psi_A), no re-derivation of sign ===")
dphi_xz = (B["phi"] - A["phi"]).mean(axis=1)  # (ngz, ngx)
dphi_col_near = dphi_xz[near_surface][:, near3].mean()
dphi_col_far = dphi_xz[near_surface][:, far3].mean()
print(f"mean (psi_dUp02 - psi_muref) in near-surface region: near-edge columns={dphi_col_near:+.4f} V  "
      f"far-from-edge columns={dphi_col_far:+.4f} V  diff={dphi_col_near-dphi_col_far:+.4f} V")
