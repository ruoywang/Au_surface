#!/usr/bin/env python3
"""Mechanism check for Step-8x1_muref vs Step-8x1_dUp02, rev 2 -- fixes 3
statistical inconsistencies from rev 1 (identified before over-interpreting
the rev-1 numbers):

1. edge_dist used the raw (sheared) a1 instead of the true perpendicular
   L_perp -- so "2/3/4 A" bands were actually ~1.73/2.60/3.46 A. Fixed:
   reuses the same L_perp-based distance as the main analysis script.
2. edge1/edge2 "far" was `~near` for THAT edge alone, which silently
   included the other edge's near-region. Fixed: far = excludes BOTH edges.
3. Accessibility-relative sampling averaged a per-column ratio (a different,
   differently-weighted statistic from the main script's whole-region
   Sigma(n)/Sigma(S)). Fixed: primary metric is the whole-region ratio,
   with the per-column-average kept as an explicitly separate, labeled
   diagnostic -- not blended into the same number.

Also implements the K_ij cross-term decomposition: K_ij mixes the
occupation function g[u_i] (from potential field i) with the accessibility
weight S_j (from state j's SION field) -- a post-processing diagnostic to
separate "response driven by potential change" from "response driven by
accessibility/geometry change." K10, K01 are NOT self-consistent physical
states -- they are cross-combinations for this decomposition only.
"""
import numpy as np
from ase.io import read

KBT = 8.6173857e-5 * 298.0  # confirmed local default (solvation.F:1703), not overridden by INCAR
N_BULK = 1.0 * 6.022e-4     # local MOLAR const (solvation.F:3161), not 6.02214076e-4
R_ION = 4.0
N_MAX = 1.0 / (2 ** (5 / 6) * R_ION) ** 3
Lz = 44.603
Z_METAL_TOP = 14.603
NEAR_SURFACE_DEPTH = 15.0
EDGE_BAND = 3.0


def read_field(path):
    with open(path) as f:
        lines = f.readlines()
    for i, l in enumerate(lines):
        parts = l.split()
        if len(parts) == 3 and all(p.isdigit() for p in parts):
            ngx, ngy, ngz = map(int, parts); dim_idx = i; break
    flat = np.fromstring("".join(lines[dim_idx + 1:]), sep=" ")[:ngx * ngy * ngz]
    return flat.reshape((ngz, ngy, ngx)), ngx, ngy, ngz


def g_occupation(u):
    """The species occupation function itself, exp(-u)/D(u) -- WITHOUT the
    S_ion or n_bulk factors, so it can be combined with a DIFFERENT state's
    S_ion for the K_ij cross decomposition."""
    u = np.clip(u, -50, 50)
    D = 1 + (2 * N_BULK / N_MAX) * (np.cosh(u) - 1)
    return np.exp(-u) / D


def load(D):
    atoms = read(f"{D}/POSCAR")
    a1 = atoms.get_cell()[0][0]
    a2v = np.array(atoms.get_cell()[1][:2])
    a1v = np.array([atoms.get_cell()[0][0], atoms.get_cell()[0][1]])
    t_hat = a2v / np.linalg.norm(a2v)
    a1_perp = a1v - np.dot(a1v, t_hat) * t_hat
    L_perp = np.linalg.norm(a1_perp)
    phi, ngx, ngy, ngz = read_field(f"{D}/PHI")
    sion, *_ = read_field(f"{D}/SION")
    z = (np.arange(ngz) / ngz) * Lz
    x = (np.arange(ngx) / ngx) * a1
    return dict(phi=phi, sion=sion, z=z, x=x, a1=a1, L_perp=L_perp, ngx=ngx, ngy=ngy, ngz=ngz)


A = load("Step-8x1_muref")
B = load("Step-8x1_dUp02")
assert A["ngx"] == B["ngx"] and A["ngz"] == B["ngz"], "grids must match for a same-position comparison"
assert abs(A["L_perp"] - B["L_perp"]) < 1e-6, "terrace width must match between the two states"
L_perp = A["L_perp"]

sion_z_A = A["sion"].mean(axis=(1, 2))
near_surface = (sion_z_A > 0.01) & (A["z"] > Z_METAL_TOP) & (A["z"] < Z_METAL_TOP + NEAR_SURFACE_DEPTH)


def edge_distance(x, a1, edges_frac, L_perp):
    u = x / a1
    d = np.full_like(u, np.inf)
    for ue in edges_frac:
        du = np.abs(((u - ue + 0.5) % 1.0) - 0.5)
        d = np.minimum(d, du)
    return d * L_perp


def region_KD(n_anion_3d, sion_3d, zmask, xmask):
    n = n_anion_3d[zmask][:, :, xmask].sum()
    s = sion_3d[zmask][:, :, xmask].sum()
    return n / (N_BULK * s) if s > 0 else float("nan")


def n_anion_field(state):
    return state["sion"] * N_BULK * g_occupation(state["phi"] / KBT)


print(f"true terrace width each side = {L_perp/2:.2f} A (perpendicular measurement)\n")

print(f"=== 1. Bandwidth sensitivity (near-edge band = 2/3/4 A, TRUE perpendicular distance, both edges lumped) ===")
d_both = np.minimum(edge_distance(A["x"], A["a1"], [0.0], L_perp),
                     edge_distance(A["x"], A["a1"], [0.5], L_perp))
nA = n_anion_field(A); nB = n_anion_field(B)
for band in [2.0, 3.0, 4.0]:
    near = d_both < band; far = ~near
    KA_n, KA_f = region_KD(nA, A["sion"], near_surface, near), region_KD(nA, A["sion"], near_surface, far)
    KB_n, KB_f = region_KD(nB, B["sion"], near_surface, near), region_KD(nB, B["sion"], near_surface, far)
    dKn, dKf = KB_n - KA_n, KB_f - KA_f
    print(f"band={band:.0f}A  near_frac={near.mean():.2f}  "
          f"K_A(near/far)={KA_n:.4f}/{KA_f:.4f}  K_B(near/far)={KB_n:.4f}/{KB_f:.4f}  "
          f"dK_near={dKn:.4f} dK_far={dKf:.4f}  ratio={dKn/dKf:.3f}")

print(f"\n=== 2. Edge1 vs edge2 separately (band={EDGE_BAND:.0f}A, 'far' EXCLUDES both edges) ===")
far_common = d_both >= EDGE_BAND
for ue, name in [(0.0, "edge1 (u=0)"), (0.5, "edge2 (u=0.5)")]:
    near = edge_distance(A["x"], A["a1"], [ue], L_perp) < EDGE_BAND
    KA_n, KA_f = region_KD(nA, A["sion"], near_surface, near), region_KD(nA, A["sion"], near_surface, far_common)
    KB_n, KB_f = region_KD(nB, B["sion"], near_surface, near), region_KD(nB, B["sion"], near_surface, far_common)
    dKn, dKf = KB_n - KA_n, KB_f - KA_f
    print(f"{name}: K_A(near/far)={KA_n:.4f}/{KA_f:.4f}  K_B(near/far)={KB_n:.4f}/{KB_f:.4f}  "
          f"dK_near={dKn:.4f} dK_far(common terrace, both edges excluded)={dKf:.4f}  ratio={dKn/dKf:.3f}")

print(f"\n=== 3. Accessibility-relative sampling (per-column local SION=0.5 boundary, "
      f"10A window ABOVE that local boundary) ===")
sion_col_A = A["sion"].mean(axis=1)
z_bound = np.full(A["ngx"], np.nan)
for i in range(A["ngx"]):
    col = sion_col_A[:, i]
    above = np.where(col > 0.5)[0]
    if len(above):
        z_bound[i] = A["z"][above[0]]
print(f"local boundary z range across columns: {np.nanmin(z_bound):.2f}-{np.nanmax(z_bound):.2f} A "
      f"(real corrugation from the step)")

WIN = 10.0
near3 = d_both < EDGE_BAND
far3 = d_both >= EDGE_BAND
# PRIMARY metric: same whole-region Sigma(n)/Sigma(S) convention as sections 1-2,
# just with a per-column z-window instead of one fixed absolute window
for label, state, nfield in [("A", A, nA), ("B", B, nB)]:
    n_near = s_near = n_far = s_far = 0.0
    for i in range(state["ngx"]):
        zlo, zhi = z_bound[i], z_bound[i] + WIN
        zmask_col = (state["z"] > zlo) & (state["z"] < zhi)
        n_i = nfield[zmask_col, :, i].sum()
        s_i = state["sion"][zmask_col, :, i].sum()
        if near3[i]:
            n_near += n_i; s_near += s_i
        else:
            n_far += n_i; s_far += s_i
    K_n = n_near / (N_BULK * s_near); K_f = n_far / (N_BULK * s_far)
    print(f"  [{label}, whole-region ratio, per-column window] K_near={K_n:.4f} K_far={K_f:.4f}")
    if label == "A":
        KA_n_acc, KA_f_acc = K_n, K_f
    else:
        KB_n_acc, KB_f_acc = K_n, K_f
dKn_acc, dKf_acc = KB_n_acc - KA_n_acc, KB_f_acc - KA_f_acc
print(f"  PRIMARY (whole-region ratio): dK_near={dKn_acc:.4f} dK_far={dKf_acc:.4f} ratio={dKn_acc/dKf_acc:.3f}")

# SEPARATE diagnostic: per-column equal-weighted average (a DIFFERENT statistic,
# reported but not blended with the primary number above)
col_KD = {}
for label, state, nfield in [("A", A, nA), ("B", B, nB)]:
    vals = []
    for i in range(state["ngx"]):
        zlo, zhi = z_bound[i], z_bound[i] + WIN
        zmask_col = (state["z"] > zlo) & (state["z"] < zhi)
        n_i = nfield[zmask_col, :, i].sum(); s_i = state["sion"][zmask_col, :, i].sum()
        vals.append(n_i / (N_BULK * s_i) if s_i > 0 else np.nan)
    col_KD[label] = np.array(vals)
KA_n_col, KA_f_col = np.nanmean(col_KD["A"][near3]), np.nanmean(col_KD["A"][far3])
KB_n_col, KB_f_col = np.nanmean(col_KD["B"][near3]), np.nanmean(col_KD["B"][far3])
dKn_col, dKf_col = KB_n_col - KA_n_col, KB_f_col - KA_f_col
print(f"  SEPARATE diagnostic (per-column equal-weighted average, NOT the primary metric): "
      f"dK_near={dKn_col:.4f} dK_far={dKf_col:.4f} ratio={dKn_col/dKf_col:.3f}")

print(f"\n=== 4. K_ij cross-term decomposition (near-surface region, band={EDGE_BAND:.0f}A) ===")
print("K_ij: occupation function g[u_i] from potential field i, combined with SION from state j.")
print("K10, K01 are post-processing cross-combinations for this decomposition ONLY -- not new self-consistent states.\n")
gA = g_occupation(A["phi"] / KBT)
gB = g_occupation(B["phi"] / KBT)
for region_name, xmask in [("near-edge", near3), ("far-from-edge (common terrace)", far3)]:
    def K(g, s):
        num = (s[near_surface][:, :, xmask] * g[near_surface][:, :, xmask]).sum()
        den = s[near_surface][:, :, xmask].sum()
        return num / den if den > 0 else float("nan")
    K00 = K(gA, A["sion"]); K11 = K(gB, B["sion"])
    K10 = K(gB, A["sion"]); K01 = K(gA, B["sion"])
    dK_phi = 0.5 * ((K10 - K00) + (K11 - K01))
    dK_S = 0.5 * ((K01 - K00) + (K11 - K10))
    print(f"[{region_name}] K00={K00:.4f} K10={K10:.4f} K01={K01:.4f} K11={K11:.4f}  "
          f"dK_phi={dK_phi:+.4f} dK_S={dK_S:+.4f}  (sum={dK_phi+dK_S:+.4f}, K11-K00={K11-K00:+.4f})")

print(f"\n=== 5. Potential difference map, same positions (psi_B - psi_A) ===")
dphi_xz = (B["phi"] - A["phi"]).mean(axis=1)
dphi_near = dphi_xz[near_surface][:, near3].mean()
dphi_far = dphi_xz[near_surface][:, far3].mean()
print(f"unweighted mean Delta-psi: near-edge={dphi_near:+.4f} V  far-from-edge={dphi_far:+.4f} V  "
      f"diff={dphi_near-dphi_far:+.4f} V")
print("NOTE: this is an UNWEIGHTED spatial average, not the S_ion-weighted average the K_D "
      "metric actually uses -- similar unweighted Delta-psi does NOT by itself rule out an "
      "electrostatic contribution to the differential response (see K_ij decomposition above).")
