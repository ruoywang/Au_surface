#!/usr/bin/env python3
"""Ion species reconstruction, rev 2 -- fixes 4 confirmed bugs from rev 1:

1. psi zero point: NOT phi minus a finite-window bulk average. Validated
   (see validate_ion_reconstruction.py) that raw PHI, with no shift at all,
   reconstructs RHOION to <2.1e-7 e/A^3 max error on the FULL 3D grid of a
   real charged reference (T_dUm02, q_e=-0.164483, not a near-zero point).
   The solver's own PHI already carries the correct physical zero point.
2. Species identity: empirically anchored, not assumed from the v2 doc's
   variable names. n_A=exp(-u) matches RHOION as n_A-n_B (to 1e-7), AND is
   the species enriched near this positively-charged electrode -- which by
   simple electrostatic screening MUST be the anion. So: n_A=anion,
   n_B=cation. (The doc's own "n_plus"/"n_minus" labels do not mean
   cation/anion -- do not assume that naming.)
3. Reconstruction done on the full 3D grid first (nonlinear relation: you
   cannot average phi/SION first and then apply the Boltzmann formula).
   Z-profiles and cross-step maps are averages of the 3D result, not the
   other way around.
4. Enrichment ratio K_D = integral(n_anion)/(, n_bulk*integral(S_ion)) over
   a region -- numerator does NOT get multiplied by S_ion again (n_anion
   already includes it). Edge distance uses proper minimum-image wrapping
   and the true perpendicular terrace width, not the raw cell vector.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read

KBT = 8.6173857e-5 * 298.0  # BOLKEV (stepver.F) and code default SolTemp (solvation.F:1703,
                             # not overridden by INCAR) -- was 298.15, confirmed wrong: at
                             # 298.15 all 4 validation points showed relL2~5e-4; at 298.0
                             # they drop to ~2.3e-5 (20x), i.e. floating-point precision
N_BULK = 1.0 * 6.022e-4  # local MOLAR const (solvation.F:3161), not 6.02214076e-4
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


def species_3d(phi3d, sion3d):
    """Full-grid reconstruction. Returns (n_anion, n_cation), both e/A^3."""
    u = np.clip(phi3d / KBT, -50, 50)
    Dfac = 1 + (2 * N_BULK / N_MAX) * (np.cosh(u) - 1)
    n_anion = sion3d * N_BULK * np.exp(-u) / Dfac
    n_cation = sion3d * N_BULK * np.exp(u) / Dfac
    return n_anion, n_cation


def edge_distance(x, a1, edges_frac, L_perp):
    """Minimum-image distance from each x (Cartesian, along a1) to the
    nearest of the given fractional-position edges, in units of the TRUE
    perpendicular terrace width L_perp (u=0..1 spans L_perp in real
    perpendicular distance, not the raw sheared lattice vector |a1|)."""
    u = x / a1
    d = np.full_like(u, np.inf)
    for ue in edges_frac:
        du = np.abs(((u - ue + 0.5) % 1.0) - 0.5)
        d = np.minimum(d, du)
    return d * L_perp


def enrichment_KD(n_anion_3d, sion_3d, xmask):
    """K_D = integral(n_anion)/( n_bulk * integral(S_ion) ) over the
    x-columns selected by xmask, no extra S_ion factor in the numerator."""
    num = n_anion_3d[:, :, xmask].sum()
    den = N_BULK * sion_3d[:, :, xmask].sum()
    return num / den if den > 0 else float("nan")


if __name__ == "__main__":
    import sys
    D = sys.argv[1] if len(sys.argv) > 1 else "Step-8x1_muref"
    atoms = read(f"{D}/POSCAR")
    V = atoms.get_volume()
    cell = atoms.get_cell()
    a1 = cell[0][0]

    phi, ngx, ngy, ngz = read_field(f"{D}/PHI")
    sion, *_ = read_field(f"{D}/SION")
    rhoion_raw, *_ = read_field(f"{D}/RHOION")
    rhoion = rhoion_raw / V

    z = (np.arange(ngz) / ngz) * Lz
    x = (np.arange(ngx) / ngx) * a1

    # --- full 3D reconstruction first ---
    n_anion, n_cation = species_3d(phi, sion)

    # sanity unit test: if psi were forced to 0 everywhere, n_anion should be
    # exactly n_bulk*SION (K_D=1 in a pure-S_ion region) -- check the formula,
    # not the real field, as a regression guard
    test_S = np.array([0.3, 0.7, 1.0])
    test_n, _ = species_3d(np.zeros_like(test_S), test_S)
    assert np.allclose(test_n, test_S * N_BULK), "unit test failed: psi=0 should give n=n_bulk*S_ion"

    # --- validation against RHOION (full 3D, not a lenient near-zero window) ---
    recon = n_anion - n_cation
    diff = recon - rhoion
    print(f"{D}: validation vs RHOION (full 3D grid): max|diff|={np.abs(diff).max():.3e} e/A^3, "
          f"RMS={np.sqrt((diff**2).mean()):.3e} e/A^3, RHOION max|.|={np.abs(rhoion).max():.3e} e/A^3")

    sion_z = sion.mean(axis=(1, 2))
    bulk_mask = sion_z > 0.9999
    n_anion_z = n_anion.mean(axis=(1, 2))
    n_cation_z = n_cation.mean(axis=(1, 2))

    # near-surface analysis region: tied to the ACTUAL metal top (from atom
    # positions, not a SION threshold), extending a fixed depth into the
    # liquid. This explicitly excludes the far artificial solvent-window
    # transition (near SOL_Z1), and does NOT stop at SION<0.98 -- the
    # response can extend into the fully-accessible (SION=1) region too
    # (confirmed: T_dUm02 shows n_anion/n_bulk=1.08 even where SION=1).
    z_metal_top = atoms.get_positions()[:, 2].max()
    NEAR_SURFACE_DEPTH = 15.0
    interface_mask = (sion_z > 0.01) & (z > z_metal_top) & (z < z_metal_top + NEAR_SURFACE_DEPTH)
    print(f"near-surface region: z={z_metal_top:.2f} (metal top) to "
          f"{z_metal_top+NEAR_SURFACE_DEPTH:.2f} A ({interface_mask.sum()} grid points along z)")

    def rel_l2(mask):
        a, b = recon[mask], rhoion[mask]
        den = np.sqrt(np.sum(b ** 2))
        return np.sqrt(np.sum((a - b) ** 2)) / den if den > 0 else float("nan")

    z3 = np.broadcast_to(z[:, None, None], recon.shape)
    m_bulk3 = np.broadcast_to(bulk_mask[:, None, None], recon.shape)
    m_iface3 = np.broadcast_to(interface_mask[:, None, None], recon.shape)
    m_all3 = np.ones_like(recon, dtype=bool)
    print(f"  relative L2 error: full-cell={rel_l2(m_all3):.4f}  "
          f"bulk(SION=1, far-field diagnostic)={rel_l2(m_bulk3):.4f}  "
          f"near-surface(metal-top to +{NEAR_SURFACE_DEPTH:.0f}A)={rel_l2(m_iface3):.4f}")

    # --- output 1: n_anion/n_cation z-profiles ---
    fig, ax = plt.subplots(figsize=(7, 5), dpi=200)
    ax.plot(z, n_anion_z / N_BULK, label="n_anion/n_bulk", color="#a83d2f")
    ax.plot(z, n_cation_z / N_BULK, label="n_cation/n_bulk", color="#2f6fa8")
    ax.plot(z, sion_z, label="SION (accessibility)", color="#888", ls="--", lw=1)
    ax.axvspan(z[bulk_mask].min(), z[bulk_mask].max(), color="#2f6fa8", alpha=0.06, label="SION=1 window")
    ax.set_xlabel("z (A)"); ax.set_ylabel("n / n_bulk, or SION"); ax.legend(fontsize=9)
    ax.set_title(f"{D}: ion species profiles (3D-reconstructed, then averaged)")
    plt.tight_layout()
    plt.savefig(f"../03_pilot/report_assets/batch1/{D}_ion_zprofile.png", dpi=200)
    plt.close()

    # --- output 2: cross-step (x,z) map ---
    n_anion_xz = n_anion.mean(axis=1)
    sion_xz = sion.mean(axis=1)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=200, sharey=True)
    for ax, field, title, cmap in [(axes[0], sion_xz, "SION", "viridis"),
                                    (axes[1], n_anion_xz / N_BULK, "n_anion/n_bulk", "magma")]:
        im = ax.pcolormesh(x, z, field, shading="auto", cmap=cmap)
        plt.colorbar(im, ax=ax, fraction=0.046)
        ax.set_title(title, fontsize=11); ax.set_xlabel("x, across-step (A)")
    axes[0].set_ylabel("z (A)")
    plt.suptitle(f"{D}: cross-step section", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"../03_pilot/report_assets/batch1/{D}_crossstep_map.png", dpi=200)
    plt.close()

    print(f"n_anion/n_bulk in SION=1 bulk window: mean={n_anion_z[bulk_mask].mean()/N_BULK:.4f}")
    print(f"n_cation/n_bulk in SION=1 bulk window: mean={n_cation_z[bulk_mask].mean()/N_BULK:.4f}")

    # --- enrichment K_D: near-edge vs far-from-edge, split into the deep-bulk
    # reservoir window (diagnostic only) and the near-interface EDL region
    # (the one that actually carries a step-induced signal, if any) ---
    if len(cell) >= 2:
        a2 = cell[1][:2]
        a1v = np.array([cell[0][0], cell[0][1]])
        t_hat = a2 / np.linalg.norm(a2)
        a1_perp = a1v - np.dot(a1v, t_hat) * t_hat
        L_perp = np.linalg.norm(a1_perp)
        d_edge = edge_distance(x, a1, [0.0, 0.5], L_perp)
        near = d_edge < 3.0
        far = ~near
        print(f"\nterrace width each side = {L_perp/2:.2f} A (true perpendicular measurement); "
              f"near-edge band = 3 A; near fraction={near.mean():.2f} far fraction={far.mean():.2f}")
        if far.mean() < 0.15:
            print("  NOTE: far-from-edge region is a thin sliver -- not a step-unaffected terrace interior; "
                  "report as a narrow-periodic-array conditional value only.")
        for label, zmask, zname in [(bulk_mask, bulk_mask, "deep-bulk reservoir window, far-field (diagnostic)"),
                                     (interface_mask, interface_mask, "near-surface region, metal-top to +15A (physical signal)")]:
            KD_near = enrichment_KD(n_anion[zmask][:, :, :], sion[zmask][:, :, :], near)
            KD_far = enrichment_KD(n_anion[zmask][:, :, :], sion[zmask][:, :, :], far)
            print(f"  [{zname}] K_D near-edge={KD_near:.4f}  K_D far-from-edge={KD_far:.4f}  "
                  f"ratio={KD_near/KD_far:.4f}")
