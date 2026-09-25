import numpy as np

with open("PHI") as f:
    lines = f.readlines()

# find the grid-dimension line (three ints)
dim_idx = None
for i, l in enumerate(lines):
    parts = l.split()
    if len(parts) == 3 and all(p.isdigit() for p in parts):
        dim_idx = i
        ngx, ngy, ngz = map(int, parts)
        break

print(f"grid dims: {ngx} {ngy} {ngz}, found at line {dim_idx}")

data_lines = lines[dim_idx + 1:]
flat = np.fromstring("".join(data_lines), sep=" ")
n_expected = ngx * ngy * ngz
print(f"read {flat.size} values, expected {n_expected}")
flat = flat[:n_expected]

# VASP CHGCAR-style order: x fastest, then y, then z
grid = flat.reshape((ngz, ngy, ngx))
phi_z = grid.mean(axis=(1, 2))

# cell c-vector length (z extent), from POSCAR header read earlier
Lz = 47.066621
z_coords = (np.arange(ngz) / ngz) * Lz

SOL_Z0, SOL_Z1 = 14.653, 37.067
np.savetxt("phi_z_profile.dat", np.column_stack([z_coords, phi_z]),
           header="z(A)  phi_planeavg(PHI file units)")

# bulk plateau window: well inside the window, away from EDL and far transition
plateau_mask = (z_coords > SOL_Z0 + 8) & (z_coords < SOL_Z1 - 4)
phi_bulk = phi_z[plateau_mask].mean()
phi_bulk_std = phi_z[plateau_mask].std()

# vacuum region (outside the window) for comparison
vac_mask = (z_coords < 5) | (z_coords > SOL_Z1 + 3)
phi_vac = phi_z[vac_mask].mean()
phi_vac_std = phi_z[vac_mask].std()

EFERMI = -4.9411  # eV, from OUTCAR

print(f"plateau window z in ({SOL_Z0+8:.2f}, {SOL_Z1-4:.2f}) A, n_pts={plateau_mask.sum()}")
print(f"phi_bulk (liquid plateau) mean={phi_bulk:.6f}, std={phi_bulk_std:.6f}")
print(f"phi_vac (vacuum region)   mean={phi_vac:.6f}, std={phi_vac_std:.6f}")
print(f"EFERMI (raw, OUTCAR) = {EFERMI} eV")
print(f"mu0 = EFERMI - phi_bulk = {EFERMI - phi_bulk:.6f} eV")
print(f"(for reference) EFERMI - phi_vac = {EFERMI - phi_vac:.6f} eV")
