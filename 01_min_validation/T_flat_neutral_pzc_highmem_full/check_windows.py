import numpy as np


def read_field(fname):
    with open(fname) as f:
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
    return grid.mean(axis=(1, 2)), ngz


Lz = 47.066621
SOL_Z0, SOL_Z1 = 14.653, 37.067
z_top_au = 17.066621

for fname in ["SION", "SDIEL", "PHI"]:
    prof, ngz = read_field(fname)
    z = (np.arange(ngz) / ngz) * Lz
    plateau_mask = (z > SOL_Z0 + 8) & (z < SOL_Z1 - 4)
    print(f"{fname}: plateau mean={prof[plateau_mask].mean():.6f} std={prof[plateau_mask].std():.2e}, "
          f"min={prof.min():.4f} max={prof.max():.4f}")
    if fname == "PHI":
        # locate largest single-step jump (|delta| between consecutive z) -> dipole jump plane
        dphi = np.abs(np.diff(prof))
        jmax = np.argmax(dphi)
        print(f"  largest jump in PHI: {dphi[jmax]:.4f} between z={z[jmax]:.3f} and z={z[jmax+1]:.3f}")
        # print full low-res profile for visual sanity check (every ~1A)
        step = max(1, ngz // 60)
        for i in range(0, ngz, step):
            print(f"    z={z[i]:7.3f}  PHI={prof[i]:10.5f}")
