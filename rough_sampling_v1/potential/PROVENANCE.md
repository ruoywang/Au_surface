# Candidate-generation potential: provenance and terms of use

This potential generates **candidate** rough Au surfaces for later constant-potential DFT. Its energies and
forces are never used as training labels, and its vacuum trajectories are not a model of the electrolyte
dynamics in this project.

## Source

| item | value |
|---|---|
| record | Materials Cloud Archive `71861-xfe58` |
| DOI | 10.24435/materialscloud:va-hx |
| title | Low-index mesoscopic surface reconstructions of Au surfaces using Bayesian force fields |
| authors | Owen et al., published 2024-02-29 |
| file taken | `Au_training.zip`, 9,975,527 bytes |
| md5 (record) | 35886d1f5eea0d8488be7fc2bd0c471c |
| md5 (downloaded) | 35886d1f5eea0d8488be7fc2bd0c471c (match) |
| sha256 (downloaded) | 90ac6825d59d345438babc1f7f22c682731ef10004f29c8e2db0975e67e4acb0 |
| also taken | `README.md` (md5 b573a13f…), `files_description.md` (md5 2c18318c…) |
| NOT taken | `Au_MD_simulations.zip`, 2.58 GB of the authors' own LAMMPS runs; not needed to use the potential |
| downloaded | 2026-10-02 from `https://archive.materialscloud.org/api/records/71861-xfe58/files/<name>/content` |
| licence field in the record | not stated in the API metadata; cite the paper and the archive record when results are used |

## Contents of `Au_training/`

| file | size | what it is |
|---|---|---|
| `lmp_t0.0001_no_bulk_vac_fix3.flare` | 394,306 B, 3,264 lines | the **mapped coefficient file** for LAMMPS `pair_style flare`; dated 2022-12-15, contributor "cjo" |
| `Au_master.xyz` | 36.9 MB, 2,966 frames | the on-the-fly training set (extxyz with energy, forces, stress) |
| `Au_master.yaml` | 1,464 B | the FLARE on-the-fly configuration the model was trained with |
| `to_lmp.py`, `fix.py` | small | the authors' scripts to rescale the energy noise and map the SGP to the coefficient file |

## Model settings, read from `Au_master.yaml` and the coefficient-file header

| setting | value |
|---|---|
| kernel | NormalizedDotProduct, sigma 2.0, power 2 |
| descriptor | B2 (ACE-like), n_max 8, l_max 4, Chebyshev radial basis, quadratic cutoff function |
| cutoff | **6.0 Å** |
| species | Au only (Z = 79); single-atom energy 0 |
| noise used in training | energy 0.1, forces 0.05, stress 0.005 (then rescaled per `fix.py`: energy noise 0.003 eV per atom) |
| coefficient-file header | `2 / chebyshev / 1 8 4 16290 / quadratic / 6.00` → 1 species, n_max 8, l_max 4, 16,290 coefficients |
| units | eV and Å (training energies ≈ −2.95 eV/atom); LAMMPS `units metal` |

## How it is called

Only through LAMMPS, with FLARE's pair style compiled in (`flare/lammps_plugins`, applied to a LAMMPS source
tree built inside this project):

```
units        metal
newton       on
pair_style   flare
pair_coeff   * * lmp_t0.0001_no_bulk_vac_fix3.flare
```

The Python `flare` package is **not** installed: the record provides the mapped coefficient file and the
training set, not the sparse-GP JSON that the Python `SGP_Calculator` would load, so the LAMMPS route is the
only one the released files support.

## Rules for its use here

* All software for it lives under `rough_sampling_v1/env/`; nothing outside the project is modified.
* Fixed atoms (bottom two layers) are not integrated and have their velocities zeroed; the thermostat acts on
  the movable atoms only.
* Its energies and forces are not added to the CP-DFT training labels.
* A low-cost behaviour check on the existing flat, step, adatom, island and pit cells is run and recorded
  before any parent surface is evolved; no new DFT is run for that check.
* The LAMMPS version, the FLARE plugin commit and the build options are recorded in `env/BUILD.md` when the
  build completes.
