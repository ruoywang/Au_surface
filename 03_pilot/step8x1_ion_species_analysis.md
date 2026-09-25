# Ion species reconstruction on Step-8x1 (post-processing only, no new DFT)

Reconstructs n₊(r)/n₋(r) from the already-computed Step-8x1_muref fields (PHI, SION,
RHOION), using v2 doc §7.2's constitutive relation:

```
u = e·psi/(kB·T),  D = 1 + (2·n_bulk/n_max)·(cosh(u)-1)
n_plus  = SION·n_bulk·exp(-u)/D
n_minus = SION·n_bulk·exp(+u)/D
psi = PHI - PHI_bulk
```

`n_max = 1/d_ion³`, `d_ion = 2^(5/6)·R_ION` (D_ION not set in INCAR → this default,
confirmed in `solvation.F:3213-3215`). T = 298.15 K.

## Validation against RHOION (required before trusting anything below)

Reconstructed n₊−n₋ vs. the volume-corrected RHOION field, in the bulk window
(z=21.11–31.06 Å): **max difference 1.5×10⁻⁵ e/Å³ — the same order as RHOION's own
bulk noise floor.** No sign flip on psi was needed (PHI is already in the convention
this formula expects). This confirms the constitutive relation, n_max, and n_bulk are
right, not assumed.

## The central caveat: this run is at ΔU=0, so a flat result is expected, not a finding

Step-8x1_muref was computed at TARGETMU = mu0 = -4.9071 eV — the *neutral/PZC-
equivalent* reference point, same as T's own reference. **There is no electrode
charge and therefore no double-layer driving force in this data.** Everything below
should be read as "the reconstruction pipeline works," not as "steps don't enrich
anions" — that question needs a Step point at nonzero ΔU, which does not exist yet.

## 1. Normal (z) profiles

`step8x1_ion_zprofile.png`: n₊/n_bulk and n₋/n_bulk both track SION closely and
converge to 1.000 ± 0.011 (std) in the SION=1 bulk window — i.e. bulk salt
concentration is fully recovered there, and n₊≈n₋ throughout (no net charge, as
expected at ΔU=0). The small (~1-2%) n₊/n₋ splitting right at the accessibility
transition (z≈18-21 Å) is a real, physically sensible Boltzmann response to the local
field gradient created by the SION transition itself, not noise.

## 2. Cross-step section

`step8x1_crossstep_map.png`: psi, SION, and n₋/n_bulk averaged along the step
direction, shown vs. (x across-step, z). SION's accessible-region boundary visibly
dips near the step foot (real atomic-scale corrugation in the excluded-volume
boundary, not a smooth mean-field cutoff) — the elevated-strip side needs the
accessible boundary to start higher than the bare-terrace side does. n₋/n_bulk tracks
SION's shape almost exactly, again because there is no differential ion response to
show at ΔU=0.

## 3. Step-adjacent vs. terrace enrichment

Near-edge (within 3 Å of either step edge) vs. far-from-edge (terrace interior),
accessible-volume-weighted, in the bulk-z window:

| region | n₋/n_bulk |
|---|---|
| near-edge | 1.0002 |
| far-from-edge (terrace interior) | 0.9999 |
| ratio | 1.0003 |

Essentially no difference — consistent with ΔU=0 giving no enrichment signal to find.
**Separately**: Step-8x1's terrace is only ~11.8 Å wide total (10.2 Å/side by the
correct perpendicular measurement); the "far-from-edge" region above is 62% of the
cell width but is a genuinely narrow periodic strip, not an isolated-step limit. Any
future enrichment ratio computed here should be reported as a narrow-periodic-array
conditional value, not generalized to an isolated step.

## What this means for the next DFT decision

The pipeline (field reconstruction, cross-step mapping, region-based enrichment) is
now validated and ready. But **the specific question "does a step enrich anions?"
cannot be answered from data that only exists at the PZC** — that's a property of the
charged double layer. Before spending budget on Step-16x1 (a width-series repeat,
still at ΔU=0 by default), it may be more directly useful to compute a Step point
(8x1 or 16x1) at a **nonzero ΔU** first — that is what would actually exercise this
analysis pipeline on a non-trivial signal. No DFT submitted this round either way;
this is a recommendation, not an action taken.
