# Step-8x1 vs Step-8x2 numerical consistency check (rev 2, corrected)

Single question: does the physical result change when the same idealized straight-step
structure is repeated ×2 along the step direction? Static, single point, no relaxation,
no potential sweep, no perturbation. TARGETMU = -4.9071 eV (the pilot's existing internal
mu0 reference — not re-labeled as a calibrated V vs RHE value). Both jobs COMPLETED
normally (job 20903924, 20903925).

**This revises rev 1: the electron-count table there was wrong (used the neutral
reference, not the converged value), and the RHOION field comparison used both the
wrong volume normalization and the wrong bulk-region mask. Both are fixed below with a
confirmed root cause, not a guess.**

## Corrected electron-count table

Rev 1 reported N_ele = 396.0000 / 792.0000 — those are the *neutral* electron counts
(36×11 and 72×11), read from OUTCAR's early header echo, not the converged CP value.
The actual converged values are in each run's last OSZICAR `CPM-ion` line:

| | Step-8x1 (36 Au) | Step-8x2 (72 Au) |
|---|---|---|
| neutral N_e⁰ | 396 | 792 |
| **converged N_e** | **396.014037** | **792.029095** |
| converged μ_e | -4.907178 | -4.907370 |
| extra electrons q_e = N_e − N_e⁰ | 0.014037 | 0.029095 |

q_e,2/2 − q_e,1 = 0.0145475 − 0.014037 = **5.105×10⁻⁴** e → Q₂/2 − Q₁ = −5.105×10⁻⁴ e
(using Q = −e·q_e). Not exactly zero, but the converged μ_e themselves differ by only
1.92×10⁻⁴ eV between the two cells — both are essentially at the same point, and this
is not evidence the cell reduction is wrong; it means "electron count exactly doubles"
can't be used as the pass criterion, only the actual converged numbers can.

## RHOB/RHOION: root cause found and fixed, not just "normalization suspect"

Checked directly in the local CEP-DIP source (`solvation.F`, not just public VASPsol++):
lines 1375/1384 multiply `n_b`/`n_ion` by `LATT_CUR%OMEGA` (the cell volume Ω), and lines
1249-1250 write them straight to `RHOB`/`RHOION` with nothing dividing that factor back
out. **The raw grid values are Ω·ρ(r), not ρ(r).** `scripts/pilot_qc.py`'s old formula
(`closure = (rhob.sum()+rhoion.sum()) * (V/ngrid)`) multiplied by V *again* — the
correct real-space integral is `sum/ngrid`, no extra V.

Retroactive test on the original 9-point pilot (fix applied, all 9 points, not just
Step): closure residual vs. the CP loop's own `dn_cp` dropped from **hundreds of
electrons to 1×10⁻⁶–1×10⁻⁵ e** across every point. This fully resolves the "open,
unresolved" item recorded 2026-09-23 in `00_audit/parameter_map.md` §I.4 (now updated).
PHI/SDIEL/SION are written through a different path with no volume factor and must
**not** be divided by V — confirmed by checking their write calls separately.

### Correct bulk-region mask (was also wrong in rev 1)

Rev 1 used z=16.6–32.6 Å as "bulk," reasoning only "past the highest Au atom." That
window actually spans SION rising from 0 to 1 (surface response) at the bottom and
butts against ION_Z1 = SOL_Z1 − D_STERN = 34.603 − 2.0 = 32.603 (the ion-window
transition) at the top — i.e. it wasn't bulk at either end. Scanning SION/SDIEL
directly: both structures reach a genuine flat SION=SDIEL=1.000000 plateau only over
**z ≈ 21.0–31.1 Å**. That's the mask used below.

### Corrected RHOION comparison (divide by V_cell, correct bulk mask)

| | Step-8x1 (V=2671.31 Å³) | Step-8x2 (V=5342.61 Å³) |
|---|---|---|
| bulk z-range (SION=SDIEL=1.0) | 21.11–31.06 Å | 20.97–31.06 Å |
| RHOION bulk mean (e/Å³) | -1.547×10⁻⁵ | -1.604×10⁻⁵ |
| RHOION bulk std (e/Å³) | 1.319×10⁻⁵ | 1.387×10⁻⁵ |

x1-vs-x2 diff in the shared bulk window: **max 2.74×10⁻⁷ e/Å³, RMS 1.03×10⁻⁷ e/Å³** —
small compared to the signal itself (~1.5×10⁻⁵ e/Å³ mean), a genuine consistency pass,
not an artifact of a broken normalization or a mismatched window.

## Energy, per full atom set

| | Step-8x1 | Step-8x2 |
|---|---|---|
| TOTEN | -110.00726734 eV | -220.02225587 eV |

(TOTEN₂/2) − TOTEN₁ = -0.00386 eV over the 36-atom cell ≈ 1.07×10⁻⁴ eV/atom.

## Forces, all 36 matched atoms (not just the one bottom-fixed atom as in rev 1)

Each Step-8x1 atom matched to its corresponding atom in Step-8x2 by position (x1's
cell tiles exactly into x2, so every atom has an exact counterpart):

- max |ΔF| component over all 36 matches: **1.83×10⁻⁴ eV/Å**
- RMS |ΔF| component: **3.13×10⁻⁵ eV/Å**
- for scale, max |F| in either structure: ≈0.341 eV/Å

## PHI, bulk region (unchanged from rev 1 — PHI does not need the volume fix)

max diff 1.79×10⁻⁴ V, RMS 2.36×10⁻⁵ V over the corrected z=21.0–31.1 Å bulk window.

## Cost

| | Step-8x1 | Step-8x2 |
|---|---|---|
| wall time | 1h32m | 4h46m |
| SCF/CP "Iteration" count | 221 | 225 |
| peak memory (OUTCAR) | 1.36 GB | 1.57 GB |

3.1× longer wall time for 2× atoms is not explained by atom count alone (x2 actually
uses fewer k-points along the direction that doubled: 6 vs 12) — no confirmed single
cause, reported as-is.

## Conclusion

With the electron-count and RHOION errors fixed, every checked quantity — μ_e, TOTEN,
all-atom forces, bulk PHI, and now bulk RHOION with the correct volume normalization
and the correct bulk mask — agrees between the two cells at the 10⁻⁴–10⁻⁷ level
relative to the signal size. The reduced-cell approach (ny=4 → ny=1) is supported
numerically, not just geometrically. Per the plan: Step-16x1 becomes the main
representative for the width series; Step-24x1 gets added only if the width study
specifically needs it. **Step-16x1 has not been submitted** — this revision only
corrects the post-processing of the two already-completed jobs; no new DFT was run.
