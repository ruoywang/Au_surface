# Step-8x1 vs Step-8x2 numerical consistency check

Single question: does the physical result change when the same idealized straight-step
structure is repeated ×2 along the step direction? Static, single point, no relaxation,
no potential sweep, no perturbation. TARGETMU = -4.9071 eV (the pilot's existing internal
mu0 reference — not re-labeled as a calibrated V vs RHE value). Both jobs COMPLETED
normally (job 20903924, 20903925).

| item | Step-8x1 (36 Au) | Step-8x2 (72 Au) | comparison |
|---|---|---|---|
| target μ_e | -4.9071 | -4.9071 | same target |
| converged E-fermi | -4.9072 | -4.9074 | Δ(target) = 0.0001 / 0.0003 eV — both inside the 0.001 eV FERMICONVERGE goal; the two runs differ from each other by 0.0002 eV |
| NELECT | 396.0000 | 792.0000 | exactly ×2 (792/2 − 396 = 0 e) |
| TOTEN | -110.00726734 eV | -220.02225587 eV | (TOTEN₂/2) − TOTEN₁ = -0.00386 eV over the 36-atom cell ≈ 1.07×10⁻⁴ eV/atom |
| force, bottom-fixed atom (0,0,5) | (0.007774, -0.005445, -0.068911) eV/Å | (0.007766, -0.005439, -0.069023) eV/Å | max component diff ≈ 1.1×10⁻⁴ eV/Å |
| PHI(z), bulk region (z=16.6–32.6 Å, clears the strip top at 14.6 Å) | mean 0.000742 V, std 9.97×10⁻⁴ V | mean 0.000750 V, std 1.02×10⁻³ V | diff: max 1.79×10⁻⁴ V, RMS 2.36×10⁻⁵ V — consistent |
| RHOION(z), same bulk region | mean -3.91×10⁻² , std 3.72×10⁻² | (interpolated to x1 grid) | diff: max 0.174, RMS 0.058 — **not clean**, see caveat below |
| SCF/CP effort ("Iteration" lines in OUTCAR) | 221 | 225 | essentially equal effort per solve |
| wall time | 1h32m | 4h46m | 3.1× longer for 2× atoms — not explained by atom count alone (x2 actually has fewer k-points in one direction, 6 vs 12); no confirmed single cause |
| peak memory (OUTCAR-reported) | 1.36 GB | 1.57 GB | |

## Reading

**E-fermi, NELECT, TOTEN/atom, forces, and bulk PHI all agree at the 10⁻⁴–10⁻³ level or
better between the two cells.** This is exactly what "same physical structure, cell
repeated ×2" should produce — the ×2 along-step repeat carries no new physics once
k-points are matched by reciprocal-space density rather than copied unchanged.

**RHOION's bulk-region planar average is not itself flat** (std ≈ its own mean, in a
region well clear of all metal) — the x1-vs-x2 difference there (0.058 RMS) is
comparable to that non-flatness, not a clean independent number. This is the same
open, previously-flagged RHOION field-normalization issue from the original 9-point
pilot (documented in `00_audit/parameter_map.md` §I), not a new problem introduced by
the cell reduction. It does not affect NELECT/TOTEN (which come from the CP loop
directly, not from integrating this field), so the electron-count-based numbers above
stay trustworthy; RHOION-derived ion-concentration/enrichment quantities specifically
still should not be trusted until that issue is traced.

## Conclusion

The reduced-cell approach (ny=4 → ny=1) checks out numerically, not just geometrically.
Per the plan: Step-16x1 becomes the main representative for the width series;
Step-24x1 gets added only if the width study specifically needs it. No expansion to
islands/pits/vicinal slabs, no potential sweep, no relaxation triggered by this result.
