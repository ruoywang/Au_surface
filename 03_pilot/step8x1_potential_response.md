# Step-8x1: potential response, muref vs. dU=+0.2V

Two points, identical geometry/cell/window/k-points/ISMEAR/ENCUT/FERMICONVERGE — only
TARGETMU differs. Both COMPLETED normally (job 20903924 muref, job 20907809 dUp02).
Both validated against RHOION on the full 3D grid at relL2=0.0005 (same as the T
references), using the corrected pipeline (no artificial zero-point shift, 3D-first
reconstruction, corrected K_D and near-surface region definitions).

## Overall charging response

| | Step-8x1_muref (TARGETMU=-4.9071) | Step-8x1_dUp02 (TARGETMU=-5.1071) |
|---|---|---|
| converged μ_e | -4.907178 | -5.106659 |
| \|μ_e - target\| | 0.000078 eV | 0.000441 eV (both within FERMICONVERGE=0.001) |
| N_e (converged) | 396.014037 | 395.923404 |
| q_e = N_e - 396 | +0.014037 | -0.076596 |
| Q_phys = -e·q_e | -0.014037 e (weakly negative electrode) | +0.076596 e (positive electrode) |
| wall time | 1h32m | 1h38m |

Lowering the target electron chemical potential by 0.2 eV moved the electrode from
weakly negative to clearly positive, as expected — the electrode lost a net 0.0906 e
(0.076596 - (-0.014037)) relative to the reference point.

## Overall (bulk-window) ion response — increased in the expected direction

| | n_anion/n_bulk | n_cation/n_bulk |
|---|---|---|
| Step-8x1_muref | 0.9873 (slightly depleted) | 1.0129 (slightly enriched) |
| Step-8x1_dUp02 | **1.0698** (clearly enriched) | **0.9347** (clearly depleted) |

The overall anion population in the bulk-accessible region rose substantially (0.987 →
1.070) as the electrode became more positive — the expected direction, and now a much
larger, more clearly resolved signal than the muref point's weak +1.3%/-1.3% split.

## Near-surface region, step-edge vs. terrace-interior selectivity — did NOT simply scale up

Same near-surface region definition (metal-top to +15 Å, excludes the far artificial
window) at both potentials; same near-edge (<3 Å) / far-from-edge split (terrace width
10.18 Å/side, true perpendicular measurement):

| | K_D near-edge | K_D far-from-edge | S(U) = K_near/K_far |
|---|---|---|---|
| Step-8x1_muref | 0.9830 | 0.9785 | **1.0046** |
| Step-8x1_dUp02 | 1.1042 | 1.1182 | **0.9875** |

**The overall anion increase did not carry through to a larger step-vs-terrace
selectivity — S(U) actually flipped sides.** At the weak/reference potential, the step
edge was very slightly *less* anion-depleted than the terrace interior (S>1). At the
more positive potential, with anions clearly enriched everywhere, the step edge became
very slightly *less* enriched than the terrace interior (S<1). Both effects are small
(0.46% and 1.25% respectively) and in a narrow, non-isolated periodic terrace (far
fraction still only 41% of the cell width) — this is not evidence of a robust,
sign-reversing step effect, but it is a real result of the corrected reconstruction,
and it directly illustrates the point raised before running this: **a bigger overall
charging response does not automatically mean a bigger (or even same-sign) local
step/terrace selectivity difference.** That has to be checked, not assumed.

## What this does and doesn't support

- The reconstruction pipeline is now validated across four independent points (T_dUm02,
  T_dUp02, Step-8x1_muref, Step-8x1_dUp02), two structures, both charge signs, all at
  relL2=0.0005 on the full 3D grid.
- The overall (bulk-region) anion-enrichment response to a more positive potential is
  clear, large, and in the expected direction.
- The step-vs-terrace *local selectivity* signal is small, and its sign is not fixed
  across potential — this needs more potential points and/or a wider terrace (to get
  clear of the narrow-periodic-array caveat) before it can support any claim about
  whether steps preferentially enrich anions relative to flat terraces under bias.

## Stopping here per plan

No further potentials, structures, or parameter scans added this round. Both figures
(`Step-8x1_dUp02_ion_zprofile.png`, `Step-8x1_dUp02_crossstep_map.png`) and the
Step-8x1_muref figures are already generated from the corrected pipeline.
