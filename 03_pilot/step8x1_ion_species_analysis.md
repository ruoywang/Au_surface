# Ion species reconstruction (rev 2 — corrected, validated on charged references)

Rev 1 is retracted where it disagreed with this. Four confirmed bugs fixed:

1. **Zero point**: rev 1 subtracted a finite SION=1-window average from PHI before
   using it as the Boltzmann driving potential. Wrong — validated below that raw PHI,
   with **no shift at all**, reconstructs RHOION to machine-precision-adjacent accuracy
   on a real charged reference. The solver's own PHI already carries the correct
   physical zero point; subtracting a window average was an unjustified extra step
   that (at Step-8x1's weak charge) partly washed out the real signal.
2. **Species identity**: not assumed from the v2 doc's variable names. Empirically
   anchored: the `exp(-u)` species matches `RHOION/V_cell` as `n_A − n_B` (confirmed to
   <2.1×10⁻⁷ e/Å³ on the full 3D grid), **and** is the species enriched near a
   positively-charged electrode — which basic electrostatic screening requires to be
   the anion. So: `n_A = n_anion`, `n_B = n_cation`. (The doc's own "n_plus"/"n_minus"
   labels do not map to cation/anion — do not assume that naming.)
3. **3D-first reconstruction**: the Boltzmann relation is nonlinear, so `⟨n(φ,S)⟩ ≠
   n(⟨φ⟩,⟨S⟩)`. Rev 1 averaged PHI/SION to 1D profiles before applying the formula.
   Fixed: reconstruct on the full (ngz,ngy,ngx) grid, average afterward.
4. **Enrichment integral**: rev 1's `weighted_mean` multiplied by SION a second time
   in the numerator (n_anion already includes it). Fixed: `K_D = ∫n_anion / (n_bulk·
   ∫S_ion)` over a region, no extra SION factor. Edge distance now uses proper
   minimum-image wrapping and the true perpendicular terrace width (10.18 Å/side for
   Step-8x1), not the raw sheared lattice vector.

## Validation: rigorous, on real charged signals — not a near-zero point

Retracted: "Step-8x1_muref is the step's PZC, so near-zero response is expected by
construction." **That was wrong.** Its converged N_e = 396.014037 vs. neutral 396 →
q_e = +0.014037, Q_phys = −0.014037 e. Small, but not zero — and a shared reference μ
(from the flat-plate T) is not automatically the step's own zero-charge potential
either. Total-cell neutrality does not imply local charge or potential are zero.

Validated instead on **T_dUm02** (q_e=−0.164483, TARGETMU=−5.1071) and **T_dUp02**
(q_e=+0.166912, TARGETMU=−4.7071) — real, oppositely-signed, non-trivial charge states:

| reference | max\|diff\| (e/Å³) | relL2, full cell | relL2, bulk (SION=1) | relL2, interface | n_anion/n_bulk (bulk) | n_cation/n_bulk (bulk) |
|---|---|---|---|---|---|---|
| T_dUm02 (Q_phys=+0.164e) | 2.07×10⁻⁷ | 0.0005 | 0.0005 | 0.0005 | 1.0806 (enriched) | 0.9265 (depleted) |
| T_dUp02 (Q_phys=−0.167e) | 2.08×10⁻⁷ | 0.0005 | 0.0005 | 0.0005 | 0.9259 (depleted) | 1.0813 (enriched) |

Relative L2 error is **0.05% in every region tested**, on a signal that is not near
zero — a real, discriminating pass, not a threshold that happens to exceed the noise
floor. Both charge signs give the physically required direction (positive electrode →
anion excess; negative electrode → cation excess), which independently confirms the
species assignment. **This is the actual pass criterion the previous "1.5×10⁻⁵ ≈ signal
scale" check in rev 1 could not provide.**

Also confirms the user's warning was correct: even in the "SION=1" window, T_dUm02
shows n_anion/n_bulk=1.08 — **the fully-accessible region has not fully recovered bulk
electroneutrality; a real double-layer tail extends into it.** "Completely accessible"
≠ "reservoir-equilibrated." This is now stated explicitly rather than assumed away.

## Corrected Step-8x1_muref result (still weak charge, but a real, honest signal)

q_e=+0.014037 (electron surplus → Q_phys=−0.014e, weakly negatively charged):

| | n_anion/n_bulk | n_cation/n_bulk |
|---|---|---|
| SION=1 bulk window | 0.9873 | 1.0129 |

Correct direction for a weakly negative electrode (anion depleted, cation enriched) —
small because the charge is small, not zero because of a methodology artifact this time.

**Near-edge vs. far-from-edge K_D** (anion, corrected formula, true 10.18 Å/side
terrace width, proper periodic edge distance):

| region | K_D near-edge (<3 Å) | K_D far-from-edge | ratio |
|---|---|---|---|
| deep-bulk reservoir window (diagnostic) | 0.9877 | 0.9866 | 1.0012 |
| near-interface EDL region (physical signal) | 0.9779 | 0.9650 | **1.0134** |

A small (~1.3%) but directionally real signal: anion depletion is *less severe* near
the step edge than in the terrace interior, in the region that actually carries the
double-layer response. Far weaker than a definitive claim, but no longer an artifact
of an over-aggressive zero-point correction washing the signal to ~0. Step-8x1's
terrace (far-from-edge fraction 41% of the cell width with a 3 Å near-edge band) is
still a narrow periodic array, not an isolated-step limit — this ratio is reported as
that conditional value, not generalized further.

## Recommendation, unchanged in substance from before

A charged Step point is the right next experiment — now for a better reason: the
corrected pipeline shows Step-8x1_muref's existing (weak) charge already produces a
small, physically consistent, directionally real signal, and a real ΔU=+0.2 V point
(TARGETMU=−5.1071 eV, same geometry/cell/window as Step-8x1_muref, no relaxation, no
expansion to other structures) would give a much larger, more clearly resolved version
of exactly this signal to check against. **Not submitted this round** — reporting the
corrected validation first, as requested.
