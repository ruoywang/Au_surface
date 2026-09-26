# Mechanism check: is the near-edge response really weaker, or an artifact?

Follows up on `step8x1_potential_response.md`'s S(U) crossing 1 between the two
potentials. Tests the simpler alternative explanation directly: **the step-edge
region's anion response to charging may just be weaker in magnitude than the terrace
interior's, at both potentials — which alone would make K_near/K_far cross 1 as the
overall sign of (K−1) flips, with no new physical mechanism required.** No new DFT;
all of this reuses the two already-completed, already-validated points.

**Also fixed while doing this**: the reconstruction used T=298.15 K; the local code's
actual default (`solvation.F:1703`, not overridden by INCAR) is **298.0 K**. Redone at
298.0 K, relL2 on all 4 validation points drops from ~5×10⁻⁴ to ~2×10⁻⁵ (20×) —
essentially floating-point/SCF precision. K_D values themselves barely move (the
298.15 K error was too small to matter for these particular numbers), so no
conclusion below changes — but the reconstruction is now exact, not 0.05%-off.

## 1. Same-explanation reframing: response increments, not point values

| | K_D near-edge | K_D far-from-edge |
|---|---|---|
| Step-8x1_muref | 0.9830 | 0.9785 |
| Step-8x1_dUp02 | 1.1043 | 1.1182 |
| ΔK (dUp02 − muref) | **0.1213** | **0.1397** |

ΔK_near/ΔK_far = 0.1213/0.1397 = **0.868**. Finite-difference susceptibility over the
0.2 V gap: χ_near = 0.607 V⁻¹, χ_far = 0.699 V⁻¹. **The near-edge region's anion
response to this charging step is about 13% smaller than the terrace interior's** —
this is what actually crosses K=1 at different potentials for the two regions, not a
sign-reversing step mechanism. Both regions' absolute occupation still only differ by
about 1% at either single potential — the 13% is a difference in *response size*, not
in absolute concentration.

## 2. Robustness checks (all reuse existing fields, no new DFT)

**Near-edge bandwidth (2/3/4 Å), both edges lumped:**

| band | ΔK_near | ΔK_far | ratio |
|---|---|---|---|
| 2 Å | 0.1159 | 0.1354 | 0.856 |
| 3 Å | 0.1191 | 0.1386 | 0.859 |
| 4 Å | 0.1235 | 0.1404 | 0.880 |

Stable across bandwidths — not a threshold artifact.

**Edge1 vs. edge2, separately (3 Å band, not lumped, still no A/B label):**

| | ΔK_near | ΔK_far | ratio |
|---|---|---|---|
| edge1 (u=0) | 0.1168 | 0.1330 | 0.878 |
| edge2 (u=0.5) | 0.1213 | 0.1316 | 0.922 |

Both edges individually show the same direction — not an artifact of one edge
dominating a lumped average.

**Accessibility-relative sampling** (per-column local SION=0.5 boundary — which
corrugates between z=18.25 and 20.64 Å near the step, confirming the step really does
locally raise/lower the accessible boundary — with a fixed 10 Å window measured
*above that local boundary*, not a fixed absolute height):

ΔK_near=0.1102, ΔK_far=0.1364, **ratio=0.807** — same direction, if anything slightly
stronger. This directly rules out "it's a fixed-height sampling artifact that ignores
the corrugated accessible boundary."

**Conclusion: the ~13-20% weaker near-edge response survives every resampling tried.**
This is now a real, if still small, model-internal result — not a sign-reversal
mechanism, and not a geometric sampling artifact.

## 3. The genuinely open part: it's not simply a weaker local potential shift

Compared psi_dUp02 − psi_muref at the same spatial positions in the near-surface
region:

| | mean Δpsi (V) |
|---|---|
| near-edge columns | −0.0048 |
| far-from-edge columns | −0.0047 |
| difference | −0.0001 |

**The driving potential shift is essentially identical near the edge and far from it**
— the simplest explanation (near-edge region "feels" a smaller potential change) does
**not** hold. The weaker ion-response increment must come from something else: the
nonlinear crowding factor D(u) in the Boltzmann relation, a difference in how SION's
own value or gradient near the corrugated boundary interacts with the response, or
some other coupling not yet identified. This is reported as an open mechanistic
question, not resolved here — the data support "near-edge response is weaker" as a
robust empirical fact, not yet a mechanistic explanation for *why*.

## What this does and doesn't support

**Can write**: at this potential step, the step-edge region's anion occupation
responds to charging by about 13-20% less than the terrace interior does, in a fixed,
narrow (10.2 Å/side), periodic Au(111) step array under this continuum electrolyte
model — robust across near-edge bandwidth, both edges individually, and an
accessibility-relative resampling. The driving electrostatic potential shift itself is
not measurably different between the two regions, so a simple "weaker local field"
explanation is ruled out.

**Cannot yet write**: that Au steps generally suppress anion response, that this
reflects a real electronic/geometric mechanism (vs. a narrow-terrace-specific
artifact), or anything about Cl⁻ specifically. The terrace is still only 10.2 Å/side —
a genuinely isolated-step comparison needs a wider terrace (Step-16x1) at the same
bias before generalizing.

## Next step, not executed

Per plan: if this trend holds up under the checks above (it does), the next DFT
budget should go to Step-16x1 at the same ΔU=+0.2V, to test whether the ~13-20%
weaker-response finding is a real local effect or specific to this narrow periodic
array. **Not submitted this round.**
