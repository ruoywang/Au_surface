# Mechanism check, rev 2 — corrected statistics, K_ij decomposition

> **Interpretation superseded 2026-09-26** (`step16x1_width_dependence.md`). The numbers
> below are unchanged and reproduce, but the near/far statistic they rest on averages
> over a bimodal across-step profile: the raised terrace responds like the flat surface
> and the *lower* terrace is enhanced 1.2–1.6×, with the countercharge to the edge charge
> displaced ~4.7 Å onto the lower terrace. §5's "different point on the nonlinear
> occupation curve" hypothesis is no longer needed. The K_ij result (response is the
> potential term, ΔK_S ≈ 0) stands.

Rev 1's claims "excludes the geometric-weighting explanation" and "excludes the
electrostatic explanation" are **retracted as overreaching** — both used inconsistent
weighting between the plain Δψ average and the S_ion-weighted K_D metric, and rev 1's
"per-column ratio" statistic wasn't actually the same quantity as the main script's
whole-region ratio. Fixed below; conclusions revised where the fix changed them.

**Three statistical bugs fixed** (identified before trusting rev 1's numbers further):
1. `edge_dist` used the raw sheared `a1`, not the true perpendicular `L_perp` — rev 1's
   "2/3/4 Å" bands were actually ~1.73/2.60/3.46 Å. Now uses the same L_perp-based
   distance as the main analysis script (confirmed: band=3Å here now matches the main
   report's 0.868 ratio exactly).
2. Edge1/edge2 "far" was `~near` for that edge alone, silently including the other
   edge's near-region. Fixed: far = excludes both edges (the common terrace region).
3. Accessibility-relative sampling blended a per-column-averaged statistic into what
   was reported as "the" result. Fixed: primary metric is the whole-region Σn/ΣS
   ratio (same convention as sections 1-2); the per-column average is reported
   separately and labeled as a different statistic (they happen to coincide here to
   4 decimals, which is a property of this particular field, not guaranteed generally).

## 1. Bandwidth sensitivity (corrected distance)

| band | ΔK_near | ΔK_far | ratio |
|---|---|---|---|
| 2 Å | 0.1168 | 0.1365 | 0.855 |
| 3 Å | 0.1213 | 0.1397 | 0.868 |
| 4 Å | 0.1258 | 0.1404 | 0.896 |

Still robustly <1 across bandwidths, magnitudes barely changed from rev 1 (the
distance-formula bug was small in absolute terms here, but is fixed for correctness).

## 2. Edge1 vs. edge2 (far now correctly excludes both edges)

| | ΔK_near | ΔK_far (common terrace) | ratio |
|---|---|---|---|
| edge1 (u=0) | 0.1185 | 0.1397 | 0.848 |
| edge2 (u=0.5) | 0.1240 | 0.1397 | 0.887 |

Both edges individually <1 against the same, correctly-shared far reference.

## 3. Accessibility-relative sampling (primary metric fixed)

Whole-region ratio (Σn/ΣS over all near/far columns together), per-column window
above the local SION=0.5 boundary (corrugating 18.25–20.64 Å): **ΔK_near=0.1129,
ΔK_far=0.1381, ratio=0.817.** The per-column equal-weighted average gives the same
numbers here (a property of this field's fairly uniform per-column accessible volume
in each region, not guaranteed in general — reported separately per the correction).

**All three checks still show ratio < 1, in the 0.82–0.90 range.** This remains a
real, robust empirical trend across every resampling tried — not yet a validated
mechanism, which is what the next section addresses properly.

## 4. K_ij decomposition — the actual mechanism split (new this revision)

`K_ij` = occupation function from potential state *i*'s φ, weighted by state *j*'s
SION. K00/K11 are the real self-consistent points; K10/K01 are post-processing
cross-combinations for this decomposition only — not new physical states.

| region | K00 | K10 | K01 | K11 | ΔK_φ | ΔK_S | ΔK_φ+ΔK_S (=K11−K00, exact identity) |
|---|---|---|---|---|---|---|---|
| near-edge | 0.9830 | 1.1037 | 0.9829 | 1.1043 | **+0.1210** | +0.0002 | 0.1213 |
| far-from-edge | 0.9785 | 1.1177 | 0.9784 | 1.1182 | **+0.1396** | +0.0002 | 0.1397 |

**The accessibility/geometry term (ΔK_S) is negligible and essentially identical in
both regions (+0.0002).** Essentially all of the response — and all of the
near/far *difference* in response — comes from the potential-driven occupation term
(ΔK_φ). This rules out "the step changes ion accessibility/packing differently than
the terrace" as the explanation. It does **not** rule out an electrostatic origin —
it relocates the question to a sharper one: why does the *same* potential
perturbation produce a smaller occupation response near the edge?

## 5. Reconciling with the potential-difference map

Section 5 (unweighted spatial Δψ) still shows near-edge and far-from-edge regions
experience essentially the same *raw* potential shift (−0.0048 V vs. −0.0047 V). But
K_ij shows their *occupation response* to that shift differs by ~13-15%. Since the
occupation function g(u) = exp(−u)/D(u) is nonlinear, a given Δu produces a different
Δg depending on the local starting point u₀ (via the crowding denominator D(u₀)) — and
K00 already differs slightly between the regions (0.9830 near vs. 0.9785 far) before
any charging. **This is now a specific, checkable hypothesis** — the differential
response comes from each region sitting at a different point on the same nonlinear
occupation curve, not from a different local potential shift or a separate
accessibility mechanism — rather than an open-ended "maybe crowding, maybe SION
gradient" guess.

## What this does and doesn't support, revised

**Can write**: the ~13-20% weaker near-edge anion response is a robust empirical
finding across bandwidth, individual edges, and accessibility-relative resampling; it
is driven almost entirely by the potential-occupation term rather than by
accessibility/geometry change; and it is consistent with the two regions sitting at
different points on the same nonlinear occupation curve, given that the raw potential
shift itself is essentially uniform between them.

**Still cannot write**: that this is confirmed as *the* mechanism (vs. correlation);
that it generalizes beyond this narrow (10.2 Å/side), single-width, two-point
comparison; or anything about a real electrode's Cl⁻ behavior.

## Next: Step-16x1 at both potentials — done, see `step16x1_width_dependence.md`

Step-16x1 was run at both potentials (same construction, only nx doubled; k-points
1×12×1 by the same reciprocal-density rule). The Accurate/ALGO=Normal dUp02 run timed
out (charge sloshing, see `00_audit/parameter_map.md` §K.6); both points were then
rerun with the adopted production config (§K.8) and converged in 43/40 min. Result:
R₁₆ = 0.936 (3 Å band, same K metric) vs R₈ = 0.868 — the sign is robust, but the
16×1 across-step profile shows that the 8×1 "far" strip sits on an *enhanced* shoulder
3–5 Å from the edges rather than on a converged terrace, so the near/far ratio is not
a width-independent number. The 2026-09-26 report re-frames the comparison on a
window-independent surface-excess footing and against the flat T pair.
