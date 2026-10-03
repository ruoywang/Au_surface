# size_complex_test_v1 — can 6×6 carry the rough environments, and which multi-layer morphologies need more?

Opened 2026-10-02 (user plan). One combined test, separate from the frozen list in `rough_sampling_v1/rough200`
(which is not rewritten): (a) 6×6 trial cuts of the existing rough centres, (b) 6×6 trial cuts of three
multi-layer morphology classes, (c) at most five size-paired CP-DFT comparisons (10 states, ≤ 8 new). All new
computations count against the rough 200; the 140 old states are untouched.

## What the test answers (and what it does not)

| question | test | not a failure criterion |
|---|---|---|
| does 6×6 keep the target atomic environment? | atom mapping, neighbourhood, declared morphology relation, seam checks | "the periphery is thinner than in 8×8" |
| does 6×6 give a complete, sensible CP-DFT label? | convergence, electronic state, energy/forces, 15 fields, charge closure | "SCF finished" alone |
| how much does the response change and how much is saved? | a few same-core, same-U 6×6 / 8×8 pairs, local comparison | total charge / energy / fields of different periodic systems being unequal |

8×8 is the larger-periodicity reference, not a proven converged truth. A geometrically legal, completely labelled
small cell is a training sample; size insensitivity is a separate conclusion.

## Status

- **First batch (2026-10-02 ~08:30):** reduced to two 8×8 references that pair with valid 6×6 cells of the same
  centre — `21015732` PA1_s51_f082500_a3617 (A, step/kink, −0.42 V, 281 Au; 6×6 = 163 Au) and `21015726`
  PB1_s51_f082500_a3929 (B, island edge, +0.46 V, 302 Au; 6×6 = 172 Au). The preferred C reference (21015735) has
  no valid 6×6 (exposed layer-1 atom and seam-created environments in every 6×6 frame), so B replaced it (both
  signs kept). Six jobs cancelled (21015725 after 38 min in its first SCF, 3 electronic steps; five pending).
  All other logical states of the frozen list are `held_for_resize` in `rough_sampling_v1/dft/queue.json`.
- Steps below in progress; new DFT only after the capped list and budget are confirmed once.

## Steps

1. Code generalised for more layers (`rough_sampling_v1`): exposed atoms searched in all non-fixed layers, not the
   top three; registration / occupancy fingerprint over every layer present; protected-atom sets and coded
   failure reasons in the cut; `md_driver`/`trajio` take a parents and MD directory.
2. `size_screen.py`: every frozen centre re-cut from its original MD frame with `CELLS = [(6,6),(8,8)]`,
   `R_CORE = 6.0` (unchanged), `target_environment` = centre + 6 Å neighbourhood, `protected_atom_ids` = that set;
   `preferred_cell = 6x6` when it passes, otherwise the coded reason and 8×8; no 10×8 escalation.
   → `size_screen.csv`.
3. `build_complex_parents.py`: six parents (two per class) on the same 32×32×4 base, every upper level from ONE
   seven-layer fcc build, each upper site supported by its three fcc hollows below (checked), four-layer base kept:
   M1 tiered islands (2-tier, 3-tier), M2 step bunching (two close steps merging locally into a double step),
   M3 open valley between two multi-layer ridges (narrow, wide; valley floor = complete original terrace).
   Same FLARE / MD protocol / 300 K frame rule; ~4 centres per parent → ~24 candidates for the 6×6 screen
   (with declared targets: upper edge + lower terrace + sidewall foot; both steps + middle terrace; both walls +
   floor). → `complex_gallery/`.
4. Vertical box: a three-level protrusion tops at ≈ 19.4 Å against the production limit SOL_Z1 − 15 = 19.6 Å;
   actual heights are checked and, if needed, BOTH members of a pair get the same extended vertical setting with a
   new configuration id (no compression, no deletion of the check).
5. Pairs (≤ 5): S1 step/kink and S2 island edge reuse the two kept 8×8 (6×6 targets their converged μ_e);
   M1–M3 one pair each (6×6 first, 8×8 matched to its actual state), |U| 0.35–0.5 V, both signs across the three;
   FERMICONVERGE 0.001 eV for pair runs, everything else at the production standard, k-mesh by the production rule.
   → `paired_dft_results.csv` (local RMS ΔF on the matched atom set C, max ΔF and where, PHI / S_ion / n⁻ / Γ⁻
   in the same region above the core, σ as a supplement; measured cost: atoms, k, FFT grid, bands, SCF steps,
   CP rounds, SCF time, field-writing time, node-hours, peak memory, output size).
6. `replacement_proposal.csv`: which unsubmitted large cells become 6×6, which quota goes to multi-layer
   environments, which large cells stay. No re-draw of potentials, no re-run of the 32 MD, no potential scans.

## Results so far (2026-10-02 ~10:00)

- **6×6 screen of the frozen 200** (`size_screen.frozen.csv`): 158 → 6×6 (112–180 Au, median 159; saves a median 119
  atoms per state), 33 → 8×8 (6×6 fails: protected set does not fit / seam contact unrepaired / seam-created
  environment / core neighbourhood changed), 9 → neither (all frozen 10×8). By class A 43/5/2, B 31/14/5, C 40/10/0,
  D 44/4/2. The two kept 8×8 references both have a valid 6×6 (163, 172 Au).
- **Multi-layer parents** after 200 ps: tiers survived (PM1b kept its third tier: 30 → 34 atoms; PM3b's level-6 cap grew
  64 → 114); tops at 19.72 Å (PM1b) and 20.23 Å (PM3b) in the DFT cell, above the 19.603 Å production limit → cells
  from them get the shared vertical-box extension in `pair_inputs.py`.
- **Complex screen** (`size_screen.complex.csv`, 24 declared-target centres): M2 step bunching 4 × 6×6 (180–190 Au: both
  steps + middle terrace kept), 1 × 8×8, 3 none; M3 valley floor with both walls → 8×8 only (320–330 Au; the protected
  set of ~70 atoms does not fit 17.6 Å), ridge-top edge → 6×6 (199 Au) on PM3a; **M1 as first built (one 0.24 ML
  tiered island, 245 lower-tier sites) cannot be cut at 6×6 or 8×8** — every frame produced seam overlaps down to
  1.9 Å because the island is larger than any cell; recorded as such, and two compact variants were built and run
  (PM1c: six 2-tier islands of 19/7 sites; PM1d: four 3-tier islands of 37/19/7) so the tiered-island environments are
  tested at a size a cell can contain.
- **Pairs** (`pair_plan.json`, `pair_budget.md`, inputs under `dft/`, nothing submitted): S1 (A step, −0.42 V, 6×6 163 Au
  vs the running 8×8 reference), S2 (B island edge, +0.46 V, 6×6 172 Au vs reference), M2 (middle terrace, +0.40 V,
  185 vs 314 Au), M3 (ridge-top edge, −0.40 V, 199 vs 345 Au); M1 to be added from PM1c/d. S1/S2 6×6 are provisional
  until the references' actual μ_e is known. Current estimate for the 6 new runs: 32 (seeded median) – 106 (cold max)
  node·h; with M1 ≈ 8 runs, ~45–160 node·h.
- **Compact tiered islands (PM1c/d, 2026-10-02 11:10):** PM1c's six small 2-tier islands (19/7 sites) flattened to
  one level during the 200 ps MD (recorded as found); PM1d's four 3-tier islands (37/19/7) lost their third tier and
  became 2-tier islands (148/76/28 → 156/96/0). Their upper-edge and sidewall-foot targets pass only at **8×8**
  (314–317 Au, seam-affected 0–4; 6×6 fails on the protected set and seam contacts for every origin). So M1 has no
  6×6 and no pair; M1 environments enter the list as 8×8 cells.
- **Final pair list: S1, S2, M2, M3 → 6 new single points** (FERMICONVERGE 0.001): S1 6×6 163 Au, S2 6×6 172 Au,
  M2 185 + 314 Au, M3 199 + 345 Au; estimate 32 (seeded median) – 106 (cold max) node·h; S1/S2 wait for the
  references' actual μ_e. **Nothing submitted; awaiting the one-time confirmation.**
- **First measured reference (A, 21015732, 8×8 281 Au, U = −0.42 V, finished 2026-10-02 16:10): COMPLETE** under the full
  acceptance — 15 fields parsed value by value, charge closure −1e-5 e, CONTCAR = POSCAR, μ_e = −4.4891 (2.0 meV
  from target). Two SCF rounds / one CP closure, 111 electronic steps, SCF 6.98 h, elapsed 7.33 h (0.35 h setup +
  field writing), MaxRSS 117.5 GB, 11 GB output. The NELECT start guess was 0.10 e from the final N_e (first-round
  μ_e 31 meV off). Cost model check: seeded-median prediction 412 min vs measured 440 min (×1.07); cold-max
  prediction 1376 min was 3.1× too high for this seeded run. Forces: mobile median 0.056, 90th pct 0.61, max
  0.88 eV/Å; the LARGE forces sit in the interior (thermal MD configuration kept as is: median 0.47 eV/Å over the
  43 interior mobile atoms), not in the seam band (FLARE-minimised: median 0.043, max 0.84) — the seam is not the
  anomaly. σ = −5.27 μC/cm² (N_e − 11N = +1.575 e), implying U_pzc ≈ +0.04 V with the dataset capacitance (±0.08).
  S1's 6×6 input is now definitive (TARGETMU = −4.4891, the reference's actual μ_e).
- **Second reference (B, 21015726, 8×8 302 Au, U = +0.46 V, finished 18:00): COMPLETE** — 15 fields parsed, closure
  −1.5e-5 e, CONTCAR = POSCAR, μ_e = −5.3746 (7.5 meV from target), N_e +0.10 e from the start guess. Two SCF
  rounds (the second oscillated at |ΔE| ~1e-6 for ~15 steps before converging), 132 steps, SCF 9.54 h, elapsed
  9.89 h, MaxRSS 122.9 GB, 11 GB. Forces: mobile median 0.105, max 1.53 eV/Å, again in the thermal interior
  (median 0.62) not the seam band (median 0.08). σ = +4.93 μC/cm² (N_e − 11N = −1.475 e).
- **Measured cost model** (`rough_sampling_v1/dft/cost_model.json`, two points): k(275 Au) = 471 (median) / 516
  (max) min, i.e. the dataset_v1 seeded extrapolation (399/563) was right and the cold envelope (1011/1346) 2–3×
  too high for seeded runs. Walltimes now k_max × (N/275)^1.5 × 2 (margin 2× until four points): 6×6 8–11 h,
  8×8 21–25 h. **Budgets with the measured model:** the 6 pair runs 37–41 node·h; the 198 remaining states under
  the proposal (141 6×6 + 15 multi-layer + 33 8×8 + 9 10×8; 112–381 Au, median 163) **932–1021 node·h** and
  ≈ 1.5 TB, against 1604–1757 node·h for the same states as frozen large cells (a 0.58× ratio based on measured
  large-cell costs, extrapolated to the resized set with N^1.5 — not yet a measured small-cell saving; the two
  numbers are the predictions of the two normalisation constants, not a confidence interval) and 757 node·h for
  all of dataset_v1.
- **Decision 2026-10-02 (user): the six pair single points are approved and SUBMITTED** (18:25; SLURM 21023280–85:
  M3 8×8 345 Au, M2 8×8 314, M3 6×6 199, M2 6×6 185 with `highmem,standard,wholenode`; S1 6×6 163 and S2 6×6 172
  with `standard,wholenode`; FERMICONVERGE 0.001; walltimes 8–25 h from the measured model). S1/S2 run at the
  references' actual μ_e, i.e. actual U = −0.4180 V and +0.4675 V (the assigned −0.42 / +0.46 remain the original
  draw; comparisons use the actual states). No M1 pair, no potential scan, no relaxation. The replacement framework is
  accepted: small cells first, 15 multi-layer states including the four pair states, the nine 10×8-only candidates
  kept but queued last and not submitted now; the remaining production waits for the pair results.
- **Unified compute list** (`merge_lists.py` → `unified_queue.csv`, 200 entries exactly): unique entries keyed by
  (final geometry, configuration, TARGETMU) with the old state ids as aliases. Three frozen size pairs whose
  members both resize to the same 6×6 collapsed into one entry each (PA2_s51 a4064, PB2_s51 a4040, PD2_s51 a4565);
  the two M2/M3 geometries that the proposal had given other potentials keep the TEST potential (+0.40 / −0.40 V;
  the proposal's +0.39 / −0.27 V are recorded, not computed); 15 multi-layer entries (M1 5 × 8×8 train, M2 2 pair
  + 3, M3 2 pair + 3; same-centre pairs share the test split); one resized single-layer entry dropped to respect
  the cap. Status: 2 complete, 6 submitted, 183 pending production (33 364 Au in total), 9 pending last (10×8).
- **Pair runs, observations while running (2026-10-03 05:30):** the 6×6 cells take 60–72 s per electronic step against
  195–300 s for the 8×8 cells on the same node type (~4× per step). CP behaviour differs with cell size: the CP
  loop's first correction is a fixed 0.10 e step; in the 8×8 cells 0.10 e moves μ_e by ≈ 26–28 meV, in the 6×6
  cells by ≈ 50 meV (smaller area, smaller cell capacitance), so S1 6×6 overshot from −4.4 meV to +45.8 meV after
  its second round and needs a third (secant) round even though its start guess was within 5 meV — the 0.001 eV
  pair tolerance will cost 6×6 runs an extra round or two. The multi-layer geometries start further from their
  target than the single-layer ones (first-round μ_e misses: M2 8×8 53 meV, M3 8×8 42 meV, M2 6×6 119 meV, S2
  6×6 46 meV vs 31–38 meV for the two references): their own PZC lies further from the dataset-median PZC used
  for the NELECT guess, which is itself a result about these morphologies.
- **Pair results, single-layer pairs S1 and S2 (2026-10-03 06:30; `paired_dft_results.md/.csv`, `profiles_S1.png`,
  `profiles_S2.png`).** Both 6×6 members are COMPLETE labels on their own (CP closed within 0.05 meV of the reference's
  μ_e, S1 3 rounds / 147 steps, S2 3 rounds / 152 steps; 15 fields value-parsed; charge closure; CONTCAR = POSCAR).
  Measured cost, 6×6 vs 8×8 of the same centre: S1 2.41 h vs 7.33 h elapsed, S2 2.63 h vs 9.89 h (0.33× and 0.27×;
  per electronic step 56–59 s vs 226–260 s); MaxRSS 50–52 GB vs 118–123 GB; output 5.3 GB vs 11.2 GB per state. The
  6×6 runs needed MORE electronic steps (147/152 vs 111/132) because of the 0.001 eV tolerance and the capacitance
  overshoot noted above, so the per-step saving is larger than the elapsed ratio.
  Matched atoms (same parent ids within R_CORE = 6 Å in both cells): S1 14, S2 17 — the core and the fixed layers under
  it. RMS ΔF over them 0.059 (S1) and 0.093 eV/Å (S2) against RMS |F| 0.21 / 0.35 eV/Å on the same atoms in the 8×8
  (ratio 0.28). The matching is by neighbour IDENTITY; the 6×6 seam lies 4.2 Å (S1) / 5.1 Å (S2) from the centre atom,
  so for most matched atoms a neighbour inside 6 Å is a seam-band atom that the repair moved or a periodic image with
  the thermal seam mismatch (relative shifts up to 0.30 / 0.57 Å). Split by that: atoms whose whole 6 Å environment is
  geometrically identical (< 0.02 Å) have RMS ΔF 0.048 eV/Å (S1, 1 atom = the centre, ΔF 0.084 on |F| 0.46) and 0.022
  eV/Å (S2, 6 atoms incl. the centre, ΔF 0.045 on |F| 0.21); atoms with a moved neighbour 0.060 / 0.114 eV/Å, the two
  largest 0.21 and 0.43 eV/Å, both on top-layer atoms 5.5–6.2 Å from the seam. Nothing in the seam band itself enters
  the comparison. Reading: the label change that is not explained by a changed geometry is ≈ 0.05–0.08 eV/Å on the
  centre atom (two cases, no error bar); atoms with the seam inside their 6 Å shell carry 0.1-eV/Å-scale differences.
  Electrolyte side (4 Å disc above the centre atom, heights ≥ 4 Å): |ΔPHI| ≤ 2.0 / 1.5 meV everywhere; the S_ion onset
  (accessible boundary) sits at 6.0 vs 6.0 Å (S1) and 7.0 vs 6.75 Å (S2), max |ΔS_ion| 0.09 / 0.15; n⁻ peak differs by
  7 % / 19 % of the 8×8 peak, Γ⁻ over the 12 Å window by −3.6 % / −5.0 %; the S1 n⁻ difference is what a 2 meV
  potential difference gives through the Boltzmann factor, the S2 one also contains the 0.25 Å cavity shift. σ (whole
  cell) −4.99 vs −5.27 and +4.75 vs +4.93 μC/cm² — different defect densities, supplement only.
  **M3 6×6 (ridge-top edge, −0.40 V, 199 Au) COMPLETE at 06:33** on its own: 3 rounds (first-round μ_e 4.8 meV off, the
  fixed 0.10 e step then overshot by 45 meV, secant round closed at −4.507122 vs −4.5071), 145 steps at 71 s/step, SCF
  2.86 h, elapsed 3.01 h, MaxRSS 54.4 GB, 5.0 GB. Its 8×8 partner is in round 3 (14.7 meV off after round 2, 292 s/step,
  11.6 h so far); M2 8×8 round 4 (13 meV off, 255 s/step); M2 6×6 round 3 (74 meV off after round 2, 73 s/step).
  **M2 6×6 (step-bunching middle terrace, +0.40 V, 185 Au) COMPLETE at 07:44** on its own: 4 rounds (first-round μ_e 119
  meV off — the largest start error of the six — then 74, 6.8 meV, closed at −5.306555 vs −5.3071 = 0.5 meV), 212
  steps at 68 s/step, SCF 4.00 h, elapsed 4.15 h, MaxRSS 55.4 GB, 5.0 GB. All four 6×6 members are therefore complete
  labels; the two 8×8 partners (M2, M3) are each in a fourth round at 255–292 s/step.
  **Cost model with 6 measured points** (`dft/cost_model.json`, 07:45): k(275 Au) median 372 / max 516 min; margin 1.5×
  (≥ 4 points). Per-step times are what the cell size sets: 59–75 s for the four 6×6 runs (163–199 Au) against 238–270 s
  for the two 8×8 (281–302 Au), a factor 3.5–4.5 for 1.5–1.8× the atoms; the number of steps is set by how many CP
  rounds the fixed-step/secant loop needs (145–212 for 6×6 at the 0.001 eV tolerance, 111–132 for the 8×8 at 0.01),
  which is why k spans 293–451 min for the 6×6 runs and 426–516 for the 8×8. Remaining 192 entries of the unified list
  (36 454 Au, median 163; 142 6×6, 39 8×8, 11 10×8): **682–922 node·h with the size-resolved k (6×6 entries at 293–451,
  larger at 426–516), 711–987 node·h with the pooled median/max**, ≈ 1.4 TB (11.2 GB per 290 Au, linear in atoms,
  measured at two sizes). The same frozen states as large cells (198 held, 55 383 Au): 1269–1757 node·h (pooled) and
  ≈ 2.1 TB. The ranges are the predictions of two normalisation constants, not confidence intervals. Production runs
  would use FERMICONVERGE 0.01 like the references, so their 6×6 round counts should sit at the low end.
- **Why the 8×8 pair runs need so many CP rounds (2026-10-03 08:45, from `CEP-HALF/src/main.F` lines 3208–3243):**
  with NESCHEME 3 the first correction is the fixed INIT_ECHANGE = 0.10 e; afterwards the loop uses the secant
  capacitance C = ΔN/Δμ, clipped to CAP_MAX = 2.0 e/eV. Measured C: 6×6 cells 2.0–2.2 e/eV (cap barely binding →
  near-full secant steps), 8×8 cells 3.7–3.9 e/eV (cap binding → every step is ≈ 52 % of the needed one, the error
  halves per round: M2 8×8 52.8 → 26.9 → 13.0 → 6.2 meV, M3 8×8 42 → 14.7 → 6.6 meV). At the production tolerance
  0.01 eV this cost the two 8×8 references nothing (2 rounds + closure); at the pair tolerance 0.001 eV an 8×8 run needs
  about three more rounds from 6 meV. Projection at 08:45: **M2 8×8 (job 21023281, ends 15:57) needs rounds 5–7 ≈ 90–125
  steps × 250 s + field writing = 6.6–9.0 h against 7.2 h left — may TIME OUT**; M3 8×8 (21023280, ends 19:57) needs ≈
  124 steps × 289 s + 0.4 h ≈ 10.4 h against 11.2 h left — should finish. A timed-out run leaves no fields (written at
  the end). Contingency prepared, NOT submitted: `dft/M2_8x8__mu-5.3071_restart/` = the same inputs with NELECT set to
  the secant estimate 3452.4289 (the measured slope 0.2600 eV/e is constant to 1 % over three intervals, so round 1
  should land within ~0.3 meV and the loop close at once: ≈ 55 steps ≈ 4 h on standard). For production nothing
  changes (0.01 eV tolerance, 6×6 cells unaffected); for any future 0.001-eV run on ≥ 8×8 cells CAP_MAX ≈ 5 would remove
  the halving without changing the converged state (CAP_MAX only shapes the path).
- **M3 pair complete (8×8 finished 14:23; `paired_dft_results.md`, `profiles_M3.png`).** M3 8×8 (345 Au, −0.40 V): 6 CP
  rounds at the 0.001 tolerance, closed 0.6 meV off, 253 steps, elapsed 19.42 h, MaxRSS 134.8 GB, 11.2 GB; at the 0.01
  standard it would have stopped after round 3 (148 steps, ≈ 11.5 h). Matched atoms 18: RMS ΔF 0.082 eV/Å against RMS
  |F| 0.30 (ratio 0.27); the centre atom (identical 6 Å shell) ΔF 0.039 on |F| 0.38; the 17 atoms with a seam-moved
  neighbour (shifts ≤ 0.31 Å) 0.084, largest 0.40 on a top-layer atom 5.4 Å from the seam. Electrolyte side: |ΔPHI| ≤
  2.0 meV, same accessible boundary (6.25 Å), max |ΔS_ion| 0.08, n⁻ peak differs 3 %, Γ⁻ identical to 3 digits. Same
  picture as S1/S2: no electronic cell-size effect above the few-meV / few-% level; the force differences follow the seam
  geometry. Cost: 6×6 3.01 h vs 8×8 19.42 h as run (0.001), 3.0 h vs ≈ 11.5 h at the 0.01 standard.
  **Cost model rebuilt on the 0.01 standard** (`recalibrate_cost.py`, 7 points): runs made at 0.001 are scaled to the
  electronic steps up to their first CP round within 0.01 (S1 52 of 147 steps, S2 122/152, M2 6×6 175/212, M3 6×6
  52/145, M3 8×8 148/253; the 0.01 references unchanged). k(275 Au): 6×6 runs 114–375 min, 8×8 runs 426–516, pooled
  median 375 / max 516. Remaining 192 entries: **503–846 node·h size-resolved, 717–986 pooled**, ≈ 1.4 TB; the 198 frozen
  large cells 1277–1757 node·h (pooled). Two normalisation constants, not an interval.
- **RULE (user, 2026-10-03 11:25): every constant-potential calculation uses FERMICONVERGE = 0.01 — production, tests,
  pairs, restarts, anything with LCEP. It is not to be changed again for any purpose.** The 0.001 used for the six pair
  runs came from the earlier Step-8×1/8×2 consistency test (2026-09-25, 36–72 Au, where the CAP_MAX clipping did not
  bite) and the pasted size-test plan's "建议", and was adopted without re-deriving its cost for 8×8 cells; that was
  my failure to review an inherited setting. Enforcement: `pair_inputs.write_member` refuses any other value,
  `rough_dft.submit` refuses to submit an INCAR that does not carry 0.01, and the prepared restart input
  `dft/M2_8x8__mu-5.3071_restart/INCAR` now carries 0.01. Where two cells must sit at the same electronic state, the
  second member targets the first member's ACTUAL converged μ_e (as S1/S2 did); the tolerance is not tightened.
  The two 8×8 runs still in progress carry 0.001 in their INCAR (cannot be changed mid-run); under the 0.01 rule they
  would already have closed (M2 at round 4, 6.2 meV; M3 at round 3, 6.6 meV).
- `replacement_proposal.csv`: 156 resize_to_6x6, 33 keep_8x8, 9 keep_10x8_no_smaller, 2 kept references; slots
  proposed for the admitted multi-layer candidates (M2 6×6 ×4, M2 8×8 ×1, M3 8×8 ×4 + 6×6 ×1, M1 8×8 ×5) — see the
  csv for the released single-layer states and the inherited U bins.
