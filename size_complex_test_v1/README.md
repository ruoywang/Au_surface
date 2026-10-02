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
- `replacement_proposal.csv`: 156 resize_to_6x6, 33 keep_8x8, 9 keep_10x8_no_smaller, 2 kept references; 10 slots
  proposed for the admitted multi-layer candidates (to be revised with M1c/d).
