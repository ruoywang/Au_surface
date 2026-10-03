# rough_sampling_v1 — old-state selection, rough parent surfaces, representative small cells

Stage opened 2026-10-01 (user plan, sections 一–九). Goal: ~140 selected old CP-DFT states + ~200 new CP-DFT
single points on small periodic cells cut from large, roughened Au(111) surfaces, so the training set gains
new local environments instead of more copies of the hand-built ones. Everything here lives inside this
folder; `dataset_v1/` and `05_production/` are read only.

**Status 2026-10-02 06:40 — list FROZEN after approval; first batch of 8 submitted; 192 pending (`CANDIDATE_REPORT.md`).**
MD 32/32 done; 1000 centres → 937 cut (1043 cells) → 200 states (40/5/5 per class, 32 parents, 2–10 cells per
parent, 177 × 8×8 + 23 × 10×8, 202–383 atoms, all k-mesh 2×2×1, 4 size pairs, seed 20261002). The review of
37a9299 required two fixes before freezing — an order-independent occupancy fingerprint registered on the frozen
bottom layer (verified on all 1043 cells and in the reviewer's 200 × 15-permutation test: 0 changes) and
value-level parsing of the PAW augmentation blocks (regression `test_fieldio.py`, 11/11) — after which the list
was re-selected: 162 states kept with their potentials, 38 replaced (bins inherited), class × split coverage
unchanged (`rough200/rough200_summary.md` lists the pairs). **First batch submitted 2026-10-02 06:3x: SLURM
21015725–21015730, 21015732, 21015735 on `wholenode` (two per class, one U > 0 and one U < 0, 3 × 10×8 + 5 × 8×8,
17–39 h walltime each). The remaining 192 are not submitted until the budget is re-approved on the measured cost
of these 8.**

## Pipeline and scripts

| step | script | output | status |
|---|---|---|---|
| 1. select old states | `select_old_states.py` | `old_selected_140.csv`, `old_excluded.csv`, `label_availability.jsonl`, `selection_report.md` | done: 140 states (80 base / 48 non-eq / 12 supplementary) over 52 geometries, every row with a reason |
| 2. parent surfaces | `build_parents.py` | `parents/*.extxyz|poscar`, `parents_manifest.json`, `parents_summary.md`, `_parents_overview.png` | done: 32 parents, 32×32×4 Au(111), 3840–4615 atoms (A 4571–4608, B 4250 / 4454, C 3994 / 3840, D 4604–4615), classes A strips / B islands / C pits / D transfer, min d 2.940 Å |
| 3. candidate-generation potential | `potential/PROVENANCE.md`, `check_potential.py` | `potential/behaviour_check.{md,json}` | done (below) |
| 4. LAMMPS + FLARE build | `env/BUILD.md` | `env/src/lammps-22Jul2025/build_mpi/lmp` | done, in-project, MPI |
| 5. MD on the parents | `md_driver.py`, `lmpio.py`, `trajio.py` | `md/<pid>/{in.lammps,data.lammps,slurm.sh,traj.lammpstrj}`, `md/md_manifest.json`, `md/md_summary.json` | done: 32/32 (SLURM 21014644–21014675, 45–47 min each on 16 ranks); sampled 300 K frames at 289–307 K; 110–160 atoms per parent changed layer |
| 6. centres | `select_centres.py` | `centres/centres.jsonl`, `centres_summary.md`, `parent_split.json`, `frames_used.json`, `rejects.json` | done: 1000 centres from 15 365 legal exposed atoms on 256 frames; CN strata 228/166/439/167; 8 close-contact rejections, 0 detached |
| 7. cut + repair | `extract_cells.py` (8 chunks + `--merge`) | `cells/<cell_id>/{POSCAR,cell.extxyz}`, `cells_manifest.jsonl`, `cells_summary.md` | done: 937 of 1000 centres pass (1043 cells; test centres in both sizes); 63 dropped, almost all because the seam would have created a (layer, CN) environment the parent does not have |
| 8. the 200 | `finalize_200.py` | `rough200/rough200_manifest.json` (**frozen = true**, 2026-10-02), `rough200_summary.md`, `rough200_cells.csv` | FROZEN after approval: 200 states; 274 occupancy duplicates and 9 seam-heavy cells excluded; fingerprint invariance verified on all 1043 cells |
| 9. DFT inputs, queue, state machine | `rough_dft.py`, `fieldio.py`, `test_fieldio.py` | `dft/<cell_id>/rough__mu<mu>/{POSCAR,INCAR,KPOINTS,POTCAR,job-run}`, `dft/queue.json`, `dft/first_batch.json`, `dft/budget.md` | 8 submitted, 192 pending, 38 withdrawn (replaced before the freeze); submit requires a frozen list AND `--confirm` |
| 10. figures + report | `render_rough_lineage.py`, `render_candidates_page.py`, `make_candidate_report.py` | `figures/rough_lineage_<class>.png`, `figures/candidates.html` (all 200, clickable), `CANDIDATE_REPORT.md` | done |

Python for this stage: `rough_sampling_v1/pyrun_rs.sh` (same interpreter as `scripts/pyrun.sh`, plus
`env/pylib` with dscribe 2.1.2 / numpy 2.0.2 / ase 3.26.0 installed with `pip --target`; nothing in `~/.local`).

## Potential: provenance, header fix, behaviour

Owen et al. 2024, Materials Cloud record 71861-xfe58 (DOI 10.24435/materialscloud:va-hx), `Au_training.zip`
md5-verified; model B2, n_max 8, l_max 4, cutoff 6.0 Å, Au only, NormalizedDotProduct kernel power 2 (from the
authors' `Au_master.yaml`). The 2022 export has `2` on its second line where the 2025 FLARE plugin expects
`<power> <kernel string>`, so a derived file `*.header2025.flare` differing in that one line is used
(`potential/PROVENANCE.md` records both checksums). The potential files are not redistributed.

Behaviour on this project's own structures (no DFT run; `potential/behaviour_check.md`): bulk minimum a0 = 4.160 Å
(project 4.158); the FLARE forces evaluated AT the CP-DFT-relaxed geometries are ≤ 0.06 (flat), 0.08 (step),
0.09/0.12 (fcc/hcp adatom), 0.18 (island-7), 0.08 eV/Å (pit-7) — the CP-DFT forces there are ~0 by construction,
so these numbers show how far the FLARE minimum sits from the CP-DFT one; they are NOT a measured FLARE-vs-CP-DFT
force error on general configurations. Minimising moves mobile atoms by RMS 0.02–0.05 Å (max 0.23 Å on the
island); signs: hcp adatom +0.047 eV above fcc, adatom 0.57 eV less bound than bulk, 8-atom step strip +0.61 eV,
7-atom island +1.89 eV, 7-atom pit +1.80 eV. Sanity signs for a candidate-generation tool; nothing here becomes a label.

## MD protocol (per parent, fixed before any run)

minimise → 300 K 20 ps → ramp to 600 K 20 ps → 600 K 100 ps → ramp to 300 K 20 ps → 300 K 40 ps; 2 fs; Nosé–Hoover
(0.1 ps) on the mobile atoms only; bottom two layers frozen (`fix setforce`, zero velocity); one frame per ps;
seed = 1000 + parent seed. Benchmark: 4580 atoms, 3.55 steps/s on 4 ranks, ≈ 33 steps/s on 16. If nothing
reconstructs or merges, that is recorded as found; no run is extended or heated further.

**Frame sampling rule (step 6):** only 300 K frames are sampled — 10, 20 ps (initial hold) and 165, 172, 179,
186, 193, 199 ps (final hold). The 600 K segment is the exploration that shapes the surface; its hot frames
are not sent to DFT, so thermal disorder in the new states stays at the 300 K level. Frames are located by
their LAMMPS TIMESTEP (`trajio.py`), checked against the protocol stage AND the logged temperature of the
mobile group (300 ± 50 K); a run without its DONE line, a dump that does not end at step 100 000, a missing
frame or an off-stage temperature stops the selection and names the parent. There is no fallback to the
initial structure or to an even spread over a short trajectory; `--test-initial` is an explicit test mode.

**Legality of a centre** is decided by geometry — exposed, connected to the substrate through the neighbour
graph, no contact < 2.5 Å — not by coordination number. CN is recorded and used for stratified sampling
(cn ≤ 5 adatom/kink, 6–7 edge, 8–9 terrace, ≥ 10), quota ∝ √count per stratum, so low-coordinated environments
are sampled deliberately. A centre with no SOAP neighbour at cosine > 0.98 in the pool is flagged rare, not deleted.

## Cutting and repair (step 7), what is and is not claimed

Core = centre atom + everything within 6 Å (protection radius, not a physical cutoff). The core is carried over
to 1e-3 Å after one translation and its 6 Å neighbourhood is the same atom set; the cut origin is chosen among
lattice-compatible origins with boundaries midway between atom rows, by the smallest seam mismatch. Repair order
is origin → size (8×8 then 10×8) → drop. A cell passes only if: no contact < 2.5 Å (only the self term is
excluded from the distance list, so two atoms on one site fail), every atom ≥ 3 neighbours AND every atom
reachable from the fixed bottom layers through the neighbour graph (no detached cluster), bottom layer
complete, no exposed atom deeper than one layer, no seam artefact — every atom's CN in the cell is compared
with ITS OWN CN in the parent frame; atoms that lost ≥ 2 neighbours were cut by the seam (unavoidable where a
strip or island larger than the cell is terminated; their number is recorded per cell as `n_seam_affected`),
and the cell is rejected when the seam CREATES a (layer, CN) environment that fewer than 3 atoms of the parent
surface have — e.g. a cn-4 fragment on a surface whose lowest top-layer CN is 6. A natural adatom (cn 3 in parent
and cell) is not an artefact. The earlier "periphery vs core minimum" rule could not tell the two apart and was
replaced. Finally the per-atom periodic mapping check: a third of the atoms moved to periodic equivalents, all
re-wrapped, order permuted — the centre's neighbour set, the CN list and the minimum distance must not change.
Contacts < 2.5 Å are classified: a pair at exactly its parent distance is a thermal pair inherited from the
frame (accepted above a 2.3 Å floor and counted), a pair created by the seam is a problem. Only seam-created
contacts are sent to the minimiser (strain can relax; occupancy cannot), and the minimiser moves only the seam
band — atoms within 2 rows of either seam in the top levels, never the core, the bottom layers or the rest of
the periphery, which keep their thermal state; the number of atoms moved is recorded in the repair note. A
failing occupancy is dropped at once. The production cut is run as 8 frame-grouped chunks
(`extract_cells.py --chunk i/8`, then `--merge`); every frame is read and its neighbour graph built once. The
label later attached to a cell belongs to that reconstructed periodic electrode, not to the parent, and its
energy is not a "core energy".

## The 200 (step 8) and their potentials — a CANDIDATE list until approved

160 train / 20 val / 20 test, 50 per class; parents split 6/1/1 by seed inside each class and every frame or cut
of a parent stays in its split. Selection is two-layered: (1) SOAP of the centre atom in the reconstructed
cell (new local environments) and (2) the whole cell's upper-surface morphology — exposed-atom CN histogram,
relief, levels, adatom-level coverage, missing-terrace and exposed-pit fractions, island and pit counts/sizes,
step density — standardised and combined in the farthest-point distance, per class × split × CN stratum, with
caps per parent (≤ 1.5 × fair share of the bucket) and per (parent, frame) (≤ 2). Cells whose layer-3/4
lattice-occupancy patterns coincide within 2 sites under translation are one geometry; one is kept. 4 size
pairs (same core in 8×8 and 10×8, same U) inside the test split. One U per geometry, 10 bins × 0.1 V over
[−0.5, +0.5] V, 20 per bin, drawn after the geometry is chosen with a recorded seed; TARGETMU = −4.9071 − U is
the request, the state is labelled by its converged μ_e. The manifest is written with `frozen = false`;
`finalize_200.py freeze` (after approval) flips the flag without re-drawing, and `rough_dft.py submit` refuses
while it is false.

## DFT budget (step 9, nothing submitted)

Production standard K.8 imported from `scripts/production.py` (same INCAR, k-mesh rule, solvent window, job
script). Added: a NELECT start guess from the dataset-wide C = 11.59 μF/cm² and U_pzc = +0.006 V (recorded in the
INCAR comment; the CP loop decides N_e). Measured basis: the >200-atom single points of dataset_v1 took 1011 min
median / 1346 max cold, 399 / 563 NELECT-seeded (4 states). For 200 cells of 256–370 atoms the estimate is
roughly **1.4–2.2 k node·h if the seeding helps as measured, 3.5–4.6 k node·h cold** (dataset_v1 in total: 757), and
≈ 2.2 TB of fields (11 GB per state; scratch at 1.9 of 100 TB). Exact numbers per size are written to
`dft/budget.md` by `rough_dft.py prepare`; the walltime is computed per cell (cold envelope × 1.2: 256 atoms →
25 h, 300 → 31 h, 370 → 43 h) and recalibrated from the first batch. This is a new computing commitment of
roughly 2–6 × dataset_v1, to be approved on the actual candidate list, not assumed.

**Acceptance of a finished state** (`rough_dft.py status`): CP closed on target (|μ_e − TARGETMU| ≤ 0.011 eV),
VASP reached its timing summary, all 15 field files (CHGCAR, LOCPOT, PHI, PHISOLV, VSOLV, RHOB, RHOION, ELOC,
P, SVDW, SION, SSOLV, SCAV, SDIEL, POT) parsed value by value — grid header, ≥ nx·ny·nz finite floats, N_ions
augmentation blocks for CHGCAR/POT (`fieldio.py`), CHGCAR integrating to the CP electron count within 0.02 e,
CONTCAR identical to POSCAR. Anything less is `partial`, never `complete`.

**First batch** (`rough_dft.py submit --first 8`): two states per class, one at U > 0 and one at U < 0, cell
sizes alternating so both 8×8 and 10×8 appear, the largest |U| within each slot (longest CP walk). Printed
before anything is sent; sent only with `--confirm` AND a frozen manifest.

## Not done / gated

- **The remaining 192 states**: not submitted until the budget is re-approved on the measured cost of the first
  8 (actual μ_e / N_e, SCF and CP convergence, field acceptance, force range incl. seam band vs core, memory,
  SCF and field-writing time → recalibrated cost model and walltimes). Follow the 8 with `rough_dft.py status`.
- Known small imperfection: layer assignment is by z rounding, so an atom caught mid-hop between levels can be
  counted in the lower level (3 of the 200 cells show `missing_terrace_ML` = −0.05, i.e. three such atoms).
  It affects the morphology descriptors and the dedup fingerprint marginally, not the geometries themselves.
- Open points carried from dataset_v1 (unchanged): egg-box force error (audit K.14) and the trainer's energy target.

## Decision log

- **2026-10-03 — RULE (user): FERMICONVERGE = 0.01 for every constant-potential calculation — production, tests, size
  pairs, restarts, anything with LCEP — never to be changed again.** Background: the six size-pair runs of 2026-10-02 were
  given 0.001 (inherited from the Step-8×1/8×2 consistency test and a pasted suggestion, adopted without review); with
  CAP_MAX = 2 e/eV the CP loop of an 8×8 cell (C ≈ 3.8 e/eV) then only halves its μ_e error per round, needing 6–7 rounds
  and risking the walltime. `rough_dft.submit` and `size_complex_test_v1/pair_inputs.write_member` now refuse any other
  value. Two cells meant to share an electronic state are matched through the ACTUAL converged μ_e of the first member.
- 2026-10-02 — NELECT start guess kept for all rough states; no cold-start subset (user). The per-cell PZC of a
  rough state is therefore not read directly from a neutral first CP round; it is estimated from the single
  (U, σ) point with the dataset-wide capacitance (≈ ±0.08 V at |U| = 0.5 V, smaller near U = 0), and reported as such.
