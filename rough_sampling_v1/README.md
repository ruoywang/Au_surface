# rough_sampling_v1 — old-state selection, rough parent surfaces, representative small cells

Stage opened 2026-10-01 (user plan, sections 一–九). Goal: ~140 selected old CP-DFT states + ~200 new CP-DFT
single points on small periodic cells cut from large, roughened Au(111) surfaces, so the training set gains
new local environments instead of more copies of the hand-built ones. Everything here lives inside this
folder; `dataset_v1/` and `05_production/` are read only.

**Status 2026-10-02 (see the bottom for what is NOT done).** No rough DFT has been submitted.

## Pipeline and scripts

| step | script | output | status |
|---|---|---|---|
| 1. select old states | `select_old_states.py` | `old_selected_140.csv`, `old_excluded.csv`, `label_availability.jsonl`, `selection_report.md` | done: 140 states (80 base / 48 non-eq / 12 supplementary) over 52 geometries, every row with a reason |
| 2. parent surfaces | `build_parents.py` | `parents/*.extxyz|poscar`, `parents_manifest.json`, `parents_summary.md`, `_parents_overview.png` | done: 32 parents, 32×32×4 Au(111), 3840–4615 atoms (A 4571–4608, B 4250 / 4454, C 3994 / 3840, D 4604–4615), classes A strips / B islands / C pits / D transfer, min d 2.940 Å |
| 3. candidate-generation potential | `potential/PROVENANCE.md`, `check_potential.py` | `potential/behaviour_check.{md,json}` | done (below) |
| 4. LAMMPS + FLARE build | `env/BUILD.md` | `env/src/lammps-22Jul2025/build_mpi/lmp` | done, in-project, MPI |
| 5. MD on the parents | `md_driver.py`, `lmpio.py` | `md/<pid>/{in.lammps,data.lammps,slurm.sh,traj.lammpstrj}`, `md/md_manifest.json` | **32 SLURM jobs running** (21014644–21014675, `shared`, 16 ranks, 8 h cap; measured ≈ 33 steps/s → ≈ 50 min each) |
| 6. centres | `select_centres.py` | `centres/centres.jsonl`, `centres_summary.md`, `parent_split.json`, `frames_used.json` | tested on the initial parents (`centres_test/`); waits for the MD |
| 7. cut + repair | `extract_cells.py` | `cells/<cell_id>/{POSCAR,cell.extxyz}`, `cells_manifest.jsonl`, `cells_summary.md` | tested (`cells_test/`: 14 of 16 centres pass; the two drops were seam-isolated atoms) |
| 8. fix the 200 | `finalize_200.py` | `rough200/rough200_manifest.json`, `rough200_summary.md` | tested; potentials drawn once, refuses to overwrite |
| 9. DFT inputs, queue, state machine | `rough_dft.py` | `dft/<cell_id>/rough__mu<mu>/{POSCAR,INCAR,KPOINTS,POTCAR,job-run}`, `dft/queue.json`, `dft/budget.md` | tested dry (`dft_test/`); **submit is gated by `--confirm`** |
| 10. figure | `render_rough_lineage.py` | `figures/rough_lineage_<class>.png` | tested (`figures_test/`) |

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
complete, no exposed atom deeper than one layer, no periphery atom less coordinated than the core minimum − 1
(this is what rejected the two test drops: a seam had isolated a single atom), and the per-atom periodic
mapping check: a third of the atoms moved to periodic equivalents, all re-wrapped, order permuted — the
centre's neighbour set, the CN list and the minimum distance must not change. The label later attached to a
cell belongs to that reconstructed periodic electrode, not to the parent, and its energy is not a "core energy".

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

- Steps 6–7 run on the real MD frames by `run_after_md.sh` (stops if any run is incomplete; logs to
  `logs/pipeline.log`); steps 8–9 are run by hand and produce a CANDIDATE list + inputs + budget.
- **Rough DFT submission: gated twice** — the manifest must be frozen by approval and `--confirm` given.
- Open points carried from dataset_v1 (unchanged): egg-box force error (audit K.14) and the trainer's energy target.
