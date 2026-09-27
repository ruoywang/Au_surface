# Complete Au surface dataset — computation plan v1 (2026-09-27)

Generated from `dataset_plan_v1.csv` by `scripts/build_dataset_plan.py`; edit the CSV, not this file. Scope is fixed here once; execution is batched by cluster limits afterwards. No DFT is submitted by this document.

## 1. Scope

**Included families** (each with main sampling and a few references, all in one table): flat Au(111); point defects / minimal clusters; straight strip steps; crystallographic vicinal step faces; kinks / edge rearrangement; single-layer islands; single-layer pits; reconstruction-related stacking environments; two composite morphologies as references.

**Explicitly excluded**: arbitrary grain boundaries, nanoparticles, multilayer spikes/towers, other low-index faces (Au(100) phase-2 stays recorded in defect_plan.md, not budgeted here), the full surface phase diagram, explicit water, any Cl species, MD trajectories. The 22×√3 herringbone cell itself (~63 Å period) is not built; R2 is its stripe/soliton approximation.

**Model, potentials, label standard**: fixed 4-layer slabs (bottom two layers fixed), a = 4.158 Å, single-sided 1 M implicit electrolyte, static CP-DFT single points; production configuration (PREC=Normal, ALGO=Fast, NPAR=16 hybrid, FERMICONVERGE=0.01) is the label standard; three common internal electron chemical potentials μ_e = −5.1071, −4.9071, −4.7071 eV (internal μ₀ reference, not V vs RHE); states are labelled by the actual converged μ_e/N_e. PREC=Accurate pilot results stay a separate reference group and are **not** counted as reusable coverage.

## 2. Sampling template (what one row means)

| column | definition |
|---|---|
| n_ref | reference configurations: the ideal bulk-truncated geometry, plus the CP-relaxed geometry when n_relax = 1 |
| n_relax | CP relaxations at μ_e = −4.9071 (movable atoms free, bottom two layers fixed), charged as 6 single points each; produce the relaxed reference |
| n_perturb | random-displacement configurations of the relaxed geometry on the movable atoms: Gaussian σ = 0.05 Å (seeds 1,2) and 0.10 Å (seeds 3,4); halved sets use seeds 1,3 |
| n_collective | collective deformations of the relaxed geometry: top interlayer spacing −3 %; in-plane strain +1 % (cell scaled); for strips the second one is an edge-row bend of ±0.15 Å along the period |
| n_path | local rearrangement images: linear interpolation between the relaxed initial and final geometries of the named move, end points excluded (they are their own rows or the relaxed reference) |
| mu_list | electron chemical potentials computed for every configuration of the row |
| n_states_target | n_configs × n_mu — the number of labelled states the row contributes |
| n_reusable_states | states already computed in the production standard with identical cell/window/k-mesh |
| cost | highmem node-hours from the measured production timings (see cost model in the script header) |

Perturbed and path configurations are single points on non-equilibrium geometries; none of them is relaxed back, so no near-duplicate minima are generated.

## 3. Master table

| ID | family | tier | coordinates | N | movable | ref | relax | pert | coll | path | μ_e (eV) | states | reusable | new DFT | cost (node·h) | purpose |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T-4x4 | flat Au(111) | main | library T.poscar (4x4x4 layers) | 64 | top 2 layers | 2 | 1 | 4 | 2 | 0 | -5.1071, -4.9071, -4.7071 | 24 | 0 | 24 | 17.6 | flat baseline; charging response of the ideal terrace; non-equilibrium flat configurations |
| Flat-16x1 | flat Au(111) | reference | library Flat-16x1.poscar (built 2026-09-27) | 64 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 2 | 1 | 0.6 | same-cell flat baseline for the Step-16x1 family; only the -4.7071 state is new |
| Flat-8x2 | flat Au(111) | reference | to build: Step-8x2 with its 8 strip atoms removed (same recipe as Flat-16x1) | 64 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 0 | 3 | 1.8 | same-cell flat baseline for the Step-8x2 main strip |
| V1 | point defect | main | library V1.poscar | 63 | top 2 layers | 2 | 1 | 4 | 2 | 0 | -5.1071, -4.9071, -4.7071 | 24 | 0 | 24 | 17.3 | single vacancy; charging response and forces around a missing top-layer atom |
| A1-fcc | point defect | main | library A1_fcc.poscar | 65 | top 2 layers + adatom | 2 | 1 | 4 | 2 | 1 | -5.1071, -4.9071, -4.7071 | 27 | 0 | 27 | 19.7 | fcc adatom; path image = adatom at the bridge site between fcc and hcp (the hcp end is A1-hcp) |
| A1-hcp | point defect | main | library A1_hcp.poscar | 65 | top 2 layers + adatom | 2 | 1 | 4 | 2 | 0 | -5.1071, -4.9071, -4.7071 | 24 | 0 | 24 | 17.9 | hcp adatom; not a reconstruction model; end point of the A1 fcc->hcp path |
| V2 | point defect | reference | library V2.poscar | 62 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 0 | 3 | 1.7 | divacancy; low weight |
| V3 | point defect | reference | library V3.poscar | 61 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 0 | 3 | 1.7 | triangular trivacancy; low weight |
| A3 | point defect | reference | library A3.poscar | 67 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 0 | 3 | 1.8 | three-atom cluster; low weight |
| Step-8x2 | strip step | main | library Step-8x2.poscar | 72 | top 2 base layers + strip (along-edge displacements independent) | 2 | 1 | 4 | 2 | 3 | -5.1071, -4.9071, -4.7071 | 33 | 0 | 33 | 25.7 | straight A/B strip edges 10 A terraces; collective = edge-row bending; path = one upper-edge atom detaching to the lower-terrace foot (3 images) |
| Step-16x2 | strip step | main | library Step-16x2.poscar | 144 | top 2 base layers + strip | 2 | 1 | 2 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 12 | 0 | 12 | 33.6 | 20 A terraces with along-edge-independent displacements; reduced sampling because of size |
| Step-8x1 | strip step | reference | library Step-8x1.poscar | 36 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 1 | 2 | 0.7 | ideal minimal-period strip; width series member |
| Step-16x1 | strip step | reference | library Step-16x1.poscar | 72 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 2 | 1 | 0.7 | ideal minimal-period strip; width series member |
| Step-24x1 | strip step | reference | library Step-24x1.poscar | 108 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 0 | 3 | 3.6 | ideal strip 30 A terraces; width extrapolation |
| Au211 | vicinal step face | main | library Au211.poscar | 48 | top 2 layers (all step-face atoms) | 2 | 1 | 4 | 2 | 0 | -5.1071, -4.9071, -4.7071 | 24 | 0 | 24 | 13.2 | A-type {100}-microfacet step face with 3-row terraces |
| Au221 | vicinal step face | main | library Au221.poscar | 72 | top 2 layers | 2 | 1 | 4 | 2 | 0 | -5.1071, -4.9071, -4.7071 | 24 | 0 | 24 | 19.8 | B-type {111}-microfacet step face representative |
| Au332 | vicinal step face | reference | library Au332.poscar | 48 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 0 | 3 | 1.3 | wider {111}-type terrace; same-family width check |
| Au554 | vicinal step face | reference | library Au554.poscar | 40 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 0 | 3 | 1.1 | widest {111}-type terrace in the library; same-family width check |
| Kink-A | kink / edge rearrangement | main | to build: Step-8x3 strip (base 8x3x4 = 96 + strip 12) with one extra edge atom on edge1 in one of the 3 rows | 109 | top 2 base layers + strip | 2 | 1 | 3 | 0 | 2 | -5.1071, -4.9071, -4.7071 | 21 | 0 | 21 | 33.2 | kink on the edge1 row; path = kink atom moving one site along the edge (2 intermediate images); needs ny=3 period - not reducible |
| Kink-B | kink / edge rearrangement | main | to build: as Kink-A with the extra atom on edge2 | 109 | top 2 base layers + strip | 2 | 1 | 3 | 0 | 2 | -5.1071, -4.9071, -4.7071 | 21 | 0 | 21 | 33.2 | kink on the other (inequivalent) edge; same path recipe |
| Step-8x2+foot-adatom | kink / edge rearrangement | reference | to build: Step-8x2 + one Au adatom at the lower-terrace foot site next to edge2 | 73 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 3 | 0 | 3 | 2.0 | edge attachment environment (attached-atom end point of the Step-8x2 detachment path) |
| Island-7-compact | single-layer island | main | library Island-7-6x6.poscar | 151 | top base layer + island | 2 | 1 | 2 | 0 | 2 | -5.1071, -4.9071, -4.7071 | 18 | 0 | 18 | 48.1 | compact 7-atom island in the 6x6 cell; path = one island-edge atom detaching to the terrace (2 images) |
| Island-7-elongated | single-layer island | main | to build: 7 atoms as a 2-row zigzag chain (3+4) on the 6x6 slab | 151 | top base layer + island | 2 | 1 | 2 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 12 | 0 | 12 | 36.1 | shape variant at equal atom count; edge-length/coordination contrast to the compact island |
| Island-19-8x8 | single-layer island | reference | library Island-19-8x8.poscar | 275 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071 | 2 | 0 | 2 | 9.9 | size extrapolation; two states only |
| Island-7-8x8 | single-layer island | reference | library Island-7-8x8.poscar | 263 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -4.9071 | 1 | 0 | 1 | 4.6 | cell-size check for Island-7 (6x6 vs 8x8) and same-cell partner of Island-19-8x8; one state |
| Pit-7-compact | single-layer pit | main | library V7.poscar / Pit-7 | 137 | top base layer + pit rim | 2 | 1 | 2 | 0 | 2 | -5.1071, -4.9071, -4.7071 | 18 | 0 | 18 | 41.6 | compact 7-vacancy pit; path = one rim atom moving into the pit (2 images) |
| Pit-7-trench | single-layer pit | main | to build: 7 vacancies as a 2-row zigzag trench in the 6x6 top layer | 137 | top base layer + pit rim | 2 | 1 | 2 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 12 | 0 | 12 | 31.2 | shape variant at equal vacancy count |
| Pit-19-8x8 | single-layer pit | reference | library Pit-19-8x8.poscar | 237 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071 | 2 | 0 | 2 | 7.9 | size extrapolation; two states only |
| Pit-7-8x8 | single-layer pit | reference | library Pit-7-8x8.poscar | 249 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -4.9071 | 1 | 0 | 1 | 4.2 | cell-size check for Pit-7 and same-cell partner of Pit-19-8x8; one state |
| R1-hcp-terminated | reconstruction-related | main | to build: T-4x4 with the whole top layer shifted to hcp registry (stacking fault at the surface) | 64 | top 2 layers | 2 | 1 | 2 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 12 | 0 | 12 | 10.6 | hcp-stacked surface region as the local model of a reconstruction domain (not an adatom) |
| R2-stripe-wall | reconstruction-related | main | to build: 16x1 cell, 4 base layers (64) + top layer of 17 atoms (one extra row, 6% compression) with fcc->hcp->fcc registry and two domain walls; validated by min distance >= 2.6 A, registry check, periodic closure | 81 | top layer + first base layer | 2 | 1 | 2 | 0 | 0 | -5.1071, -4.9071, -4.7071 | 12 | 0 | 12 | 14.2 | stripe-reconstruction approximation (herringbone soliton walls); ny=1 keeps the walls straight - valid only without along-wall perturbations |
| C1-island-near-step | composite | reference | to build: Step-8x3 strip (108) + compact 7-atom island on the lower terrace 1 row from edge2 | 115 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071 | 2 | 0 | 2 | 2.7 | island-step proximity environment; two states |
| C2-island+pit | composite | reference | to build: 6x6 slab with the compact 7-atom island and the compact 7-vacancy pit 2 rows apart (144 atoms net) | 144 | none (ideal only) | 1 | 0 | 0 | 0 | 0 | -5.1071, -4.9071 | 2 | 0 | 2 | 3.7 | island-pit pair at fixed separation; two states |

## 4. Totals and budget

| | structures | configs | target states | reusable | new DFT states | relaxations | cost (highmem node·h) |
|---|---|---|---|---|---|---|---|
| main | 16 | 106 | 318 | 0 | 318 | 16 | 413 |
| reference | 17 | 17 | 43 | 5 | 38 | 0 | 50 |
| **all** | 33 | 123 | 361 | 5 | 356 | 16 | 463 |

Per family (full scenario):

| family | structures | new DFT states | cost (node·h) | share |
|---|---|---|---|---|
| flat Au(111) | 3 | 28 | 20 | 4 % |
| point defect | 6 | 84 | 60 | 13 % |
| strip step | 5 | 51 | 64 | 14 % |
| vicinal step face | 4 | 54 | 35 | 8 % |
| kink / edge rearrangement | 3 | 45 | 68 | 15 % |
| single-layer island | 4 | 33 | 99 | 21 % |
| single-layer pit | 4 | 33 | 85 | 18 % |
| reconstruction-related | 2 | 24 | 25 | 5 % |
| composite | 2 | 4 | 6 | 1 % |

**Lean scenario** (perturbations halved on main rows; 2 potentials instead of 3 for perturbation/path configurations of rows with N > 100; references unchanged): 266 new DFT states, 367 node·h (79 % of the full scenario). Rows with the largest savings: Island-7-compact (12 h), Kink-A (11 h), Kink-B (11 h), Pit-7-compact (10 h), Island-7-elongated (8 h), Step-16x2 (7 h).

**Throughput**: highmem allows 2 running jobs per user; at 1 job = 1 node the full scenario is ≈ 10 days of continuous submission (8 days lean). Highmem is billed at 4× (audit §F); node·h × 128 cores × 4 ≈ 237 k SU full, 188 k SU lean. Warm-starting each perturbed configuration from its parent CHGCAR (ICHARG=1) would cut SCF steps by roughly a third; not assumed in the numbers above.

**Cost-model limits**: measured points cover 36–72 atoms only; the N^1.5 growth used for 108–275 atoms is an assumption, so the island/pit/kink rows (≈ half of the budget) carry the largest uncertainty. The first completed run of each size class recalibrates the CSV.

## 5. Structures that still have to be built (concrete specifications)

All builds use the existing generators' conventions (a = 4.158 Å, 4 base layers at z = 5.0/7.4/9.8/12.2 Å, Lz = 44.603 Å, SOL_Z0 = 9.801, SOL_Z1 = 34.603 Å, DIPOL 0.5 0.5 0.5, k-mesh by the N·|a| ≈ 35 Å rule) and pass the existing checks (minimum Au–Au distance ≥ 2.6 Å after perturbation, ≥ 2.9 Å for ideal builds; periodic-image isolation; registry check of every added layer).

| ID | build recipe | validation |
|---|---|---|
| Flat-8x2 | remove the 8 strip atoms from Step-8x2 (as done for Flat-16x1) | base atoms unchanged (bijective match), 64 atoms |
| Kink-A / Kink-B | Step-8x3 (base 96, strip 12) plus one extra edge atom on edge1 (A) or edge2 (B) in one of the three rows; the ny = 3 period is the minimum for one kink per period and is **not** reducible | added atom in a hollow of the layer below (registry), CN of kink atom = 5, min distance 2.94 Å |
| Step-8x2+foot-adatom | Step-8x2 plus one adatom in the fcc hollow of the lower terrace adjacent to edge2 | registry, min distance |
| Island-7-elongated | 7 atoms as a 3+4 two-row zigzag chain in fcc hollows on the 6x6 slab (same cell as Island-7-6x6) | connected island (every atom ≥ 2 island neighbours), image separation ≥ 8 Å |
| Pit-7-trench | 7 vacancies forming a 3+4 two-row zigzag trench in the 6x6 top layer | connected vacancy footprint, image separation as Pit-7 |
| R1-hcp-terminated | T-4x4 with the top layer translated by (a/√3)·[1̄10]-type shift into hcp registry | every top atom above a second-layer atom (hcp), min distance 2.94 Å |
| R2-stripe-wall | 16x1 cell; top layer with 17 atoms on the 16-site row: two domains (fcc, hcp) joined by two walls of ~4 atoms each along a1, compression 1/16 spread over the walls; positions from a 1-D soliton profile, then constrained relaxation | min distance ≥ 2.6 Å before relaxation, registry fcc/hcp in the domain centres, periodic closure of the displacement field (total shift = one row spacing) |
| C1-island-near-step | Step-8x3 strip (108) + compact 7-atom island on the lower terrace, nearest island atom one row from edge2 | island–edge distance recorded; no island–image contact |
| C2-island+pit | 6x6 slab, compact 7-atom island and compact 7-vacancy pit with rims two rows apart | 144 atoms net; separation recorded |
| path images | A1: bridge site (1); Step-8x2 detachment: 3 images between the relaxed strip and the relaxed foot-adatom structure; kinks: 2 images of the kink atom moving one site along the edge; Island-7 / Pit-7: 2 images of one rim atom leaving/entering | interpolated on movable atoms only; end points are their own rows |

## 6. Execution batches (resource management only — the scope above does not change)

1. **Batch A (≤ 72 atoms, anchors)**: T-4x4, V1, A1-fcc, A1-hcp, Flat-16x1 (+1), Flat-8x2, Step-8x1 (+2), Step-16x1 (+1), Step-8x2, Au211, Au221, R1, V2/V3/A3, Au332/Au554, Step-8x2+foot-adatom. Relaxations first, then perturbations/paths.
2. **Batch B (73–151 atoms)**: R2, Step-24x1, Kink-A/B, C1, Pit-7-compact, Pit-7-trench, Step-16x2, C2, Island-7-compact, Island-7-elongated. The first completed run of each size recalibrates the cost model.
3. **Batch C (> 200 atoms, references)**: Pit-19-8x8, Pit-7-8x8, Island-7-8x8, Island-19-8x8 — one or two states each; submitted last, cancelled first if the budget tightens.

Batch order does not decide scope; if the budget is cut, apply the lean scenario (or a further global reduction of perturbation counts) across all rows rather than dropping a family.

## 7. Not decided or not required here

- The trainer's energy target and input keys (settled against the training code on its own machine, not here); this plan only guarantees that every state carries TOTEN, E_without_entropy, E(σ→0), GCE, forces, actual N_e/μ_e and the indexed 3D fields.
- Model interface tests or small-sample training are not prerequisites for starting Batch A.
- The closed step-mechanism study is not reopened; Step-16x1/8x1 results enter only as reusable states.
