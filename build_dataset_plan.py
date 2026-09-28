#!/usr/bin/env python3
"""Generate dataset_plan_v1.md from dataset_plan_v1.csv (the single source of truth for the
complete Au dataset computation plan). Every multiplier is a CSV column; totals and costs are
computed here, never hand-summed.  Rev 2 (2026-09-27, after the user's decisions).

Sampling rule (frozen):
  reference configs (n_ref: ideal + relaxed)        -> all potentials in mu_list (3 for main rows)
  other configs (perturb + collective + path)        -> 3 potentials if N <= 100, else the 2 end points (-5.1071, -4.7071)
  n_configs        = n_ref + n_perturb + n_collective + n_path
  n_states_target  = n_ref*n_mu_ref + (n_perturb+n_collective+n_path)*n_mu_var
  n_new_dft        = n_states_target - n_reusable_states     (reuse = same production standard only)
Cost model c(N) [min per static CP single point, 1 highmem node = 16 MPI x 8 OpenMP, production config]:
  c(N) = 0.55 * N * max(1, N/72)**0.5   fitted to measured 36 -> 17.9 min, 64 -> 27-39 min, 72 -> 40-43 min;
  growth beyond 72 atoms is an assumption (untested above 72 atoms in this configuration).
Relaxation cost = f_relax * c(N) with f_relax in {4.3 (measured mean), 5.5 (WORKING BUDGET), 7.5 (measured max)},
  recalibrated 2026-09-28 from all sixteen relaxations (the a-priori 6/10/15 assumption is retired).

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_dataset_plan.py
"""
import csv

# Measured on ALL 16 relaxations (2026-09-27/28, 28-151 atoms, EDIFFG=-0.02, 3-57 ionic steps), wall / model single point:
# Batch A: 2.2, 2.4, 2.7, 1.0, 1.3, 3.9, 3.7, 7.3, 5.4;  Batch B: 2.9, 5.4, 4.9, 7.2, 5.6, 6.2, 6.9  -> mean 4.3x, median 4.4x, max 7.3x.
F_RELAX = {"4.3x (measured mean, n=16)": 4.3, "5.5x (working budget)": 5.5, "7.5x (measured max)": 7.5}
F_WORK = 5.5
MU_ENDPOINTS = ["-5.1071", "-4.7071"]


def c_min(n):
    return 0.55 * n * max(1.0, n / 72.0) ** 0.5


rows = list(csv.DictReader(open("dataset_plan_v1.csv")))
for r in rows:
    for k in ("n_atoms", "n_ref", "n_relax", "n_perturb", "n_collective", "n_path", "n_reusable_states"):
        r[k] = int(r[k])
    r["mu"] = [m for m in r["mu_list"].split(";") if m]
    r["n_mu_ref"] = len(r["mu"])
    r["n_var"] = r["n_perturb"] + r["n_collective"] + r["n_path"]
    r["n_mu_var"] = (2 if r["n_atoms"] > 100 else len(r["mu"])) if r["n_var"] else 0
    r["mu_var"] = (MU_ENDPOINTS if r["n_atoms"] > 100 else r["mu"]) if r["n_var"] else []
    r["n_configs"] = r["n_ref"] + r["n_var"]
    r["n_states_target"] = r["n_ref"] * r["n_mu_ref"] + r["n_var"] * r["n_mu_var"]
    r["n_new_dft"] = r["n_states_target"] - r["n_reusable_states"]
    r["c_min"] = c_min(r["n_atoms"])
    r["cost_sp_h"] = r["n_new_dft"] * r["c_min"] / 60.0
    r["cost_relax_h"] = {k: r["n_relax"] * f * r["c_min"] / 60.0 for k, f in F_RELAX.items()}
    r["cost_work_h"] = r["cost_sp_h"] + r["n_relax"] * F_WORK * r["c_min"] / 60.0

fam_order = []
for r in rows:
    if r["family"] not in fam_order:
        fam_order.append(r["family"])


def tot(key, sel=lambda r: True):
    return sum(r[key] for r in rows if sel(r))


def tot_relax(k, sel=lambda r: True):
    return sum(r["cost_relax_h"][k] for r in rows if sel(r))


L = []
L += ["# Complete Au surface dataset — computation plan v1 (rev 2, 2026-09-27)", "",
      "Generated from `dataset_plan_v1.csv` by `scripts/build_dataset_plan.py`; edit the CSV, not this file. Scope and sampling "
      "weights are frozen here once; execution is batched by cluster limits afterwards. No DFT is submitted by this document. "
      "Rev 2 applies the user's decisions of 2026-09-27: nine families kept, two-scale lean sampling, relaxation cost shown at "
      "4.3/5.5/7.5 single points (measured mean / working budget / measured max, recalibrated from all sixteen relaxations; the a-priori 6/10/15 assumption is retired), and four definition fixes (vicinal faces rebuilt from the correct "
      "basis, atom-conserving step path, R1/R2 registry and layer count, kink/C1 naming and periodicity).", "",
      "## 1. Scope", "",
      "**Included families** (each with main sampling and a few references, all in one table): flat Au(111); point defects / minimal "
      "clusters; straight strip steps; crystallographic vicinal step faces; kinks / edge rearrangement; single-layer islands; "
      "single-layer pits; reconstruction-related stacking environments; two composite morphologies as references. No further family is added in v1.", "",
      "**Explicitly excluded**: arbitrary grain boundaries, nanoparticles, multilayer spikes/towers, other low-index faces "
      "(Au(100) phase-2 stays recorded in defect_plan.md, not budgeted here), the full surface phase diagram, explicit water, "
      "any Cl species, MD trajectories. The 22×√3 herringbone cell itself (~63 Å period) is not built; R2 is a constrained stripe/domain-wall approximation.", "",
      "**Model, potentials, label standard**: fixed 4-layer (111) slabs (bottom two layers fixed) or ~9 Å-thick vicinal slabs, a = 4.158 Å, "
      "single-sided 1 M implicit electrolyte, static CP-DFT single points; production configuration (PREC=Normal, ALGO=Fast, NPAR=16 hybrid, "
      "FERMICONVERGE=0.01, mixing untouched) is the label standard; three common internal electron chemical potentials μ_e = −5.1071, −4.9071, "
      "−4.7071 eV (internal μ₀ reference, not V vs RHE); states are labelled by the actual converged μ_e/N_e. PREC=Accurate pilot results stay a "
      "separate reference group and are **not** counted as reusable coverage. The 33 rows are plan entries; the final dataset counts actual "
      "coordinates and states — duplicates that appear after building or relaxing are merged and recorded, not kept to preserve a count.", "",
      "## 2. Sampling rule (frozen)", "",
      "| item | decision |", "|---|---|",
      "| reference configs (n_ref = 2 for main rows) | the ideal bulk-truncated geometry and the CP-relaxed geometry; the second one *is* the product of n_relax = 1, not an extra relaxation |",
      "| relaxation | one per main row at μ_e = −4.9071 eV, movable atoms free; if the relaxed state leaves the initial registry (A1-hcp, R1) it is recorded under its actual geometry, and coincidences with other rows are de-duplicated and logged |",
      "| random perturbations | two per main row: one configuration at σ = 0.05 Å and one at σ = 0.10 Å (Gaussian, movable atoms, fixed seeds), applied to the relaxed geometry |",
      "| collective deformations | kept where the table already had them (2): top interlayer spacing −3 %; in-plane strain +1 % (strips: edge-row bend ±0.15 Å instead); rows with 0 stay at 0 |",
      "| path images | kept at the low counts in the table; images are constructed non-equilibrium geometries by linear interpolation on the movable atoms with **atom number and mapping conserved**; where the end point has no row of its own, the last image *is* the end point (unrelaxed) — no hidden end-point relaxations |",
      "| potentials, reference configs | all three μ_e for every row that lists them |",
      "| potentials, other configs | three μ_e if N ≤ 100; the two end points −5.1071 and −4.7071 eV if N > 100 |",
      "| reference rows (tier = reference) | the one to three states already listed, no expansion |",
      "| relaxed final step as a state | reused as the relaxed-reference single point only if its energy, forces and every grid belong to the same final geometry; otherwise a static export is run — never counted twice, never omitted |", "",
      "## 3. Master table", "",
      "| ID | family | tier | coordinates | N | movable | ref | relax | pert | coll | path | μ_e ref configs | μ_e other configs | states | reusable | new DFT | single-point cost (node·h) | purpose |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for fam in fam_order:
    for r in [x for x in rows if x["family"] == fam]:
        L.append(f"| {r['structure_id']} | {r['family']} | {r['tier']} | {r['coordinate_source']} | {r['n_atoms']} | {r['movable_atoms']} | "
                 f"{r['n_ref']} | {r['n_relax']} | {r['n_perturb']} | {r['n_collective']} | {r['n_path']} | {', '.join(r['mu'])} | "
                 f"{', '.join(r['mu_var']) if r['mu_var'] else '—'} | {r['n_states_target']} | {r['n_reusable_states']} | {r['n_new_dft']} | "
                 f"{r['cost_sp_h']:.1f} | {r['purpose']} |")

L += ["", "## 4. Totals and budget", "",
      "| | structures | configs | target states | reusable | new DFT single points | relaxations | single-point cost (node·h) | + relax 4.3× (measured mean, n=16) | + relax **5.5× (working)** | + relax 7.5× (measured max) |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
for name, sel in (("main", lambda r: r["tier"] == "main"), ("reference", lambda r: r["tier"] == "reference"), ("**all**", lambda r: True)):
    sp = tot("cost_sp_h", sel)
    L.append(f"| {name} | {sum(1 for r in rows if sel(r))} | {tot('n_configs', sel)} | {tot('n_states_target', sel)} | {tot('n_reusable_states', sel)} | "
             f"{tot('n_new_dft', sel)} | {tot('n_relax', sel)} | {sp:.0f} | {sp + tot_relax('4.3x (measured mean, n=16)', sel):.0f} | "
             f"**{sp + tot_relax('5.5x (working budget)', sel):.0f}** | {sp + tot_relax('7.5x (measured max)', sel):.0f} |")
W = tot("cost_work_h")
L += ["", "Per family (working budget, relaxation at 5.5×):", "", "| family | structures | new DFT single points | relaxations | cost (node·h) | share |", "|---|---|---|---|---|---|"]
for fam in fam_order:
    sel = lambda r, f=fam: r["family"] == f
    L.append(f"| {fam} | {sum(1 for r in rows if sel(r))} | {tot('n_new_dft', sel)} | {tot('n_relax', sel)} | {tot('cost_work_h', sel):.0f} | {tot('cost_work_h', sel)/W*100:.0f} % |")
big = lambda r: r["n_atoms"] > 200
L += ["", f"The four > 200-atom references (Island-19-8x8, Island-7-8x8, Pit-19-8x8, Pit-7-8x8) cost {tot('cost_work_h', big):.0f} node·h "
      f"({tot('cost_work_h', big)/W*100:.0f} % of the working budget) for {tot('n_new_dft', big)} states; they are kept, without perturbations or paths.", "",
      "**Reading the numbers.** Node·h are highmem node-hours (16 MPI × 8 OpenMP per node). The single-point cost model is fitted on 36–72-atom "
      "runs only; its N^1.5 growth for 108–275 atoms and the relaxation factor are assumptions — the 5.5× column is the working budget, 4.3× and 7.5× "
      "are the measured mean and maximum over all sixteen relaxations, none is an upper bound (failed/re-started runs are not included). Highmem is billed at 4 SU per "
      f"core-hour of *requested* resources (audit §F): working budget ≈ {W*128*4/1000:.0f} k SU, range {(tot('cost_sp_h')+tot_relax('4.3x (measured mean, n=16)'))*128*4/1000:.0f}–"
      f"{(tot('cost_sp_h')+tot_relax('7.5x (measured max)'))*128*4/1000:.0f} k SU. With the 2-concurrent-job limit the *ideal full-load* time is node·h/48 days "
      f"(≈ {W/48:.0f} days at the working budget) — before queueing, dependencies (relax → perturb) and re-runs; it is not a completion promise. "
      "Warm-starting perturbed configurations from the parent CHGCAR is not assumed anywhere (not measured).", "",
      "**Recalibration**: the first completed Batch-A relaxation of each size class and the first > 100-atom single points update `c_min` and "
      "`F_RELAX` in this script and regenerate this document; scope and weights do not change with them.", "",
      "## 5. Structures that still have to be built, and the four corrected definitions", "",
      "All builds use the existing generators' conventions (a = 4.158 Å; (111) slabs: 4 layers at z = 5.0/7.4/9.8/12.2 Å; vicinal slabs: ~9 Å metal "
      "thickness like Au211; Lz = 44.603 Å, SOL_Z0 = 9.801, SOL_Z1 = 34.603 Å, DIPOL 0.5 0.5 0.5, k-mesh by the N·|a| ≈ 35 Å rule) and pass the "
      "existing checks (minimum Au–Au distance ≥ 2.6 Å after perturbation, ≥ 2.9 Å for ideal builds; periodic-image isolation; registry check of every added layer).", "",
      "| ID | build recipe | validation / acceptance |", "|---|---|---|",
      "| **Au221 / Au332 / Au554 (DONE, corrected)** | `scripts/build_vicinal_fixed.py`: fcc primitive cell with indices transformed to the primitive basis ((2,2,1)→(3,3,4), (3,3,2)→(5,5,6), (5,5,4)→(9,9,10)), Gauss-reduced surface cell, 14/21/36 atomic planes (0.693/0.443/0.256 Å) for ~9 Å thickness; Au221 repeated ×2 along the 2.94 Å step vector | surface normal recovered from the coordinates ∥ (h,k,l), plane spacing = a/(2√(h²+k²+l²)), min distance 2.940 Å — recorded in `report_assets/batch1/vicinal_rebuild_validation.json`. The previous files were cubic (113), a second (211) cell and (223): kept as `Au113_retired`, `Au211b_retired`, `Au223_retired`, not in the plan; the old generator block is disabled with a note |",
      "| Flat-8x2 (DONE 2026-09-27) | remove the 8 strip atoms from Step-8x2 (as done for Flat-16x1) | base atoms unchanged (bijective match), 64 atoms |",
      "| Kink-edge1 / Kink-edge2 (DONE 2026-09-27, 109 atoms each, kink CN 5) | Step-8x3 (base 96, strip 12) plus one extra edge atom on edge1 or edge2 in one of the three rows; ny = 3 is the minimum period for one kink per period and is **not** reducible | added atom in a hollow of the layer below (registry); the actual edge contour is traced and the number of corners per period recorded — the name states the edge only, no A/B assignment |",
      "| Step-8x2_edge-vacancy_plus_foot-adatom (DONE 2026-09-27) | Step-8x2 with one edge2 strip atom removed and placed in the fcc hollow of the lower terrace adjacent to edge2 — **72 atoms** | registry, min distance; same cell, atom count and atom mapping as Step-8x2 so the 3 detachment images interpolate between the two rows |",
      "| Island-7-elongated (DONE 2026-09-27, image separation 8.82 Å) | 7 atoms as a 3+4 two-row zigzag chain in fcc hollows on the 6x6 slab (same cell as Island-7-6x6) | connected island (every atom ≥ 2 island neighbours), image separation ≥ 8 Å |",
      "| Pit-7-trench (DONE 2026-09-27) | 7 vacancies forming a 3+4 two-row zigzag trench in the 6x6 top layer | connected vacancy footprint, image separation as Pit-7 |",
      "| R1-hcp-terminated (DONE 2026-09-27; relaxed state stays hcp) | T-4x4 with the whole top layer translated by the fcc→hcp registry vector taken from the slab itself (difference between the fcc and hcp hollow positions of the second layer; length a/√6 = 1.70 Å, an in-plane a/6⟨112⟩-type direction — **not** a/√3 [1̄10]) | every top atom in a threefold hollow of layer n−1 **and** vertically above an atom of layer n−2 (hcp stacking); no atom atop a layer n−1 atom; min distance 2.94 Å |",
      "| R2-stripe-wall (DONE 2026-09-27, 65 atoms, 5 fcc / 4 hcp / 8 wall) | 16x1 cell; **3 base layers (48 atoms) + a 17-atom top layer** on the 16-site row = 65 atoms, 4 layers in total; two domains (fcc, hcp registry) joined by two walls of ~4 atoms along a1 carrying the 1/16 compression, positions from a 1-D soliton profile; constrained relaxation follows | min distance ≥ 2.6 Å before relaxation, fcc/hcp registry in the domain centres, periodic closure of the displacement field (total shift = one row spacing); described as a constrained domain-wall approximation, not as the herringbone |",
      "| C1-island-near-step (DONE 2026-09-27, 151 atoms) | **Step-8x4** strip (base 128 + strip 16 = 144) + compact 7-atom island on the lower terrace = 151 atoms (Step-8x3 rejected: a 7-atom hexagon spans 3 rows, so in a 3-row period it touches its own image) | built: the 10 Å lower terrace has 4 rows and the hexagon 3, so attachment to one foot row is unavoidable — placed on the edge2 side, 3 island–strip contacts → **step-attached**; island–image separation 5.88 Å (one vacant site) |",
      "| C2-island+pit (DONE 2026-09-27, **adjusted: 8x8 cell, 256 atoms**) | the 6x6 cell cannot hold two 3-row features with a 2-row rim gap (3+2+3 > 6 rows); built in the 8x8 slab: compact island rows 0–2, compact pit rows 4–6, one vacant row between the rims on both sides, 4-site a2 offset | 256 atoms net; rim-to-rim distance 6.12 Å; island/pit image separation 17.6 Å (plan adjustment under §6 rule a, cost re-estimated) |",
      "| path images | A1: 1 image (bridge site); Step-8x2: 3 images between the relaxed strip and the Step-8x2_edge-vacancy_plus_foot-adatom geometry; kinks: 2 images of the kink atom moving one site along the edge (image 2 = translated end point); Island-7 / Pit-7: 2 images of one rim atom leaving / entering (image 2 = end point) | interpolated on movable atoms only; atom number, order and mapping identical at both ends; no end point is relaxed unless it has its own row |", "",
      "## 6. Execution batches (resource management only — the frozen scope and weights do not change)", "",
      "1. **Batch A (≤ 72 atoms, anchors)**: T-4x4, V1, A1-fcc, A1-hcp, Au221, Au211, Flat-16x1 (+1), Flat-8x2, Step-8x1 (+2), Step-16x1 (+1), Step-8x2, R1, R2, V2/V3/A3, Au332/Au554, Step-8x2_edge-vacancy_plus_foot-adatom. Relaxations first (they recalibrate the relaxation factor), then perturbations/collective/paths.",
      "2. **Batch B (73–151 atoms)**: Step-24x1, Kink-edge1/edge2, Pit-7-compact, Pit-7-trench, Step-16x2, C2, Island-7-compact, Island-7-elongated, C1. The first completed single point of each size recalibrates the cost model.",
      "3. **Batch C (> 200 atoms, references)**: Pit-19-8x8, Pit-7-8x8, Island-7-8x8, Island-19-8x8 — one or two states each; submitted last.", "",
      "Adjustments after freezing are made only for (a) a build that fails its acceptance test, (b) an actual duplicate found after building/relaxing, "
      "(c) an explicit budget cut applied globally to perturbation counts — each recorded in the CSV with the reason. Batches are not decision points.", "",
      "## 7. Not decided or not required here", "",
      "- The trainer's energy target and input keys (settled against the training code on its own machine); every state carries TOTEN, E_without_entropy, E(σ→0), GCE, forces, actual N_e/μ_e and the indexed 3D fields.",
      "- Model interface tests or small-sample training are not prerequisites for starting Batch A.",
      "- The closed step-mechanism study is not reopened; Step-16x1/8x1 results enter only as reusable states. The old-fork interface example in `dataset_v0/mace_input/` is unrelated to this production plan.",
      "- Follow-up outside this plan: the structure-gallery figures for the three replaced vicinal faces still show the retired geometries and need re-rendering.", ""]
open("dataset_plan_v1.md", "w").write("\n".join(L))
print(f"rows={len(rows)}  main={sum(1 for r in rows if r['tier']=='main')}  reference={sum(1 for r in rows if r['tier']=='reference')}")
print(f"configs={tot('n_configs')}  target states={tot('n_states_target')}  reusable={tot('n_reusable_states')}  new DFT single points={tot('n_new_dft')}  relaxations={tot('n_relax')}")
sp = tot("cost_sp_h")
print(f"single points {sp:.0f} node-h; +relax 4.3x {sp+tot_relax('4.3x (measured mean, n=16)'):.0f}; 5.5x (working) {sp+tot_relax('5.5x (working budget)'):.0f}; 7.5x {sp+tot_relax('7.5x (measured max)'):.0f}")
for fam in fam_order:
    sel = lambda r, f=fam: r["family"] == f
    print(f"  {fam:28s} new={tot('n_new_dft', sel):4d}  relax={tot('n_relax', sel)}  work={tot('cost_work_h', sel):6.0f} h")
