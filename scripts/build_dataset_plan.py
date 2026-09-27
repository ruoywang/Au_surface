#!/usr/bin/env python3
"""Generate dataset_plan_v1.md from dataset_plan_v1.csv (the single source of truth for the
complete Au dataset computation plan). Every multiplier is a CSV column; totals and costs are
computed here, never hand-summed.

Per row:
  n_configs         = n_ref + n_perturb + n_collective + n_path
  n_states_target   = n_configs * n_mu
  n_new_dft         = n_states_target - n_reusable_states           (reuse = same production standard only)
  cost_node_h       = [ n_new_dft * c(N) + n_relax * RELAX_FACTOR * c(N) ] / 60
Cost model c(N) [min per static CP single point, 1 highmem node = 16 MPI x 8 OpenMP, production config]:
  c(N) = 0.55 * N * max(1, N/72)**0.5      fitted to measured 36 -> 17.9 min, 64 -> 27-39 min, 72 -> 40-43 min;
  the N^1.5 growth beyond 72 atoms is an assumption (untested above 72 atoms in this config).
  A CP relaxation is charged as RELAX_FACTOR = 6 single points (15-20 ionic steps with warm-started SCF).
Scenario "lean" (computed alongside): perturbation configs halved (floor, min 1) for main rows, and
  rows with N > 100 use 2 potentials instead of 3 for the perturbation/path configs.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_dataset_plan.py
"""
import csv
import math

RELAX_FACTOR = 6.0
MU_MAIN = ["-5.1071", "-4.9071", "-4.7071"]


def c_min(n):
    return 0.55 * n * max(1.0, n / 72.0) ** 0.5


rows = list(csv.DictReader(open("dataset_plan_v1.csv")))
for r in rows:
    for k in ("n_atoms", "n_ref", "n_relax", "n_perturb", "n_collective", "n_path", "n_reusable_states"):
        r[k] = int(r[k])
    r["mu"] = [m for m in r["mu_list"].split(";") if m]
    r["n_mu"] = len(r["mu"])
    r["n_configs"] = r["n_ref"] + r["n_perturb"] + r["n_collective"] + r["n_path"]
    r["n_states_target"] = r["n_configs"] * r["n_mu"]
    r["n_new_dft"] = r["n_states_target"] - r["n_reusable_states"]
    r["c_min"] = c_min(r["n_atoms"])
    r["cost_h"] = (r["n_new_dft"] * r["c_min"] + r["n_relax"] * RELAX_FACTOR * r["c_min"]) / 60.0
    # lean scenario
    n_pert_lean = max(1, r["n_perturb"] // 2) if r["n_perturb"] > 0 else 0
    n_mu_lean_pert = 2 if r["n_atoms"] > 100 else r["n_mu"]
    ref_states = r["n_ref"] * r["n_mu"]
    var_states = (n_pert_lean + r["n_collective"] + r["n_path"]) * n_mu_lean_pert
    r["n_states_lean"] = ref_states + var_states
    r["n_new_lean"] = r["n_states_lean"] - r["n_reusable_states"]
    r["cost_lean_h"] = (r["n_new_lean"] * r["c_min"] + r["n_relax"] * RELAX_FACTOR * r["c_min"]) / 60.0

fam_order = []
for r in rows:
    if r["family"] not in fam_order:
        fam_order.append(r["family"])


def tot(key, sel=lambda r: True):
    return sum(r[key] for r in rows if sel(r))


L = []
L += ["# Complete Au surface dataset — computation plan v1 (2026-09-27)", "",
      "Generated from `dataset_plan_v1.csv` by `scripts/build_dataset_plan.py`; edit the CSV, not this file. "
      "Scope is fixed here once; execution is batched by cluster limits afterwards. No DFT is submitted by this document.", "",
      "## 1. Scope", "",
      "**Included families** (each with main sampling and a few references, all in one table): flat Au(111); point defects / minimal "
      "clusters; straight strip steps; crystallographic vicinal step faces; kinks / edge rearrangement; single-layer islands; "
      "single-layer pits; reconstruction-related stacking environments; two composite morphologies as references.", "",
      "**Explicitly excluded**: arbitrary grain boundaries, nanoparticles, multilayer spikes/towers, other low-index faces "
      "(Au(100) phase-2 stays recorded in defect_plan.md, not budgeted here), the full surface phase diagram, explicit water, "
      "any Cl species, MD trajectories. The 22×√3 herringbone cell itself (~63 Å period) is not built; R2 is its stripe/soliton approximation.", "",
      "**Model, potentials, label standard**: fixed 4-layer slabs (bottom two layers fixed), a = 4.158 Å, single-sided 1 M implicit "
      "electrolyte, static CP-DFT single points; production configuration (PREC=Normal, ALGO=Fast, NPAR=16 hybrid, FERMICONVERGE=0.01) "
      "is the label standard; three common internal electron chemical potentials μ_e = −5.1071, −4.9071, −4.7071 eV (internal μ₀ "
      "reference, not V vs RHE); states are labelled by the actual converged μ_e/N_e. PREC=Accurate pilot results stay a separate "
      "reference group and are **not** counted as reusable coverage.", "",
      "## 2. Sampling template (what one row means)", "",
      "| column | definition |", "|---|---|",
      "| n_ref | reference configurations: the ideal bulk-truncated geometry, plus the CP-relaxed geometry when n_relax = 1 |",
      "| n_relax | CP relaxations at μ_e = −4.9071 (movable atoms free, bottom two layers fixed), charged as 6 single points each; produce the relaxed reference |",
      "| n_perturb | random-displacement configurations of the relaxed geometry on the movable atoms: Gaussian σ = 0.05 Å (seeds 1,2) and 0.10 Å (seeds 3,4); halved sets use seeds 1,3 |",
      "| n_collective | collective deformations of the relaxed geometry: top interlayer spacing −3 %; in-plane strain +1 % (cell scaled); for strips the second one is an edge-row bend of ±0.15 Å along the period |",
      "| n_path | local rearrangement images: linear interpolation between the relaxed initial and final geometries of the named move, end points excluded (they are their own rows or the relaxed reference) |",
      "| mu_list | electron chemical potentials computed for every configuration of the row |",
      "| n_states_target | n_configs × n_mu — the number of labelled states the row contributes |",
      "| n_reusable_states | states already computed in the production standard with identical cell/window/k-mesh |",
      "| cost | highmem node-hours from the measured production timings (see cost model in the script header) |", "",
      "Perturbed and path configurations are single points on non-equilibrium geometries; none of them is relaxed back, so no near-duplicate minima are generated.", "",
      "## 3. Master table", "",
      "| ID | family | tier | coordinates | N | movable | ref | relax | pert | coll | path | μ_e (eV) | states | reusable | new DFT | cost (node·h) | purpose |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for fam in fam_order:
    for r in [x for x in rows if x["family"] == fam]:
        L.append(f"| {r['structure_id']} | {r['family']} | {r['tier']} | {r['coordinate_source']} | {r['n_atoms']} | {r['movable_atoms']} | "
                 f"{r['n_ref']} | {r['n_relax']} | {r['n_perturb']} | {r['n_collective']} | {r['n_path']} | {', '.join(r['mu'])} | "
                 f"{r['n_states_target']} | {r['n_reusable_states']} | {r['n_new_dft']} | {r['cost_h']:.1f} | {r['purpose']} |")

L += ["", "## 4. Totals and budget", "",
      "| | structures | configs | target states | reusable | new DFT states | relaxations | cost (highmem node·h) |", "|---|---|---|---|---|---|---|---|"]
for name, sel in (("main", lambda r: r["tier"] == "main"), ("reference", lambda r: r["tier"] == "reference"), ("**all**", lambda r: True)):
    L.append(f"| {name} | {sum(1 for r in rows if sel(r))} | {tot('n_configs', sel)} | {tot('n_states_target', sel)} | {tot('n_reusable_states', sel)} | "
             f"{tot('n_new_dft', sel)} | {tot('n_relax', sel)} | {tot('cost_h', sel):.0f} |")
L += ["", "Per family (full scenario):", "", "| family | structures | new DFT states | cost (node·h) | share |", "|---|---|---|---|---|"]
T = tot("cost_h")
for fam in fam_order:
    sel = lambda r, f=fam: r["family"] == f
    L.append(f"| {fam} | {sum(1 for r in rows if sel(r))} | {tot('n_new_dft', sel)} | {tot('cost_h', sel):.0f} | {tot('cost_h', sel)/T*100:.0f} % |")
L += ["", f"**Lean scenario** (perturbations halved on main rows; 2 potentials instead of 3 for perturbation/path configurations of rows with N > 100; "
      f"references unchanged): {tot('n_new_lean')} new DFT states, {tot('cost_lean_h'):.0f} node·h "
      f"({tot('cost_lean_h')/T*100:.0f} % of the full scenario). Rows with the largest savings: " +
      ", ".join(f"{r['structure_id']} ({r['cost_h']-r['cost_lean_h']:.0f} h)" for r in sorted(rows, key=lambda r: r["cost_lean_h"] - r["cost_h"])[:6]) + ".", "",
      "**Throughput**: highmem allows 2 running jobs per user; at 1 job = 1 node the full scenario is ≈ "
      f"{T/2/24:.0f} days of continuous submission ({tot('cost_lean_h')/2/24:.0f} days lean). Highmem is billed at 4× (audit §F); "
      f"node·h × 128 cores × 4 ≈ {T*128*4/1000:.0f} k SU full, {tot('cost_lean_h')*128*4/1000:.0f} k SU lean. Warm-starting each perturbed "
      "configuration from its parent CHGCAR (ICHARG=1) would cut SCF steps by roughly a third; not assumed in the numbers above.", "",
      "**Cost-model limits**: measured points cover 36–72 atoms only; the N^1.5 growth used for 108–275 atoms is an assumption, so the "
      "island/pit/kink rows (≈ half of the budget) carry the largest uncertainty. The first completed run of each size class recalibrates the CSV.", "",
      "## 5. Structures that still have to be built (concrete specifications)", "",
      "All builds use the existing generators' conventions (a = 4.158 Å, 4 base layers at z = 5.0/7.4/9.8/12.2 Å, Lz = 44.603 Å, "
      "SOL_Z0 = 9.801, SOL_Z1 = 34.603 Å, DIPOL 0.5 0.5 0.5, k-mesh by the N·|a| ≈ 35 Å rule) and pass the existing checks "
      "(minimum Au–Au distance ≥ 2.6 Å after perturbation, ≥ 2.9 Å for ideal builds; periodic-image isolation; registry check of every added layer).", "",
      "| ID | build recipe | validation |", "|---|---|---|",
      "| Flat-8x2 | remove the 8 strip atoms from Step-8x2 (as done for Flat-16x1) | base atoms unchanged (bijective match), 64 atoms |",
      "| Kink-A / Kink-B | Step-8x3 (base 96, strip 12) plus one extra edge atom on edge1 (A) or edge2 (B) in one of the three rows; the ny = 3 period is the minimum for one kink per period and is **not** reducible | added atom in a hollow of the layer below (registry), CN of kink atom = 5, min distance 2.94 Å |",
      "| Step-8x2+foot-adatom | Step-8x2 plus one adatom in the fcc hollow of the lower terrace adjacent to edge2 | registry, min distance |",
      "| Island-7-elongated | 7 atoms as a 3+4 two-row zigzag chain in fcc hollows on the 6x6 slab (same cell as Island-7-6x6) | connected island (every atom ≥ 2 island neighbours), image separation ≥ 8 Å |",
      "| Pit-7-trench | 7 vacancies forming a 3+4 two-row zigzag trench in the 6x6 top layer | connected vacancy footprint, image separation as Pit-7 |",
      "| R1-hcp-terminated | T-4x4 with the top layer translated by (a/√3)·[1̄10]-type shift into hcp registry | every top atom above a second-layer atom (hcp), min distance 2.94 Å |",
      "| R2-stripe-wall | 16x1 cell; top layer with 17 atoms on the 16-site row: two domains (fcc, hcp) joined by two walls of ~4 atoms each along a1, compression 1/16 spread over the walls; positions from a 1-D soliton profile, then constrained relaxation | min distance ≥ 2.6 Å before relaxation, registry fcc/hcp in the domain centres, periodic closure of the displacement field (total shift = one row spacing) |",
      "| C1-island-near-step | Step-8x3 strip (108) + compact 7-atom island on the lower terrace, nearest island atom one row from edge2 | island–edge distance recorded; no island–image contact |",
      "| C2-island+pit | 6x6 slab, compact 7-atom island and compact 7-vacancy pit with rims two rows apart | 144 atoms net; separation recorded |",
      "| path images | A1: bridge site (1); Step-8x2 detachment: 3 images between the relaxed strip and the relaxed foot-adatom structure; kinks: 2 images of the kink atom moving one site along the edge; Island-7 / Pit-7: 2 images of one rim atom leaving/entering | interpolated on movable atoms only; end points are their own rows |", "",
      "## 6. Execution batches (resource management only — the scope above does not change)", "",
      "1. **Batch A (≤ 72 atoms, anchors)**: T-4x4, V1, A1-fcc, A1-hcp, Flat-16x1 (+1), Flat-8x2, Step-8x1 (+2), Step-16x1 (+1), Step-8x2, Au211, Au221, R1, V2/V3/A3, Au332/Au554, Step-8x2+foot-adatom. Relaxations first, then perturbations/paths.",
      "2. **Batch B (73–151 atoms)**: R2, Step-24x1, Kink-A/B, C1, Pit-7-compact, Pit-7-trench, Step-16x2, C2, Island-7-compact, Island-7-elongated. The first completed run of each size recalibrates the cost model.",
      "3. **Batch C (> 200 atoms, references)**: Pit-19-8x8, Pit-7-8x8, Island-7-8x8, Island-19-8x8 — one or two states each; submitted last, cancelled first if the budget tightens.", "",
      "Batch order does not decide scope; if the budget is cut, apply the lean scenario (or a further global reduction of perturbation counts) across all rows rather than dropping a family.", "",
      "## 7. Not decided or not required here", "",
      "- The trainer's energy target and input keys (settled against the training code on its own machine, not here); this plan only guarantees that every state carries TOTEN, E_without_entropy, E(σ→0), GCE, forces, actual N_e/μ_e and the indexed 3D fields.",
      "- Model interface tests or small-sample training are not prerequisites for starting Batch A.",
      "- The closed step-mechanism study is not reopened; Step-16x1/8x1 results enter only as reusable states.", ""]
open("dataset_plan_v1.md", "w").write("\n".join(L))
print(f"rows={len(rows)}  main={sum(1 for r in rows if r['tier']=='main')}  reference={sum(1 for r in rows if r['tier']=='reference')}")
print(f"configs={tot('n_configs')}  target states={tot('n_states_target')}  reusable={tot('n_reusable_states')}  new DFT={tot('n_new_dft')}  relaxations={tot('n_relax')}")
print(f"cost full={T:.0f} node-h ({T/2/24:.0f} days at 2 jobs)   lean: new={tot('n_new_lean')} cost={tot('cost_lean_h'):.0f} node-h")
for fam in fam_order:
    sel = lambda r, f=fam: r["family"] == f
    print(f"  {fam:28s} new={tot('n_new_dft', sel):4d}  cost={tot('cost_h', sel):6.0f} h  lean={tot('cost_lean_h', sel):6.0f} h")
