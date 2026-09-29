# Dataset v1 (production states) — summary

Generated 2026-09-29 17:48 from 05_production/queue.json. **439 accepted states** (107 distinct geometries, 33 structures), 0 skipped. Labels are by the actual converged state; no `energy` key is assigned (see dataset_v0/README.md).

| family | states |
|---|---|
| point defect | 110 |
| vicinal step face | 70 |
| strip step | 70 |
| kink / edge rearrangement | 51 |
| reconstruction-related | 40 |
| flat Au(111) | 40 |
| single-layer pit | 27 |
| single-layer island | 27 |
| composite | 4 |

| config type | states |
|---|---|
| ideal | 135 |
| coll | 70 |
| pert05 | 61 |
| pert10 | 61 |
| relaxed | 52 |
| path | 44 |
| relax | 16 |

- μ_e more than 5 meV from TARGETMU (inside FERMICONVERGE = 0.01): 70 states; max |Δμ| = 9.9 meV.
- |total drift_z|: median 0.270, max 2.899 eV/Å (PREC=Normal aliasing; see dataset_v0/README.md).
- relaxations included as final-geometry states: 16 (all with 'reached required accuracy'); 16 queue records are aliases of an exported state (the reused relaxed@-4.9071 record = the relax task's final step) and are listed under `aliases`, not exported twice (before 2026-09-28 20:50 they were).
- campaigns: {'dataset_plan_v1_rev2': 291, 'mu05_extension': 148}; TARGETMU coverage: {'-4.4071': 74, '-4.7071': 101, '-4.9071': 85, '-5.1071': 105, '-5.4071': 74}; states per geometry: {1: 2, 2: 20, 3: 11, 4: 6, 5: 68}.
- QC (mu closure <= 0.011 eV, CHGCAR written after task start, 15 field files, geometry = registry within 0.005 Å, relaxation accuracy): 439 accepted, 0 flagged.
- geometry deviation from the registry source: max 0.0030 Å (3 states above 1e-4 Å, the pilot-reused ideal states of Step-8x1/Step-16x1); states in the relaxation directory (relax task and reused relaxed@-4.9071) carry the CONTCAR geometry.
- fields: every state indexes CHGCAR, LOCPOT, PHI, PHISOLV, VSOLV, RHOB, RHOION, ELOC, P, SVDW, SION, SSOLV, SCAV, SDIEL, POT in its run directory (not copied).
