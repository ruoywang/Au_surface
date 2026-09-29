# Dataset v1 (production states) — summary

Generated 2026-09-28 from 05_production/queue.json. **268 accepted states** (97 distinct geometries, 26 structures), 0 skipped. Labels are by the actual converged state; no `energy` key is assigned (see dataset_v0/README.md).

| family | states |
|---|---|
| point defect | 66 |
| vicinal step face | 42 |
| strip step | 40 |
| kink / edge rearrangement | 31 |
| reconstruction-related | 24 |
| flat Au(111) | 24 |
| single-layer pit | 23 |
| single-layer island | 18 |

| config type | states |
|---|---|
| ideal | 68 |
| coll | 42 |
| pert05 | 41 |
| pert10 | 41 |
| relaxed | 32 |
| path | 28 |
| relax | 16 |

- μ_e more than 5 meV from TARGETMU (inside FERMICONVERGE = 0.01): 20 states; max |Δμ| = 9.6 meV.
- |total drift_z|: median 0.313, max 2.899 eV/Å (PREC=Normal aliasing; see dataset_v0/README.md).
- relaxations included as final-geometry states: 16 (all with 'reached required accuracy'); 16 queue records are aliases of an exported state (the reused relaxed@-4.9071 record = the relax task's final step) and are listed under `aliases`, not exported twice (before 2026-09-28 20:50 they were).
- campaigns: {'dataset_plan_v1_rev2': 268}; TARGETMU coverage: {'-4.7071': 97, '-4.9071': 75, '-5.1071': 96}; states per geometry: {2: 23, 3: 74}.
- QC (mu closure <= 0.011 eV, CHGCAR written after task start, 15 field files, geometry = registry within 0.005 Å, relaxation accuracy): 268 accepted, 0 flagged.
- geometry deviation from the registry source: max 0.0030 Å (3 states above 1e-4 Å, the pilot-reused ideal states of Step-8x1/Step-16x1); states in the relaxation directory (relax task and reused relaxed@-4.9071) carry the CONTCAR geometry.
- fields: every state indexes CHGCAR, LOCPOT, PHI, PHISOLV, VSOLV, RHOB, RHOION, ELOC, P, SVDW, SION, SSOLV, SCAV, SDIEL, POT in its run directory (not copied).
