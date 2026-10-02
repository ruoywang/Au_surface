# Dataset v1 (production states) — summary

Generated 2026-10-01 06:57 from 05_production/queue.json. **505 accepted states** (107 distinct geometries, 33 structures), 0 skipped. Labels are by the actual converged state; no `energy` key is assigned (see dataset_v0/README.md).

| family | states |
|---|---|
| point defect | 110 |
| strip step | 78 |
| vicinal step face | 70 |
| kink / edge rearrangement | 57 |
| single-layer pit | 51 |
| single-layer island | 51 |
| reconstruction-related | 40 |
| flat Au(111) | 40 |
| composite | 8 |

| config type | states |
|---|---|
| ideal | 157 |
| pert05 | 73 |
| pert10 | 73 |
| coll | 70 |
| relaxed | 64 |
| path | 52 |
| relax | 16 |

- μ_e more than 5 meV from TARGETMU (inside FERMICONVERGE = 0.01): 89 states; max |Δμ| = 9.9 meV.
- |total drift_z|: median 0.312, max 2.900 eV/Å (PREC=Normal aliasing; see dataset_v0/README.md).
- relaxations included as final-geometry states: 16 (all with 'reached required accuracy'); 16 queue records are aliases of an exported state (the reused relaxed@-4.9071 record = the relax task's final step) and are listed under `aliases`, not exported twice (before 2026-09-28 20:50 they were).
- campaigns: {'dataset_plan_v1_rev2': 291, 'mu05_extension': 214}; TARGETMU coverage: {'-4.4071': 107, '-4.7071': 101, '-4.9071': 85, '-5.1071': 105, '-5.4071': 107}; states per geometry: {3: 2, 4: 26, 5: 79}.
- QC (mu closure <= 0.011 eV, CHGCAR written after task start, 15 field files, geometry = registry within 0.005 Å, relaxation accuracy): 504 accepted, 1 flagged: C2-island+pit__ideal__mu-5.4071 {'fields_complete': False}
- geometry deviation from the registry source: max 0.0030 Å (3 states above 1e-4 Å, the pilot-reused ideal states of Step-8x1/Step-16x1); states in the relaxation directory (relax task and reused relaxed@-4.9071) carry the CONTCAR geometry.
- fields: every state indexes CHGCAR, LOCPOT, PHI, PHISOLV, VSOLV, RHOB, RHOION, ELOC, P, SVDW, SION, SSOLV, SCAV, SDIEL, POT in its run directory (not copied).
