# Dataset v1 (production states) — summary

Generated 2026-09-28 from 05_production/queue.json. **225 accepted states** (90 distinct geometries, 25 structures), 0 skipped. Labels are by the actual converged state; no `energy` key is assigned (see dataset_v0/README.md).

| family | states |
|---|---|
| point defect | 69 |
| vicinal step face | 44 |
| strip step | 36 |
| reconstruction-related | 26 |
| flat Au(111) | 25 |
| kink / edge rearrangement | 19 |
| single-layer pit | 4 |
| single-layer island | 2 |

| config type | states |
|---|---|
| ideal | 63 |
| coll | 42 |
| relaxed | 33 |
| pert05 | 27 |
| pert10 | 27 |
| path | 18 |
| relax | 15 |

- μ_e more than 5 meV from TARGETMU (inside FERMICONVERGE = 0.01): 20 states; max |Δμ| = 9.6 meV.
- |total drift_z|: median 0.248, max 1.600 eV/Å (PREC=Normal aliasing; see dataset_v0/README.md).
- relaxations included as final-geometry states: 15 (all with 'reached required accuracy').
- fields: every state indexes CHGCAR, LOCPOT, PHI, PHISOLV, VSOLV, RHOB, RHOION, ELOC, P, SVDW, SION, SSOLV, SCAV, SDIEL, POT in its run directory (not copied).

## Force-label quality: egg-box error of the production configuration (found 2026-09-28)

`total_drift` (sum of all forces, should vanish) is −1.26 to −1.60 eV/Å along z for every ideal single point in the
4×4 cell (T, V1, V2, V3, A1, A3, R1: 0.02 eV/Å per atom, coherent), −0.27 in the 16×1 cells, −0.5 to −0.7 for
perturbed/relaxed 4×4 configurations, −0.7 for the kink cells and −1.0 for the 6×6 pit/island cells. Cause, verified
on the OUTCARs: PREC=Normal puts exactly 16.000 coarse FFT points per 2.94 Å atom spacing in the 4×4 cell, so the
LREAL=Auto real-space projection ("egg-box") error is identical for every atom in a layer and adds up; in the 16×1
cell (16.875 points) the phases differ and the error largely cancels in the sum; PREC=Accurate (22.5 points, tighter
ROPT) gives 0.017 eV/Å total. The per-atom error is therefore of order 0.02 eV/Å (z) throughout the production set,
visible in the drift only where the cell is commensurate. Options (none applied here): (a) subtract drift/N per atom
— removes the coherent part only; (b) quantify directly by recomputing 2–3 states with LREAL=.FALSE. (~1–2 node·h);
(c) change the label standard (ADDGRID or LREAL=.FALSE.) for a future relabelling. The v1 labels are exported as
computed, with `total_drift_eV_per_A` per state.
