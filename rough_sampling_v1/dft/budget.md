# Rough CP-DFT budget (estimate, nothing submitted)

200 single points. Production standard K.8 from scripts/production.py; one node (16 MPI x 8 OpenMP) per task on `wholenode`.
NELECT start guess from C = 11.59 uF/cm2, U_pzc = +0.0056 V (dataset-wide medians; the per-geometry PZC spread is 254 mV, so the guess is off by up to ~0.5 e in an 8x8 cell).

Measured cost basis (dataset_v1, 231-300 atoms, minutes per single point): cold start median 1011 / max 1346 (13 states); NELECT-seeded median 399 / max 563 (4 states, not a controlled comparison). Scaled by (n/275)^1.5 above that range: an ASSUMPTION.

| cell | tasks | atoms | node-h if seeded median | if seeded max | if cold median | if cold max | walltime requested |
|---|---|---|---|---|---|---|---|
| 10x8 | 23 | 284-383 | 213 | 300 | 539 | 718 | 29:00:00-45:00:00 |
| 8x8 | 177 | 202-318 | 1161 | 1638 | 2942 | 3917 | 17:00:00-34:00:00 |
| **all** | 200 | | **1374** | 1938 | 3481 | **4634** | |

For scale: dataset_v1 (505 states) cost 757 node-h in total.
Disk: ~11 GB per state with all fields (measured on the 8x8 references) -> ~2.2 TB for 200 states; scratch usage 1.9 TB of 100 TB (2026-10-02).

Plan: submit 8 first (`rough_dft.py submit --first 8 --confirm`), recalibrate COST and the walltime from them, then the remaining 192.
