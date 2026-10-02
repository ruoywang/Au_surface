# Rough states: CANDIDATE list of 8 (seed 3, frozen = False)

Nothing here is submitted or frozen. Approval = `finalize_200.py freeze`, then `rough_dft.py submit --first 8 --confirm`.

## Where they come from

| class | train | val | test |
|---|---|---|---|
| A | 0 | 0 | 2 |
| B | 0 | 0 | 2 |
| C | 0 | 0 | 2 |
| D | 0 | 0 | 2 |

Parents used: 4 of 32; cells per parent min/median/max 2/2/2 (caps: {}); cells per (parent, frame) max 2 (cap 2).
Frames (ps): 0 ps: 8

| CN stratum of the centre | states |
|---|---|
| cn<=5 | 8 |
| cn6-7 | 0 |
| cn8-9 | 0 |
| cn>=10 | 0 |

## Size, atoms, k-mesh

Cell sizes: {'10x8': 4, '8x8': 4}
Atoms: 217: 1, 262: 1, 271: 1, 273: 1, 281: 1, 340: 1, 343: 1, 351: 1
k-mesh: {'2x2x1': 8}
Size pairs (test, same U): pair1 A PA2_s51_f000000_a4390, pair2 B PB2_s51_f000000_a4279, pair3 C PC2_s51_f000000_a3486, pair4 D PD2_s51_f000000_a4171
Repair: {'no repair needed': 8}

## Morphology coverage of the chosen cells (medians and ranges)

| descriptor | min | median | max |
|---|---|---|---|
| f_cn<=5 | 0.022 | 0.055 | 0.081 |
| f_cn6-7 | 0.108 | 0.145 | 0.237 |
| f_cn8-9 | 0.300 | 0.463 | 0.584 |
| f_cn>=10 | 0.236 | 0.332 | 0.412 |
| relief_layers | 1.000 | 1.000 | 1.000 |
| n_levels | 2.000 | 2.000 | 2.000 |
| adatom_level_ML | 0.000 | 0.258 | 0.391 |
| missing_terrace_ML | 0.000 | 0.000 | 0.725 |
| exposed_pit_ML | 0.000 | 0.000 | 0.912 |
| n_islands | 0.000 | 2.000 | 4.000 |
| max_island_frac | 0.000 | 0.228 | 0.391 |
| n_pits | 0.000 | 0.000 | 1.000 |
| step_density | 0.187 | 0.217 | 0.300 |

## Potentials

| U bin (V) | states |
|---|---|
| -0.5 to -0.4 | 0 |
| -0.4 to -0.3 | 2 |
| -0.3 to -0.2 | 0 |
| -0.2 to -0.1 | 2 |
| -0.1 to +0.0 | 0 |
| +0.0 to +0.1 | 2 |
| +0.1 to +0.2 | 0 |
| +0.2 to +0.3 | 0 |
| +0.3 to +0.4 | 2 |
| +0.4 to +0.5 | 0 |

## Rejections

Centres with no passing cell: 1 of 20. Reasons over all failed trials:
- periphery has an atom less coordinated (3) than the core minimum (6) minus one: 9
- periphery has an atom less coordinated (4) than the core minimum (6) minus one: 1

Geometry duplicates dropped (occupancy patterns within 2 sites under translation, same size): 0

SHORTFALL (bucket: missing): {"A/train": 8, "A/val": 1, "B/train": 8, "B/val": 1, "C/train": 8, "C/val": 1, "D/train": 8, "D/val": 1} -- fewer eligible cells than the quota under the caps; not filled from other buckets.

Novelty of the chosen cells' centres vs the 52 old geometries (1 - max cosine): median 0.0098; rare-flagged centres chosen: 0.
