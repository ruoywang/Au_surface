# mu05 extension: +-0.5 V end points for every distinct geometry

Generated 2026-09-28 20:50 by `scripts/extend_mu05.py --report`. mu_target = mu0 - dU, mu0 = -4.9071 eV (internal reference). TARGETMU -5.4071 (dU = +0.5 V) and -4.4071 (dU = -0.5 V); existing -5.1071/-4.9071/-4.7071 kept; no intermediate points added. Static CP single points (IBRION=-1, NSW=0) on the geometry actually computed (relaxed = accepted CONTCAR). Warm start = ICHARG=1 from a copy of the same-side neighbour's CHGCAR.

Distinct geometries: 107. Extension tasks: 214 -> status {'pending': 214}. Estimated 243 node-h (measured same-structure single-point means); actual so far 0 node-h over 0 completed.

Status legend: C completed/accepted, R running, P pending, U reused, F failed, - not planned. `w:` warm start source.

| structure | config | N | family | -5.1071 | -4.9071 | -4.7071 | -5.4071 (+0.5 V) | -4.4071 (-0.5 V) | est min | actual min |
|---|---|---|---|---|---|---|---|---|---|---|
| A1-fcc | coll_strain+1pct | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 52 | 0 |
| A1-fcc | coll_topspacing-3pct | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 52 | 0 |
| A1-fcc | ideal | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 52 | 0 |
| A1-fcc | path_bridge | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 52 | 0 |
| A1-fcc | pert05 | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 52 | 0 |
| A1-fcc | pert10 | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 52 | 0 |
| A1-fcc | relaxed | 65 | point defect | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 52 | 0 |
| A1-hcp | coll_strain+1pct | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 54 | 0 |
| A1-hcp | coll_topspacing-3pct | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 54 | 0 |
| A1-hcp | ideal | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 54 | 0 |
| A1-hcp | pert05 | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 54 | 0 |
| A1-hcp | pert10 | 65 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 54 | 0 |
| A1-hcp | relaxed | 65 | point defect | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 54 | 0 |
| A3 | ideal | 67 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 52 | 0 |
| Au211 | coll_strain+1pct | 48 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 31 | 0 |
| Au211 | coll_topspacing-3pct | 48 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 31 | 0 |
| Au211 | ideal | 48 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 31 | 0 |
| Au211 | pert05 | 48 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 31 | 0 |
| Au211 | pert10 | 48 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 31 | 0 |
| Au211 | relaxed | 48 | vicinal step face | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 31 | 0 |
| Au221 | coll_strain+1pct | 28 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 22 | 0 |
| Au221 | coll_topspacing-3pct | 28 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 22 | 0 |
| Au221 | ideal | 28 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 22 | 0 |
| Au221 | pert05 | 28 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 22 | 0 |
| Au221 | pert10 | 28 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 22 | 0 |
| Au221 | relaxed | 28 | vicinal step face | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 22 | 0 |
| Au332 | ideal | 21 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 17 | 0 |
| Au554 | ideal | 36 | vicinal step face | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 40 | 0 |
| C1-island-near-step | ideal | 151 | composite | P | P | - | P (cold) | P (cold) | 195 | 0 |
| C2-island+pit | ideal | 256 | composite | P | P | - | P (cold) | P (cold) | 430 | 0 |
| Flat-16x1 | ideal | 64 | flat Au(111) | U | U | C | P (w:-5.1071) | P (w:-4.7071) | 81 | 0 |
| Flat-8x2 | ideal | 64 | flat Au(111) | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 77 | 0 |
| Island-19-8x8 | ideal | 275 | single-layer island | P | P | - | P (cold) | P (cold) | 479 | 0 |
| Island-7-8x8 | ideal | 263 | single-layer island | - | P | - | P (cold) | P (cold) | 448 | 0 |
| Island-7-compact | ideal | 151 | single-layer island | R | R | R | P (cold) | P (cold) | 286 | 0 |
| Island-7-compact | path_detach1 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-compact | path_detach2 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-compact | pert05 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-compact | pert10 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-compact | relaxed | 151 | single-layer island | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-elongated | ideal | 151 | single-layer island | R | R | R | P (cold) | P (cold) | 295 | 0 |
| Island-7-elongated | pert05 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 295 | 0 |
| Island-7-elongated | pert10 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 295 | 0 |
| Island-7-elongated | relaxed | 151 | single-layer island | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 295 | 0 |
| Kink-edge1 | ideal | 109 | kink / edge rearrangement | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 177 | 0 |
| Kink-edge1 | path_kinkmove1 | 109 | kink / edge rearrangement | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 177 | 0 |
| Kink-edge1 | path_kinkmove2 | 109 | kink / edge rearrangement | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 177 | 0 |
| Kink-edge1 | pert05 | 109 | kink / edge rearrangement | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 177 | 0 |
| Kink-edge1 | pert10 | 109 | kink / edge rearrangement | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 177 | 0 |
| Kink-edge1 | relaxed | 109 | kink / edge rearrangement | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 177 | 0 |
| Kink-edge2 | ideal | 109 | kink / edge rearrangement | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 173 | 0 |
| Kink-edge2 | path_kinkmove1 | 109 | kink / edge rearrangement | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 173 | 0 |
| Kink-edge2 | path_kinkmove2 | 109 | kink / edge rearrangement | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 173 | 0 |
| Kink-edge2 | pert05 | 109 | kink / edge rearrangement | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 173 | 0 |
| Kink-edge2 | pert10 | 109 | kink / edge rearrangement | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 173 | 0 |
| Kink-edge2 | relaxed | 109 | kink / edge rearrangement | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 173 | 0 |
| Pit-19-8x8 | ideal | 237 | single-layer pit | P | P | - | P (cold) | P (cold) | 383 | 0 |
| Pit-7-8x8 | ideal | 249 | single-layer pit | - | P | - | P (cold) | P (cold) | 413 | 0 |
| Pit-7-compact | ideal | 137 | single-layer pit | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-compact | path_rimin1 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-compact | path_rimin2 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-compact | pert05 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-compact | pert10 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-compact | relaxed | 137 | single-layer pit | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-trench | ideal | 137 | single-layer pit | R | C | C | P (w:-4.9071) | P (w:-4.7071) | 268 | 0 |
| Pit-7-trench | pert05 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 268 | 0 |
| Pit-7-trench | pert10 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 268 | 0 |
| Pit-7-trench | relaxed | 137 | single-layer pit | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 268 | 0 |
| R1-hcp-terminated | ideal | 64 | reconstruction-related | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 46 | 0 |
| R1-hcp-terminated | pert05 | 64 | reconstruction-related | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 46 | 0 |
| R1-hcp-terminated | pert10 | 64 | reconstruction-related | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 46 | 0 |
| R1-hcp-terminated | relaxed | 64 | reconstruction-related | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 46 | 0 |
| R2-stripe-wall | ideal | 65 | reconstruction-related | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 81 | 0 |
| R2-stripe-wall | pert05 | 65 | reconstruction-related | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 81 | 0 |
| R2-stripe-wall | pert10 | 65 | reconstruction-related | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 81 | 0 |
| R2-stripe-wall | relaxed | 65 | reconstruction-related | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 81 | 0 |
| Step-16x1 | ideal | 72 | strip step | U | U | C | P (w:-5.1071) | P (w:-4.7071) | 97 | 0 |
| Step-16x2 | ideal | 144 | strip step | R | R | R | P (cold) | P (cold) | 260 | 0 |
| Step-16x2 | pert05 | 144 | strip step | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 260 | 0 |
| Step-16x2 | pert10 | 144 | strip step | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 260 | 0 |
| Step-16x2 | relaxed | 144 | strip step | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 260 | 0 |
| Step-24x1 | ideal | 108 | strip step | P | R | R | P (cold) | P (cold) | 118 | 0 |
| Step-8x1 | ideal | 36 | strip step | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 39 | 0 |
| Step-8x2 | coll_edgebend0.15 | 72 | strip step | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 92 | 0 |
| Step-8x2 | coll_topspacing-3pct | 72 | strip step | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 92 | 0 |
| Step-8x2 | ideal | 72 | strip step | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 92 | 0 |
| Step-8x2 | path_detach1 | 72 | strip step | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 92 | 0 |
| Step-8x2 | path_detach2 | 72 | strip step | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 92 | 0 |
| Step-8x2 | path_detach3 | 72 | strip step | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 92 | 0 |
| Step-8x2 | pert05 | 72 | strip step | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 92 | 0 |
| Step-8x2 | pert10 | 72 | strip step | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 92 | 0 |
| Step-8x2 | relaxed | 72 | strip step | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 92 | 0 |
| Step-8x2_edge-vacancy_plus_foot-adatom | ideal | 72 | kink / edge rearrangement | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 93 | 0 |
| T-4x4 | coll_strain+1pct | 64 | flat Au(111) | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 50 | 0 |
| T-4x4 | coll_topspacing-3pct | 64 | flat Au(111) | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 50 | 0 |
| T-4x4 | ideal | 64 | flat Au(111) | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 50 | 0 |
| T-4x4 | pert05 | 64 | flat Au(111) | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 50 | 0 |
| T-4x4 | pert10 | 64 | flat Au(111) | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 50 | 0 |
| T-4x4 | relaxed | 64 | flat Au(111) | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 50 | 0 |
| V1 | coll_strain+1pct | 63 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 51 | 0 |
| V1 | coll_topspacing-3pct | 63 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 51 | 0 |
| V1 | ideal | 63 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 51 | 0 |
| V1 | pert05 | 63 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 51 | 0 |
| V1 | pert10 | 63 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 51 | 0 |
| V1 | relaxed | 63 | point defect | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 51 | 0 |
| V2 | ideal | 62 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 44 | 0 |
| V3 | ideal | 61 | point defect | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 51 | 0 |

## Per-family coverage of the new end points (completed / planned)

- composite: 0/4
- flat Au(111): 0/16
- kink / edge rearrangement: 0/26
- point defect: 0/44
- reconstruction-related: 0/16
- single-layer island: 0/24
- single-layer pit: 0/24
- strip step: 0/32
- vicinal step face: 0/28
