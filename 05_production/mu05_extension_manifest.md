# mu05 extension: +-0.5 V end points for every distinct geometry

Generated 2026-09-29 22:59 by `scripts/extend_mu05.py --report`. mu_target = mu0 - dU, mu0 = -4.9071 eV (internal reference). TARGETMU -5.4071 (dU = +0.5 V) and -4.4071 (dU = -0.5 V); existing -5.1071/-4.9071/-4.7071 kept; no intermediate points added. Static CP single points (IBRION=-1, NSW=0) on the geometry actually computed (relaxed = accepted CONTCAR). Warm start = ICHARG=1 from a copy of the same-side neighbour's CHGCAR.

Distinct geometries: 107. Extension tasks: 214 -> status {'running': 7, 'pending': 59, 'completed': 148}. Estimated 243 node-h (measured same-structure single-point means); actual so far 109 node-h over 148 completed.

Status legend: C completed/accepted, R running, P pending, U reused, F failed, - not planned. `w:` warm start source.

| structure | config | N | family | -5.1071 | -4.9071 | -4.7071 | -5.4071 (+0.5 V) | -4.4071 (-0.5 V) | est min | actual min |
|---|---|---|---|---|---|---|---|---|---|---|
| A1-fcc | coll_strain+1pct | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4063 | C (w:-4.7071) mu_e=-4.4012 | 52 | 74 |
| A1-fcc | coll_topspacing-3pct | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4058 | C (w:-4.7071) mu_e=-4.4011 | 52 | 63 |
| A1-fcc | ideal | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4061 | C (w:-4.7071) mu_e=-4.4036 | 52 | 64 |
| A1-fcc | path_bridge | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4059 | C (w:-4.7071) mu_e=-4.4009 | 52 | 63 |
| A1-fcc | pert05 | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4058 | C (w:-4.7071) mu_e=-4.4022 | 52 | 63 |
| A1-fcc | pert10 | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4058 | C (w:-4.7071) mu_e=-4.4027 | 52 | 71 |
| A1-fcc | relaxed | 65 | point defect | C | U | C | C (w:-5.1071) mu_e=-5.4059 | C (w:-4.7071) mu_e=-4.4022 | 52 | 66 |
| A1-hcp | coll_strain+1pct | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4062 | C (w:-4.7071) mu_e=-4.4014 | 54 | 75 |
| A1-hcp | coll_topspacing-3pct | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4059 | C (w:-4.7071) mu_e=-4.4022 | 54 | 58 |
| A1-hcp | ideal | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4061 | C (w:-4.7071) mu_e=-4.4034 | 54 | 66 |
| A1-hcp | pert05 | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4057 | C (w:-4.7071) mu_e=-4.4023 | 54 | 76 |
| A1-hcp | pert10 | 65 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4063 | C (w:-4.7071) mu_e=-4.4010 | 54 | 64 |
| A1-hcp | relaxed | 65 | point defect | C | U | C | C (w:-5.1071) mu_e=-5.4059 | C (w:-4.7071) mu_e=-4.4021 | 54 | 66 |
| A3 | ideal | 67 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4062 | C (w:-4.7071) mu_e=-4.4001 | 52 | 60 |
| Au211 | coll_strain+1pct | 48 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4149 | C (w:-4.7071) mu_e=-4.4079 | 31 | 43 |
| Au211 | coll_topspacing-3pct | 48 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4074 | 31 | 41 |
| Au211 | ideal | 48 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4002 | 31 | 37 |
| Au211 | pert05 | 48 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.3982 | 31 | 35 |
| Au211 | pert10 | 48 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4069 | C (w:-4.7071) mu_e=-4.4075 | 31 | 39 |
| Au211 | relaxed | 48 | vicinal step face | C | U | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4074 | 31 | 39 |
| Au221 | coll_strain+1pct | 28 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.3998 | 22 | 23 |
| Au221 | coll_topspacing-3pct | 28 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4009 | 22 | 20 |
| Au221 | ideal | 28 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4069 | C (w:-4.7071) mu_e=-4.3974 | 22 | 20 |
| Au221 | pert05 | 28 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4009 | 22 | 20 |
| Au221 | pert10 | 28 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4144 | C (w:-4.7071) mu_e=-4.4075 | 22 | 22 |
| Au221 | relaxed | 28 | vicinal step face | C | U | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.3996 | 22 | 21 |
| Au332 | ideal | 21 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4125 | C (w:-4.7071) mu_e=-4.3977 | 17 | 16 |
| Au554 | ideal | 36 | vicinal step face | C | C | C | C (w:-5.1071) mu_e=-5.4139 | C (w:-4.7071) mu_e=-4.4081 | 40 | 38 |
| C1-island-near-step | ideal | 151 | composite | C | C | - | P (w:-5.1071) | P (w:-4.9071) | 195 | 0 |
| C2-island+pit | ideal | 256 | composite | C | C | - | R (w:-4.9071) | R (w:-4.9071) | 430 | 0 |
| Flat-16x1 | ideal | 64 | flat Au(111) | U | U | C | C (w:-5.1071) mu_e=-5.4165 | C (w:-4.7071) mu_e=-4.3977 | 81 | 82 |
| Flat-8x2 | ideal | 64 | flat Au(111) | C | C | C | C (w:-5.1071) mu_e=-5.4165 | C (w:-4.7071) mu_e=-4.3979 | 77 | 76 |
| Island-19-8x8 | ideal | 275 | single-layer island | C | C | - | P (w:-4.9071) | P (w:-4.9071) | 479 | 0 |
| Island-7-8x8 | ideal | 263 | single-layer island | - | C | - | P (w:-4.9071) | R (w:-4.9071) | 448 | 0 |
| Island-7-compact | ideal | 151 | single-layer island | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-compact | path_detach1 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-compact | path_detach2 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-compact | pert05 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-compact | pert10 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-compact | relaxed | 151 | single-layer island | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 286 | 0 |
| Island-7-elongated | ideal | 151 | single-layer island | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 295 | 0 |
| Island-7-elongated | pert05 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 295 | 0 |
| Island-7-elongated | pert10 | 151 | single-layer island | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 295 | 0 |
| Island-7-elongated | relaxed | 151 | single-layer island | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 295 | 0 |
| Kink-edge1 | ideal | 109 | kink / edge rearrangement | C | C | C | C (w:-5.1071) mu_e=-5.4066 | C (w:-4.7071) mu_e=-4.4080 | 177 | 231 |
| Kink-edge1 | path_kinkmove1 | 109 | kink / edge rearrangement | C | - | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4081 | 177 | 234 |
| Kink-edge1 | path_kinkmove2 | 109 | kink / edge rearrangement | C | - | C | C (w:-5.1071) mu_e=-5.4066 | C (w:-4.7071) mu_e=-4.4080 | 177 | 233 |
| Kink-edge1 | pert05 | 109 | kink / edge rearrangement | C | - | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4080 | 177 | 223 |
| Kink-edge1 | pert10 | 109 | kink / edge rearrangement | C | - | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4080 | 177 | 236 |
| Kink-edge1 | relaxed | 109 | kink / edge rearrangement | C | U | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4080 | 177 | 254 |
| Kink-edge2 | ideal | 109 | kink / edge rearrangement | C | C | C | C (w:-5.1071) mu_e=-5.4066 | C (w:-4.7071) mu_e=-4.4079 | 173 | 246 |
| Kink-edge2 | path_kinkmove1 | 109 | kink / edge rearrangement | C | - | C | C (w:-5.1071) mu_e=-5.4065 | C (w:-4.7071) mu_e=-4.4078 | 173 | 228 |
| Kink-edge2 | path_kinkmove2 | 109 | kink / edge rearrangement | C | - | C | C (w:-5.1071) mu_e=-5.4065 | C (w:-4.7071) mu_e=-4.4079 | 173 | 230 |
| Kink-edge2 | pert05 | 109 | kink / edge rearrangement | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 173 | 0 |
| Kink-edge2 | pert10 | 109 | kink / edge rearrangement | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 173 | 0 |
| Kink-edge2 | relaxed | 109 | kink / edge rearrangement | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 173 | 0 |
| Pit-19-8x8 | ideal | 237 | single-layer pit | C | C | - | R (w:-4.9071) | R (w:-4.9071) | 383 | 0 |
| Pit-7-8x8 | ideal | 249 | single-layer pit | - | C | - | R (cold) | R (cold) | 413 | 0 |
| Pit-7-compact | ideal | 137 | single-layer pit | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-compact | path_rimin1 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-compact | path_rimin2 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-compact | pert05 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-compact | pert10 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-compact | relaxed | 137 | single-layer pit | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 242 | 0 |
| Pit-7-trench | ideal | 137 | single-layer pit | C | C | C | P (w:-4.9071) | P (w:-4.7071) | 268 | 0 |
| Pit-7-trench | pert05 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 268 | 0 |
| Pit-7-trench | pert10 | 137 | single-layer pit | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 268 | 0 |
| Pit-7-trench | relaxed | 137 | single-layer pit | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 268 | 0 |
| R1-hcp-terminated | ideal | 64 | reconstruction-related | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4039 | 46 | 54 |
| R1-hcp-terminated | pert05 | 64 | reconstruction-related | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4032 | 46 | 60 |
| R1-hcp-terminated | pert10 | 64 | reconstruction-related | C | C | C | C (w:-5.1071) mu_e=-5.4069 | C (w:-4.7071) mu_e=-4.4022 | 46 | 66 |
| R1-hcp-terminated | relaxed | 64 | reconstruction-related | C | U | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4031 | 46 | 60 |
| R2-stripe-wall | ideal | 65 | reconstruction-related | C | C | C | C (w:-5.1071) mu_e=-5.4166 | C (w:-4.7071) mu_e=-4.4050 | 81 | 86 |
| R2-stripe-wall | pert05 | 65 | reconstruction-related | C | C | C | C (w:-5.1071) mu_e=-5.4138 | C (w:-4.7071) mu_e=-4.3995 | 81 | 70 |
| R2-stripe-wall | pert10 | 65 | reconstruction-related | C | C | C | C (w:-5.1071) mu_e=-5.4135 | C (w:-4.7071) mu_e=-4.3976 | 81 | 74 |
| R2-stripe-wall | relaxed | 65 | reconstruction-related | C | U | C | C (w:-5.1071) mu_e=-5.4138 | C (w:-4.7071) mu_e=-4.3998 | 81 | 71 |
| Step-16x1 | ideal | 72 | strip step | U | U | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4077 | 97 | 90 |
| Step-16x2 | ideal | 144 | strip step | C | C | C | P (w:-5.1071) | P (w:-4.7071) | 260 | 0 |
| Step-16x2 | pert05 | 144 | strip step | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 260 | 0 |
| Step-16x2 | pert10 | 144 | strip step | C | - | C | P (w:-5.1071) | P (w:-4.7071) | 260 | 0 |
| Step-16x2 | relaxed | 144 | strip step | C | U | C | P (w:-5.1071) | P (w:-4.7071) | 260 | 0 |
| Step-24x1 | ideal | 108 | strip step | C | C | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4078 | 118 | 300 |
| Step-8x1 | ideal | 36 | strip step | C | U | C | C (w:-5.1071) mu_e=-5.4150 | C (w:-4.7071) mu_e=-4.4078 | 39 | 37 |
| Step-8x2 | coll_edgebend0.15 | 72 | strip step | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4079 | 92 | 112 |
| Step-8x2 | coll_topspacing-3pct | 72 | strip step | C | C | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4077 | 92 | 103 |
| Step-8x2 | ideal | 72 | strip step | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4079 | 92 | 108 |
| Step-8x2 | path_detach1 | 72 | strip step | C | C | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4079 | 92 | 109 |
| Step-8x2 | path_detach2 | 72 | strip step | C | C | C | C (w:-5.1071) mu_e=-5.4066 | C (w:-4.7071) mu_e=-4.4078 | 92 | 115 |
| Step-8x2 | path_detach3 | 72 | strip step | C | C | C | C (w:-5.1071) mu_e=-5.4065 | C (w:-4.7071) mu_e=-4.4077 | 92 | 124 |
| Step-8x2 | pert05 | 72 | strip step | C | C | C | C (w:-5.1071) mu_e=-5.4067 | C (w:-4.7071) mu_e=-4.4078 | 92 | 114 |
| Step-8x2 | pert10 | 72 | strip step | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4080 | 92 | 119 |
| Step-8x2 | relaxed | 72 | strip step | C | U | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4079 | 92 | 115 |
| Step-8x2_edge-vacancy_plus_foot-adatom | ideal | 72 | kink / edge rearrangement | C | C | C | C (w:-5.1071) mu_e=-5.4066 | C (w:-4.7071) mu_e=-4.4076 | 93 | 110 |
| T-4x4 | coll_strain+1pct | 64 | flat Au(111) | C | C | C | C (w:-5.1071) mu_e=-5.4069 | C (w:-4.7071) mu_e=-4.3995 | 50 | 64 |
| T-4x4 | coll_topspacing-3pct | 64 | flat Au(111) | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4028 | 50 | 55 |
| T-4x4 | ideal | 64 | flat Au(111) | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4033 | 50 | 52 |
| T-4x4 | pert05 | 64 | flat Au(111) | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4023 | 50 | 53 |
| T-4x4 | pert10 | 64 | flat Au(111) | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4011 | 50 | 58 |
| T-4x4 | relaxed | 64 | flat Au(111) | C | U | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.4022 | 50 | 68 |
| V1 | coll_strain+1pct | 63 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.3979 | 51 | 70 |
| V1 | coll_topspacing-3pct | 63 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4066 | C (w:-4.7071) mu_e=-4.4010 | 51 | 56 |
| V1 | ideal | 63 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.3991 | 51 | 57 |
| V1 | pert05 | 63 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.3979 | 51 | 56 |
| V1 | pert10 | 63 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.3981 | 51 | 52 |
| V1 | relaxed | 63 | point defect | C | U | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.3987 | 51 | 56 |
| V2 | ideal | 62 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4068 | C (w:-4.7071) mu_e=-4.3973 | 44 | 58 |
| V3 | ideal | 61 | point defect | C | C | C | C (w:-5.1071) mu_e=-5.4170 | C (w:-4.7071) mu_e=-4.4077 | 51 | 57 |

## Per-family coverage of the new end points (completed / planned)

- composite: 0/4
- flat Au(111): 16/16
- kink / edge rearrangement: 20/26
- point defect: 44/44
- reconstruction-related: 16/16
- single-layer island: 0/24
- single-layer pit: 0/24
- strip step: 24/32
- vicinal step face: 28/28
