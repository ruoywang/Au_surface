# Rough states: CANDIDATE list of 200 (seed 20261002, frozen = False)

Nothing here is submitted or frozen. Approval = `finalize_200.py freeze`, then `rough_dft.py submit --first 8 --confirm`.

## Where they come from

| class | train | val | test |
|---|---|---|---|
| A | 40 | 5 | 5 |
| B | 40 | 5 | 5 |
| C | 40 | 5 | 5 |
| D | 40 | 5 | 5 |

Parents used: 32 of 32; cells per parent min/median/max 2/5/10 (caps: {"A/train": {"parents": 6, "cap_per_parent": 10}, "A/val": {"parents": 1, "cap_per_parent": 5}, "A/test": {"parents": 1, "cap_per_parent": 3}, "B/train": {"parents": 6, "cap_per_parent": 10}, "B/val": {"parents": 1, "cap_per_parent": 5}, "B/test": {"parents": 1, "cap_per_parent": 3}, "C/train": {"parents": 6, "cap_per_parent": 10}, "C/val": {"parents": 1, "cap_per_parent": 5}, "C/test": {"parents": 1, "cap_per_parent": 3}, "D/train": {"parents": 6, "cap_per_parent": 10}, "D/val": {"parents": 1, "cap_per_parent": 5}, "D/test": {"parents": 1, "cap_per_parent": 3}}); cells per (parent, frame) max 2 (cap 2).
Frames (ps): 10 ps: 34, 20 ps: 37, 165 ps: 29, 172 ps: 26, 179 ps: 21, 186 ps: 15, 193 ps: 23, 199 ps: 15

| CN stratum of the centre | states |
|---|---|
| cn<=5 | 59 |
| cn6-7 | 45 |
| cn8-9 | 57 |
| cn>=10 | 39 |

## Size, atoms, k-mesh

Cell sizes: {'8x8': 177, '10x8': 23}
Atoms: 201: 1, 202: 1, 203: 1, 205: 1, 206: 1, 209: 2, 211: 1, 212: 1, 213: 1, 216: 2, 217: 1, 223: 1, 226: 1, 228: 1, 229: 1, 233: 1, 234: 2, 235: 1, 238: 2, 239: 2, 240: 1, 241: 1, 242: 3, 243: 5, 244: 3, 245: 2, 247: 1, 248: 4, 254: 1, 256: 1, 258: 1, 260: 1, 261: 1, 263: 2, 264: 2, 265: 2, 266: 1, 267: 4, 268: 3, 270: 1, 271: 6, 272: 4, 273: 1, 274: 4, 275: 2, 276: 1, 277: 4, 278: 5, 279: 2, 280: 5, 281: 4, 282: 2, 283: 2, 284: 3, 285: 3, 286: 1, 287: 4, 288: 2, 289: 5, 290: 3, 291: 2, 292: 3, 294: 2, 295: 1, 297: 5, 300: 2, 301: 3, 302: 5, 303: 2, 304: 4, 305: 2, 306: 3, 307: 4, 308: 2, 309: 1, 311: 4, 312: 1, 315: 2, 316: 2, 317: 2, 318: 1, 321: 1, 327: 1, 328: 1, 340: 1, 343: 1, 344: 1, 345: 2, 348: 1, 349: 1, 350: 1, 352: 1, 355: 1, 363: 1, 366: 2, 373: 1, 375: 1, 377: 1
k-mesh: {'2x2x1': 200}
Size pairs (test, same U): pair1 A PA2_s51_f010000_a4064, pair2 B PB2_s51_f096500_a4040, pair3 C PC2_s51_f086000_a2769, pair4 D PD2_s51_f005000_a4565
Repair: {'no repair needed': 53, 'seam band (128 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (110 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 4, 'seam band (60 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (95 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 5, 'seam band (93 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (77 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (133 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (106 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (115 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 5, 'seam band (86 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (124 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 4, 'seam band (101 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (114 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 6, 'seam band (116 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (143 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (122 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (107 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (103 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (117 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (94 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (87 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (97 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (123 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 4, 'seam band (113 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (135 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (127 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (89 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 5, 'seam band (88 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (92 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 4, 'seam band (104 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (84 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (91 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (83 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (145 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (102 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (98 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 5, 'seam band (120 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (105 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (108 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 4, 'seam band (144 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (90 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (118 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (111 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (85 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (73 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (71 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (74 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (80 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (72 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (68 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (70 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 3, 'seam band (76 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (62 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (61 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (69 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (48 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (50 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (56 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (75 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (51 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (78 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (99 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (100 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (129 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (96 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 2, 'seam band (121 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (139 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (52 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1, 'seam band (64 atoms) minimised; core, bottom layers and the rest of the periphery fixed; occupancy unchanged': 1}

## Morphology coverage of the chosen cells (medians and ranges)

| descriptor | min | median | max |
|---|---|---|---|
| f_cn<=5 | 0.000 | 0.026 | 0.123 |
| f_cn6-7 | 0.000 | 0.173 | 0.312 |
| f_cn8-9 | 0.274 | 0.576 | 1.000 |
| f_cn>=10 | 0.000 | 0.221 | 0.427 |
| relief_layers | 0.235 | 1.179 | 2.273 |
| n_levels | 1.000 | 2.000 | 3.000 |
| adatom_level_ML | 0.000 | 0.344 | 1.000 |
| missing_terrace_ML | -0.050 | 0.016 | 0.859 |
| exposed_pit_ML | 0.000 | 0.000 | 0.984 |
| n_islands | 0.000 | 1.000 | 3.000 |
| max_island_frac | 0.000 | 0.344 | 1.000 |
| n_pits | 0.000 | 0.000 | 4.000 |
| step_density | 0.074 | 0.319 | 0.588 |

## Potentials

| U bin (V) | states |
|---|---|
| -0.5 to -0.4 | 19 |
| -0.4 to -0.3 | 20 |
| -0.3 to -0.2 | 20 |
| -0.2 to -0.1 | 21 |
| -0.1 to +0.0 | 20 |
| +0.0 to +0.1 | 20 |
| +0.1 to +0.2 | 19 |
| +0.2 to +0.3 | 20 |
| +0.3 to +0.4 | 21 |
| +0.4 to +0.5 | 20 |

## Rejections

Centres with no passing cell: 63 of 1000. Reasons over all failed trials:
- seam created environments absent from the parent: (layer, cn) = [(4, 4)]: 121
- seam created environments absent from the parent: (layer, cn) = [(3, 6)]: 52
- seam created environments absent from the parent: (layer, cn) = [(3, 5)]: 32
- seam created environments absent from the parent: (layer, cn) = [(3, 4)]: 28
- seam created environments absent from the parent: (layer, cn) = [(4, 3)]: 27
- seam created environments absent from the parent: (layer, cn) = [(3, 5), (3, 6)]: 17
- seam created environments absent from the parent: (layer, cn) = [(3, 4), (3, 5)]: 13
- seam created environments absent from the parent: (layer, cn) = [(3, 7), (4, 3)]: 10
- contact 2.44 A < 2.5 created by the seam (parent atoms 3561, 3786): 10
- contact 2.04 A < 2.5 created by the seam (parent atoms 3563, 4160): 10

Geometry duplicates dropped (occupancy patterns within 2 sites under translation, same size): 295; e.g. PA2_s37_f082500_a3984_8x8 = PD1_s51_f089500_a3744_8x8; PA1_s51_f082500_a3840_8x8 = PD1_s51_f089500_a3744_8x8; PC2_s37_f096500_a2298_8x8 = PC2_s51_f086000_a2769_8x8; PC2_s23_f099500_a2323_8x8 = PC1_s51_f005000_a2677_8x8; PC1_s11_f093000_a3527_8x8 = PC1_s51_f093000_a2504_8x8
Seam-affected atoms (CN reduced by >= 2 vs the parent) per chosen cell: median 5, max 16 (all passing cells: median 5, max 19; cells above the cap of 16 excluded: 9).

Novelty of the chosen cells' centres vs the 52 old geometries (1 - max cosine): median 0.0040; rare-flagged centres chosen: 0.
