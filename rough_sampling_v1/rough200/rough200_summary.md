# Rough states: CANDIDATE list of 200 (seed 20261002, frozen = False)

Nothing here is submitted or frozen. Approval = `finalize_200.py freeze`, then `rough_dft.py submit --first 8 --confirm`.

Occupancy fingerprint: registered on the frozen bottom layer and the measured stacking step; invariance (permutation, periodic re-representation, translations) verified on all 1043 passing cells, 0 failures.

## Changes against the previous candidate list

Kept from `rough200/rough200_manifest.v1_37a9299.json`: 162 states with their potentials unchanged. Replaced: 38 out, 38 in (new cells took the bins vacated by the replaced ones).

| out (previous) | in (new) |
|---|---|
| PA2_s51_f005000_a4132_8x8 | PA1_s37_f082500_a4544_8x8 |
| PA2_s51_f082500_a3576_10x8 | PA2_s23_f082500_a3799_8x8 |
| PC2_s51_f082500_a3776_8x8 | PA1_s23_f005000_a4199_8x8 |
| PA1_s11_f089500_a4346_8x8 | PA1_s37_f082500_a4154_8x8 |
| PA1_s11_f096500_a3497_8x8 | PA1_s11_f082500_a3340_8x8 |
| PA1_s23_f005000_a3455_8x8 | PA1_s51_f096500_a4100_8x8 |
| PA1_s23_f010000_a4304_8x8 | PA2_s23_f096500_a4385_8x8 |
| PA1_s23_f086000_a3270_8x8 | PA1_s37_f099500_a3366_8x8 |
| PA1_s23_f086000_a4535_8x8 | PA1_s11_f096500_a3733_8x8 |
| PA1_s37_f099500_a4124_8x8 | PA1_s23_f089500_a3835_8x8 |
| PA2_s11_f010000_a4422_8x8 | PA2_s37_f089500_a3499_10x8 |
| PA2_s23_f082500_a3256_8x8 | PA2_s51_f005000_a4132_10x8 |
| PA2_s23_f096500_a4582_8x8 | PA2_s51_f005000_a3520_8x8 |
| PB1_s11_f005000_a4143_8x8 | PB2_s11_f093000_a4444_8x8 |
| PB1_s11_f010000_a4143_10x8 | PB1_s51_f010000_a4156_8x8 |
| PB1_s23_f010000_a4171_8x8 | PB1_s11_f096500_a4201_8x8 |
| PB1_s23_f086000_a4071_8x8 | PB1_s51_f093000_a4156_8x8 |
| PB1_s37_f005000_a4096_8x8 | PB2_s23_f005000_a4364_8x8 |
| PB1_s51_f089500_a4218_8x8 | PB1_s37_f096500_a3696_8x8 |
| PB1_s51_f093000_a4042_8x8 | PB1_s37_f089500_a3920_8x8 |
| PB1_s51_f093000_a4158_8x8 | PB1_s51_f093000_a3395_8x8 |
| PB2_s11_f086000_a3632_10x8 | PB2_s23_f089500_a3887_8x8 |
| PB2_s11_f093000_a3658_8x8 | PB1_s51_f089500_a4164_8x8 |
| PB2_s23_f010000_a4354_8x8 | PB1_s37_f005000_a3375_8x8 |
| PB2_s23_f089500_a4158_8x8 | PB2_s11_f005000_a3202_10x8 |
| PB2_s23_f096500_a3504_8x8 | PB1_s37_f082500_a2643_10x8 |
| PC1_s11_f082500_a3094_8x8 | PB2_s37_f096500_a4030_8x8 |
| PC1_s11_f089500_a2132_8x8 | PC1_s37_f010000_a3458_8x8 |
| PC1_s23_f093000_a3581_8x8 | PC2_s23_f082500_a2482_8x8 |
| PC1_s37_f005000_a3288_8x8 | PC2_s23_f010000_a3567_8x8 |
| PC1_s37_f010000_a2367_8x8 | PC1_s11_f082500_a3556_8x8 |
| PC1_s51_f082500_a2746_8x8 | PC1_s51_f099500_a2757_8x8 |
| PC2_s11_f005000_a2078_8x8 | PC2_s23_f096500_a3536_8x8 |
| PC2_s23_f005000_a3537_10x8 | PC2_s11_f010000_a2078_8x8 |
| PC2_s23_f010000_a3615_8x8 | PC2_s11_f086000_a3780_8x8 |
| PC2_s23_f082500_a3593_8x8 | PC1_s51_f099500_a2714_8x8 |
| PA2_s37_f089500_a3376_8x8 | PC1_s51_f093000_a2316_8x8 |
| PB2_s37_f005000_a4420_8x8 | PC2_s51_f005000_a3646_8x8 |

Coverage before -> after: class x split unchanged; CN strata {'cn8-9': 57, 'cn<=5': 59, 'cn>=10': 39, 'cn6-7': 45} -> {'cn>=10': 41, 'cn8-9': 58, 'cn<=5': 58, 'cn6-7': 43}; cells per size {'8x8': 177, '10x8': 23} -> {'8x8': 177, '10x8': 23}.

## Where they come from

| class | train | val | test |
|---|---|---|---|
| A | 40 | 5 | 5 |
| B | 40 | 5 | 5 |
| C | 40 | 5 | 5 |
| D | 40 | 5 | 5 |

Parents used: 32 of 32; cells per parent min/median/max 2/5/10 (caps: {"A/train": {"parents": 6, "cap_per_parent": 10}, "A/val": {"parents": 1, "cap_per_parent": 5}, "A/test": {"parents": 1, "cap_per_parent": 3}, "B/train": {"parents": 6, "cap_per_parent": 10}, "B/val": {"parents": 1, "cap_per_parent": 5}, "B/test": {"parents": 1, "cap_per_parent": 3}, "C/train": {"parents": 6, "cap_per_parent": 10}, "C/val": {"parents": 1, "cap_per_parent": 5}, "C/test": {"parents": 1, "cap_per_parent": 3}, "D/train": {"parents": 6, "cap_per_parent": 10}, "D/val": {"parents": 1, "cap_per_parent": 5}, "D/test": {"parents": 1, "cap_per_parent": 3}}); cells per (parent, frame) max 2 (cap 2).
Frames (ps): 10 ps: 33, 20 ps: 34, 165 ps: 30, 172 ps: 23, 179 ps: 21, 186 ps: 15, 193 ps: 27, 199 ps: 17

| CN stratum of the centre | states |
|---|---|
| cn<=5 | 58 |
| cn6-7 | 43 |
| cn8-9 | 58 |
| cn>=10 | 41 |

## Size, atoms, k-mesh

Cell sizes: {'8x8': 177, '10x8': 23}
Atoms: 202: 1, 203: 1, 205: 1, 206: 2, 209: 2, 211: 2, 212: 2, 214: 1, 216: 2, 217: 1, 223: 1, 226: 1, 228: 1, 230: 1, 231: 1, 234: 2, 238: 2, 239: 2, 240: 1, 242: 1, 243: 5, 244: 3, 245: 2, 246: 1, 247: 1, 248: 5, 251: 1, 254: 1, 255: 1, 256: 1, 258: 1, 259: 1, 260: 1, 261: 2, 264: 2, 265: 2, 267: 6, 268: 3, 269: 1, 270: 1, 271: 4, 272: 4, 273: 1, 274: 3, 275: 2, 276: 1, 277: 5, 278: 5, 279: 1, 280: 6, 281: 4, 282: 1, 283: 2, 284: 3, 285: 3, 286: 1, 287: 5, 288: 1, 289: 4, 290: 4, 291: 2, 292: 4, 293: 1, 294: 1, 295: 1, 297: 2, 300: 2, 301: 3, 302: 5, 303: 4, 304: 3, 305: 3, 306: 4, 307: 4, 308: 2, 309: 1, 311: 3, 312: 1, 315: 2, 316: 1, 317: 1, 318: 1, 321: 1, 327: 1, 328: 1, 336: 1, 338: 1, 339: 1, 340: 1, 343: 1, 344: 1, 345: 1, 348: 1, 349: 1, 350: 1, 352: 1, 355: 1, 363: 1, 366: 1, 373: 1, 375: 1, 383: 1
k-mesh: {'2x2x1': 200}
Size pairs (test, same U): pair1 A PA2_s51_f010000_a4064, pair2 B PB2_s51_f096500_a4040, pair3 C PC2_s51_f086000_a2769, pair4 D PD2_s51_f005000_a4565
Repair: {'no repair': 52, 'seam band minimised': 148}

## Morphology coverage of the chosen cells (medians and ranges)

| descriptor | min | median | max |
|---|---|---|---|
| f_cn<=5 | 0.000 | 0.026 | 0.123 |
| f_cn6-7 | 0.000 | 0.173 | 0.312 |
| f_cn8-9 | 0.293 | 0.580 | 1.000 |
| f_cn>=10 | 0.000 | 0.217 | 0.427 |
| relief_layers | 0.235 | 1.182 | 2.273 |
| n_levels | 1.000 | 2.000 | 3.000 |
| adatom_level_ML | 0.000 | 0.344 | 1.000 |
| terrace_vacancy_ML | 0.000 | 0.031 | 1.000 |
| exposed_pit_ML | 0.000 | 0.000 | 0.984 |
| n_islands | 0.000 | 1.000 | 3.000 |
| max_island_frac | 0.000 | 0.328 | 1.000 |
| n_pits | 0.000 | 0.000 | 4.000 |
| step_density | 0.074 | 0.318 | 0.588 |
| n_unregistered | 0.000 | 0.000 | 36.000 |
| layer3_deficit | -0.062 | 0.016 | 0.844 |

`terrace_vacancy_ML` counts truly unoccupied layer-3 sites (>= 0); `layer3_deficit` is the count-based 1 - N3/sites and goes negative when an atom caught between levels is counted in layer 3; `n_unregistered` counts atoms farther than 0.3 of a site spacing from every site of their layer.

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

Occupancy duplicates compressed (layer-3/4 site occupancy within 2 sites under translation, same size; similar morphology, NOT a claim of identical geometry): 274; e.g. PA2_s37_f089500_a3376_8x8 ~ PD1_s51_f089500_a3744_8x8; PA2_s37_f082500_a3984_8x8 ~ PD1_s51_f089500_a3744_8x8; PA1_s51_f082500_a3840_8x8 ~ PD1_s51_f089500_a3744_8x8; PC2_s37_f096500_a2298_8x8 ~ PC2_s51_f086000_a2769_8x8; PC2_s23_f099500_a2323_8x8 ~ PC1_s51_f005000_a2677_8x8
Seam-affected atoms (CN reduced by >= 2 vs the parent) per chosen cell: median 5, max 16 (all passing cells: median 5, max 19; cells above the cap of 16 excluded: 9).

Novelty of the chosen cells' centres vs the 52 old geometries (1 - max cosine): median 0.0042; rare-flagged centres chosen: 0.
