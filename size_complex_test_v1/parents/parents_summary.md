# Multi-layer parents for size_complex_test_v1

6 parents on the 32 x 32 x 4 base; every upper level from one 7-layer fcc build; every upper atom supported by its three fcc hollows (unsupported sites removed and counted).

| id | class | description | atoms | L4 | L5 | L6 | top in DFT cell (A) | min Au-Au (A) | status |
|---|---|---|---|---|---|---|---|---|---|
| PM1a_s11 | M1 | tiered island: wide lower tier, smaller upper tier(s) | 4433 | 245 | 92 | 0 | 17.003 | 2.940 | PASS |
| PM1b_s11 | M1 | tiered island: wide lower tier, smaller upper tier(s) | 4483 | 245 | 112 | 30 | 19.404 | 2.940 | PASS |
| PM2a_s11 | M2 | step bunching: two step edges locally merged into a double step | 4876 | 512 | 268 | 0 | 17.003 | 2.940 | PASS |
| PM2b_s11 | M2 | step bunching: two step edges locally merged into a double step | 4801 | 513 | 192 | 0 | 17.003 | 2.940 | PASS |
| PM3a_s11 | M3 | open valley between two multi-level ridges; floor = complete original terrace | 5123 | 643 | 384 | 0 | 17.003 | 2.940 | PASS |
| PM3b_s11 | M3 | open valley between two multi-level ridges; floor = complete original terrace | 5182 | 638 | 384 | 64 | 19.404 | 2.940 | PASS |

The production solvent window requires the highest Au below SOL_Z1 - 15 = 19.603 A (metal bottom at 5.0 A): three-level parents top at 19.40 A before thermal motion, so cells cut from them will need the vertical-box check (step 4 of the README).
