# Size-pair inputs (nothing submitted)

| pair | member | source | atoms | k | U (V) | TARGETMU (eV) | mu source | config | FERMICONVERGE | walltime | status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| S1 | 6x6 | new | 163 | 2x2x1 | -0.42 | -4.4891 | actual mu_e of the reference PA1_s51_f08 | production | 0.001 | 13:00:00 | pending (not submitted) |
| S1 | 8x8 | reference | - | - | -0.42 | -4.4871 | its own | production | 0.01 | - | finished |
| S2 | 6x6 | new | 172 | 2x2x1 | +0.46 | -5.3671 | PROVISIONAL: reference not finished; reg | production | 0.001 | 14:00:00 | provisional |
| S2 | 8x8 | reference | - | - | +0.46 | -5.3671 | its own | production | 0.01 | - | running/queued |
| M2 | 6x6 | new | 185 | 2x2x1 | +0.40 | -5.3071 | planned TARGETMU | production | 0.001 | 15:00:00 | pending (not submitted) |
| M2 | 8x8 | new | 314 | 2x2x1 | +0.40 | -5.3071 | planned TARGETMU | production | 0.001 | 33:00:00 | pending (not submitted) |
| M3 | 6x6 | new | 199 | 2x2x1 | -0.40 | -4.5071 | planned TARGETMU | production | 0.001 | 17:00:00 | pending (not submitted) |
| M3 | 8x8 | new | 345 | 2x2x1 | -0.40 | -4.5071 | planned TARGETMU | production | 0.001 | 38:00:00 | pending (not submitted) |

6 new single points; node-hour estimate (rough_dft cost model, seeded median / cold max): 32 / 106 node-h.
All at FERMICONVERGE 0.001; the production references stay at 0.01 (their actual mu_e is the pair's target).
