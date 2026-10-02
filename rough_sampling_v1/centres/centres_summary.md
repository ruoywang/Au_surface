# Candidate centres (MD frames)

1000 centres from a pool of 15365 legal exposed atoms sampled on 256 frames (32 parents x [10.0, 20.0, 165.0, 172.0, 179.0, 186.0, 193.0, 199.0]). Surface atoms rejected before sampling: {'detached': 0, 'close_contact': 8}. Rare-flagged (no pool neighbour at cosine > 0.98): 0 in pool, 0 chosen (kept: geometry legal).

| class | train | val | test |
|---|---|---|---|
| A | 188 | 31 | 31 |
| B | 188 | 31 | 31 |
| C | 188 | 31 | 31 |
| D | 188 | 31 | 31 |

| CN stratum | in pool | chosen |
|---|---|---|
| cn<=5 | 2831 | 228 |
| cn6-7 | 1438 | 166 |
| cn8-9 | 9658 | 439 |
| cn>=10 | 1438 | 167 |

Novelty vs the 52 old geometries (1 - max cosine): median 0.0030, 10th-90th pct 0.0009-0.0086.
Context (within 10 A), medians over chosen: cn 8.00, relief_10A 2.82, n_under_10A 14.00, n_levels_10A 2.00
Detached atoms seen (parent@step: count): none

These are candidates for extraction; the final list is fixed only after repair, with the second-layer (whole-cell) stratification.
