# Replacement proposal (200 states of the frozen list) -- a PROPOSAL for one-time confirmation

| action | states |
|---|---|
| resize_to_6x6 | 146 |
| keep_8x8 | 33 |
| release_for_complex | 10 |
| keep_10x8_no_smaller | 9 |
| kept_reference_8x8 | 2 |

Atoms summed over the 200 states: 55966 (frozen) -> 37490 (proposed); with the cost model ~N^1.5 the single-point cost sum scales by 0.57x (an extrapolation; the pairs will measure it).

Potentials: every resized state keeps its U; complex candidates inherit the released U bin. Nothing is re-drawn.

