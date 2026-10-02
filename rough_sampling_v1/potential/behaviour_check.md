# Behaviour check of the Au FLARE potential on this project's structures

Bulk fcc: energy minimum at a0 = 4.160 A (project uses 4.158 A); E/atom -3.2169 eV.

| structure | atoms | geometry | max force at DFT geometry (eV/A) | energy drop on minimising (eV) | RMS / max displacement of mobile atoms (A) |
|---|---|---|---|---|---|
| T-4x4 | 64 | relax | 0.057 | -0.012 | 0.021 / 0.024 |
| Step-8x2 | 72 | relax | 0.076 | -0.023 | 0.051 / 0.115 |
| A1-fcc | 65 | relax | 0.092 | -0.020 | 0.025 / 0.068 |
| A1-hcp | 65 | relax | 0.119 | -0.025 | 0.031 / 0.083 |
| Island-7-compact | 151 | relax | 0.178 | -0.054 | 0.034 / 0.225 |
| Pit-7-compact | 137 | relax | 0.077 | -0.019 | 0.018 / 0.040 |

- adatom hcp minus fcc: +0.047 eV
- adatom binding relative to bulk: +0.567 eV
- 8-atom step strip (two 8-atom edges) formation: +0.610 eV, +0.038 eV per edge atom
- 7-atom island formation: +1.886 eV; 7-atom pit formation: +1.795 eV

These are sanity signs for a candidate-generation tool. None of them is compared with the CP-DFT data and none enters a label.
