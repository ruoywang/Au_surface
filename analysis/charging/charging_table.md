# Morphology -> charging response (from the converged CP records)

Generated 2026-09-29 23:34 from `dataset_v1/states.json`. 107 distinct geometries. U = mu0 - mu_e(actual), mu0 = -4.9071 eV (internal reference, not vs RHE, not a per-structure PZC). sigma = -e (N_e - N_e^0) / A_proj with N_e^0 = 11 per Au atom. Every value uses the ACTUAL converged mu_e and N_e, never TARGETMU.

`C_sec` is the secant slope between adjacent sampled points, not a differential capacitance fit. A PZC is given only where the sampled points bracket sigma = 0; otherwise the row says how far the nearest point is.

| family | structure | config | N | A_proj (A^2) | pts | sigma at U (uC/cm^2) | C_sec (uF/cm^2) | U(sigma=0) |
|---|---|---|---|---|---|---|---|---|
| flat Au(111) | Flat-16x1 | ideal | 64 | 119.8 | 5 | -0.509V: -6.52 ; -0.199V: -3.10 ; -0.000V: -0.82 ; +0.191V: +1.34 ; +0.509V: +4.80 | 11.0@-0.35 ; 11.5@-0.10 ; 11.3@+0.10 ; 10.9@+0.35 | +0.0726 V |
| flat Au(111) | Flat-8x2 | ideal | 64 | 119.8 | 5 | -0.509V: -6.50 ; -0.199V: -3.10 ; -0.000V: -0.83 ; +0.192V: +1.34 ; +0.509V: +4.78 | 11.0@-0.35 ; 11.5@-0.10 ; 11.2@+0.10 ; 10.9@+0.35 | +0.0732 V |
| flat Au(111) | T-4x4 | coll_strain+1pct | 64 | 122.2 | 5 | -0.508V: -5.86 ; -0.199V: -2.45 ; -0.000V: -0.17 ; +0.201V: +2.08 ; +0.500V: +5.32 | 11.1@-0.35 ; 11.4@-0.10 ; 11.2@+0.10 ; 10.9@+0.35 | +0.0152 V |
| flat Au(111) | T-4x4 | coll_topspacing-3pct | 64 | 119.8 | 5 | -0.504V: -5.47 ; -0.199V: -2.06 ; -0.000V: +0.21 ; +0.201V: +2.45 ; +0.500V: +5.68 | 11.2@-0.35 ; 11.4@-0.10 ; 11.1@+0.10 ; 10.8@+0.35 | -0.0188 V |
| flat Au(111) | T-4x4 | ideal | 64 | 119.8 | 5 | -0.504V: -5.59 ; -0.199V: -2.23 ; +0.000V: -0.00 ; +0.201V: +2.22 ; +0.500V: +5.40 | 11.0@-0.35 ; 11.2@-0.10 ; 11.0@+0.10 ; 10.7@+0.35 | +0.0000 V |
| flat Au(111) | T-4x4 | pert05 | 64 | 119.8 | 5 | -0.505V: -5.59 ; -0.199V: -2.21 ; -0.004V: -0.00 ; +0.201V: +2.27 ; +0.500V: +5.48 | 11.1@-0.35 ; 11.3@-0.10 ; 11.1@+0.10 ; 10.7@+0.35 | -0.0044 V |
| flat Au(111) | T-4x4 | pert10 | 64 | 119.8 | 5 | -0.506V: -5.85 ; -0.199V: -2.47 ; -0.000V: -0.21 ; +0.201V: +2.01 ; +0.500V: +5.23 | 11.0@-0.35 ; 11.3@-0.10 ; 11.1@+0.10 ; 10.8@+0.35 | +0.0193 V |
| flat Au(111) | T-4x4 | relaxed | 64 | 119.8 | 5 | -0.505V: -5.71 ; -0.199V: -2.34 ; +0.001V: -0.08 ; +0.201V: +2.13 ; +0.500V: +5.33 | 11.0@-0.35 ; 11.3@-0.10 ; 11.0@+0.10 ; 10.7@+0.35 | +0.0081 V |
| point defect | A1-fcc | coll_strain+1pct | 65 | 122.2 | 5 | -0.506V: -4.69 ; -0.209V: -1.31 ; -0.000V: +1.11 ; +0.204V: +3.42 ; +0.499V: +6.64 | 11.4@-0.36 ; 11.6@-0.10 ; 11.3@+0.10 ; 10.9@+0.35 | -0.0960 V |
| point defect | A1-fcc | coll_topspacing-3pct | 65 | 119.8 | 5 | -0.506V: -4.41 ; -0.200V: -0.89 ; -0.008V: +1.34 ; +0.205V: +3.75 ; +0.499V: +6.92 | 11.5@-0.35 ; 11.6@-0.10 ; 11.3@+0.10 ; 10.8@+0.35 | -0.1233 V |
| point defect | A1-fcc | ideal | 65 | 119.8 | 5 | -0.504V: -5.24 ; -0.200V: -1.62 ; -0.000V: +0.78 ; +0.203V: +3.18 ; +0.499V: +6.52 | 11.9@-0.35 ; 12.0@-0.10 ; 11.8@+0.10 ; 11.3@+0.35 | -0.0653 V |
| point defect | A1-fcc | path_bridge | 65 | 119.8 | 5 | -0.506V: -4.53 ; -0.200V: -1.00 ; +0.001V: +1.34 ; +0.205V: +3.65 ; +0.499V: +6.84 | 11.5@-0.35 ; 11.6@-0.10 ; 11.3@+0.10 ; 10.8@+0.35 | -0.1139 V |
| point defect | A1-fcc | pert05 | 65 | 119.8 | 5 | -0.505V: -4.45 ; -0.200V: -0.94 ; -0.004V: +1.34 ; +0.205V: +3.70 ; +0.499V: +6.88 | 11.5@-0.35 ; 11.6@-0.10 ; 11.3@+0.10 ; 10.8@+0.35 | -0.1188 V |
| point defect | A1-fcc | pert10 | 65 | 119.8 | 5 | -0.504V: -4.35 ; -0.200V: -0.85 ; +0.000V: +1.48 ; +0.205V: +3.80 ; +0.499V: +6.98 | 11.5@-0.35 ; 11.6@-0.10 ; 11.3@+0.10 ; 10.8@+0.35 | -0.1272 V |
| point defect | A1-fcc | relaxed | 65 | 119.8 | 5 | -0.505V: -4.51 ; -0.200V: -1.00 ; -0.000V: +1.33 ; +0.205V: +3.65 ; +0.499V: +6.84 | 11.5@-0.35 ; 11.7@-0.10 ; 11.3@+0.10 ; 10.9@+0.35 | -0.1142 V |
| point defect | A1-hcp | coll_strain+1pct | 65 | 122.2 | 5 | -0.506V: -4.73 ; -0.207V: -1.31 ; -0.000V: +1.08 ; +0.204V: +3.40 ; +0.499V: +6.62 | 11.4@-0.36 ; 11.6@-0.10 ; 11.3@+0.10 ; 10.9@+0.35 | -0.0936 V |
| point defect | A1-hcp | coll_topspacing-3pct | 65 | 119.8 | 5 | -0.505V: -4.44 ; -0.200V: -0.93 ; -0.006V: +1.34 ; +0.205V: +3.73 ; +0.499V: +6.90 | 11.5@-0.35 ; 11.7@-0.10 ; 11.3@+0.10 ; 10.8@+0.35 | -0.1204 V |
| point defect | A1-hcp | ideal | 65 | 119.8 | 5 | -0.504V: -5.23 ; -0.200V: -1.61 ; -0.000V: +0.80 ; +0.203V: +3.19 ; +0.499V: +6.53 | 11.9@-0.35 ; 12.0@-0.10 ; 11.8@+0.10 ; 11.3@+0.35 | -0.0665 V |
| point defect | A1-hcp | pert05 | 65 | 119.8 | 5 | -0.505V: -4.43 ; -0.200V: -0.92 ; -0.006V: +1.34 ; +0.205V: +3.73 ; +0.499V: +6.91 | 11.5@-0.35 ; 11.6@-0.10 ; 11.3@+0.10 ; 10.8@+0.35 | -0.1213 V |
| point defect | A1-hcp | pert10 | 65 | 119.8 | 5 | -0.506V: -5.12 ; -0.200V: -1.62 ; +0.002V: +0.71 ; +0.203V: +3.02 ; +0.499V: +6.26 | 11.4@-0.35 ; 11.5@-0.10 ; 11.5@+0.10 ; 10.9@+0.35 | -0.0592 V |
| point defect | A1-hcp | relaxed | 65 | 119.8 | 5 | -0.505V: -4.57 ; -0.200V: -1.04 ; +0.001V: +1.30 ; +0.205V: +3.62 ; +0.499V: +6.82 | 11.6@-0.35 ; 11.7@-0.10 ; 11.4@+0.10 ; 10.9@+0.35 | -0.1107 V |
| point defect | A3 | ideal | 67 | 119.8 | 5 | -0.507V: -5.30 ; -0.192V: -1.34 ; -0.000V: +1.13 ; +0.202V: +3.70 ; +0.499V: +7.32 | 12.6@-0.35 ; 12.9@-0.10 ; 12.7@+0.10 ; 12.2@+0.35 | -0.0878 V |
| point defect | V1 | coll_strain+1pct | 63 | 122.2 | 5 | -0.509V: -5.95 ; -0.200V: -2.46 ; -0.000V: -0.14 ; +0.201V: +2.16 ; +0.500V: +5.48 | 11.3@-0.35 ; 11.7@-0.10 ; 11.5@+0.10 ; 11.1@+0.35 | +0.0121 V |
| point defect | V1 | coll_topspacing-3pct | 63 | 119.8 | 5 | -0.506V: -5.55 ; -0.200V: -2.05 ; -0.000V: +0.28 ; +0.202V: +2.59 ; +0.500V: +5.89 | 11.4@-0.35 ; 11.7@-0.10 ; 11.4@+0.10 ; 11.1@+0.35 | -0.0242 V |
| point defect | V1 | ideal | 63 | 119.8 | 5 | -0.508V: -5.84 ; -0.200V: -2.42 ; -0.000V: -0.15 ; +0.201V: +2.11 ; +0.500V: +5.38 | 11.1@-0.35 ; 11.4@-0.10 ; 11.2@+0.10 ; 10.9@+0.35 | +0.0127 V |
| point defect | V1 | pert05 | 63 | 119.8 | 5 | -0.509V: -5.80 ; -0.200V: -2.32 ; +0.001V: -0.00 ; +0.201V: +2.27 ; +0.500V: +5.56 | 11.2@-0.35 ; 11.6@-0.10 ; 11.4@+0.10 ; 11.0@+0.35 | +0.0012 V |
| point defect | V1 | pert10 | 63 | 119.8 | 5 | -0.509V: -6.07 ; -0.200V: -2.60 ; -0.000V: -0.29 ; +0.201V: +2.00 ; +0.500V: +5.32 | 11.2@-0.35 ; 11.6@-0.10 ; 11.4@+0.10 ; 11.1@+0.35 | +0.0250 V |
| point defect | V1 | relaxed | 63 | 119.8 | 5 | -0.508V: -5.84 ; -0.199V: -2.39 ; -0.000V: -0.09 ; +0.201V: +2.18 ; +0.500V: +5.47 | 11.2@-0.35 ; 11.5@-0.10 ; 11.3@+0.10 ; 11.0@+0.35 | +0.0079 V |
| point defect | V2 | ideal | 62 | 119.8 | 5 | -0.510V: -5.93 ; -0.200V: -2.43 ; +0.008V: -0.00 ; +0.201V: +2.23 ; +0.500V: +5.57 | 11.3@-0.35 ; 11.7@-0.10 ; 11.5@+0.10 ; 11.2@+0.35 | +0.0077 V |
| point defect | V3 | ideal | 61 | 119.8 | 5 | -0.499V: -6.03 ; -0.200V: -2.54 ; -0.000V: -0.12 ; +0.201V: +2.27 ; +0.510V: +5.87 | 11.7@-0.35 ; 12.1@-0.10 ; 11.9@+0.10 ; 11.6@+0.36 | +0.0104 V |
| reconstruction-related | R1-hcp-terminated | ideal | 64 | 119.8 | 5 | -0.503V: -5.52 ; -0.199V: -2.15 ; -0.009V: -0.00 ; +0.201V: +2.32 ; +0.500V: +5.51 | 11.1@-0.35 ; 11.3@-0.10 ; 11.0@+0.10 ; 10.7@+0.35 | -0.0093 V |
| reconstruction-related | R1-hcp-terminated | pert05 | 64 | 119.8 | 5 | -0.504V: -5.57 ; -0.199V: -2.18 ; -0.006V: -0.00 ; +0.201V: +2.29 ; +0.500V: +5.49 | 11.1@-0.35 ; 11.3@-0.10 ; 11.1@+0.10 ; 10.7@+0.35 | -0.0064 V |
| reconstruction-related | R1-hcp-terminated | pert10 | 64 | 119.8 | 5 | -0.505V: -5.80 ; -0.199V: -2.43 ; -0.000V: -0.18 ; +0.201V: +2.04 ; +0.500V: +5.24 | 11.0@-0.35 ; 11.3@-0.10 ; 11.0@+0.10 ; 10.7@+0.35 | +0.0162 V |
| reconstruction-related | R1-hcp-terminated | relaxed | 64 | 119.8 | 5 | -0.504V: -5.68 ; -0.199V: -2.30 ; +0.000V: -0.05 ; +0.201V: +2.17 ; +0.500V: +5.37 | 11.1@-0.35 ; 11.3@-0.10 ; 11.0@+0.10 ; 10.7@+0.35 | +0.0045 V |
| reconstruction-related | R2-stripe-wall | ideal | 65 | 119.8 | 5 | -0.502V: -6.28 ; -0.198V: -2.94 ; -0.000V: -0.69 ; +0.200V: +1.53 ; +0.509V: +4.85 | 11.0@-0.35 ; 11.3@-0.10 ; 11.1@+0.10 ; 10.7@+0.35 | +0.0622 V |
| reconstruction-related | R2-stripe-wall | pert05 | 65 | 119.8 | 5 | -0.508V: -6.85 ; -0.198V: -3.53 ; -0.003V: -1.34 ; +0.200V: +0.90 ; +0.507V: +4.17 | 10.7@-0.35 ; 11.2@-0.10 ; 11.0@+0.10 ; 10.7@+0.35 | +0.1187 V |
| reconstruction-related | R2-stripe-wall | pert10 | 65 | 119.8 | 5 | -0.510V: -6.94 ; -0.198V: -3.62 ; +0.005V: -1.34 ; +0.201V: +0.81 ; +0.506V: +4.09 | 10.7@-0.35 ; 11.2@-0.10 ; 11.0@+0.10 ; 10.7@+0.35 | +0.1271 V |
| reconstruction-related | R2-stripe-wall | relaxed | 65 | 119.8 | 5 | -0.507V: -6.78 ; -0.198V: -3.45 ; -0.000V: -1.22 ; +0.200V: +0.98 ; +0.507V: +4.26 | 10.8@-0.35 ; 11.2@-0.10 ; 11.0@+0.10 ; 10.7@+0.35 | +0.1109 V |
| strip step | Step-16x1 | ideal | 72 | 119.8 | 5 | -0.499V: -6.31 ; -0.200V: -2.88 ; -0.000V: -0.52 ; +0.200V: +1.84 ; +0.500V: +5.24 | 11.4@-0.35 ; 11.8@-0.10 ; 11.7@+0.10 ; 11.4@+0.35 | +0.0440 V |
| strip step | Step-16x2 | ideal | 144 | 239.6 | 3 | -0.199V: -2.88 ; -0.000V: -0.50 ; +0.201V: +1.86 | 11.9@-0.10 ; 11.7@+0.10 | +0.0427 V |
| strip step | Step-16x2 | pert05 | 144 | 239.6 | 2 | -0.199V: -2.91 ; +0.201V: +1.75 | 11.7@+0.00 | +0.0508 V |
| strip step | Step-16x2 | pert10 | 144 | 239.6 | 2 | -0.199V: -3.03 ; +0.201V: +1.65 | 11.7@+0.00 | +0.0603 V |
| strip step | Step-16x2 | relaxed | 144 | 239.6 | 3 | -0.199V: -3.00 ; -0.000V: -0.66 ; +0.201V: +1.66 | 11.7@-0.10 ; 11.6@+0.10 | +0.0572 V |
| strip step | Step-24x1 | ideal | 108 | 179.7 | 5 | -0.499V: -6.51 ; -0.199V: -3.11 ; -0.000V: -0.76 ; +0.201V: +1.57 ; +0.500V: +4.92 | 11.3@-0.35 ; 11.8@-0.10 ; 11.6@+0.10 ; 11.2@+0.35 | +0.0657 V |
| strip step | Step-8x1 | ideal | 36 | 59.9 | 5 | -0.499V: -6.35 ; -0.200V: -2.83 ; -0.000V: -0.39 ; +0.200V: +2.06 ; +0.508V: +5.72 | 11.8@-0.35 ; 12.2@-0.10 ; 12.3@+0.10 ; 11.9@+0.35 | +0.0319 V |
| strip step | Step-8x2 | coll_edgebend0.15 | 72 | 119.8 | 5 | -0.499V: -6.31 ; -0.201V: -2.86 ; -0.000V: -0.44 ; +0.201V: +1.97 ; +0.500V: +5.43 | 11.6@-0.35 ; 12.1@-0.10 ; 12.0@+0.10 ; 11.6@+0.35 | +0.0365 V |
| strip step | Step-8x2 | coll_topspacing-3pct | 72 | 119.8 | 5 | -0.499V: -6.05 ; -0.200V: -2.57 ; +0.000V: -0.15 ; +0.201V: +2.26 ; +0.500V: +5.71 | 11.6@-0.35 ; 12.1@-0.10 ; 12.0@+0.10 ; 11.5@+0.35 | +0.0124 V |
| strip step | Step-8x2 | ideal | 72 | 119.8 | 5 | -0.499V: -6.37 ; -0.201V: -2.86 ; -0.000V: -0.39 ; +0.201V: +2.06 ; +0.500V: +5.61 | 11.8@-0.35 ; 12.3@-0.10 ; 12.2@+0.10 ; 11.8@+0.35 | +0.0319 V |
| strip step | Step-8x2 | path_detach1 | 72 | 119.8 | 5 | -0.499V: -6.24 ; -0.201V: -2.75 ; -0.000V: -0.30 ; +0.201V: +2.13 ; +0.500V: +5.63 | 11.7@-0.35 ; 12.2@-0.10 ; 12.1@+0.10 ; 11.7@+0.35 | +0.0249 V |
| strip step | Step-8x2 | path_detach2 | 72 | 119.8 | 5 | -0.499V: -6.07 ; -0.200V: -2.54 ; +0.006V: -0.00 ; +0.201V: +2.39 ; +0.500V: +5.92 | 11.8@-0.35 ; 12.3@-0.10 ; 12.2@+0.10 ; 11.8@+0.35 | +0.0056 V |
| strip step | Step-8x2 | path_detach3 | 72 | 119.8 | 5 | -0.499V: -5.82 ; -0.200V: -2.26 ; +0.000V: +0.23 ; +0.201V: +2.71 ; +0.499V: +6.25 | 11.9@-0.35 ; 12.4@-0.10 ; 12.3@+0.10 ; 11.9@+0.35 | -0.0187 V |
| strip step | Step-8x2 | pert05 | 72 | 119.8 | 5 | -0.499V: -6.20 ; -0.200V: -2.73 ; -0.000V: -0.30 ; +0.201V: +2.11 ; +0.500V: +5.57 | 11.6@-0.35 ; 12.1@-0.10 ; 12.0@+0.10 ; 11.6@+0.35 | +0.0252 V |
| strip step | Step-8x2 | pert10 | 72 | 119.8 | 5 | -0.499V: -6.55 ; -0.201V: -3.11 ; -0.000V: -0.68 ; +0.200V: +1.74 ; +0.500V: +5.23 | 11.6@-0.35 ; 12.1@-0.10 ; 12.1@+0.10 ; 11.7@+0.35 | +0.0564 V |
| strip step | Step-8x2 | relaxed | 72 | 119.8 | 5 | -0.499V: -6.32 ; -0.201V: -2.87 ; +0.001V: -0.45 ; +0.200V: +1.95 ; +0.500V: +5.41 | 11.5@-0.35 ; 12.1@-0.10 ; 12.0@+0.10 ; 11.6@+0.35 | +0.0377 V |
| vicinal step face | Au211 | coll_strain+1pct | 48 | 86.4 | 5 | -0.499V: -5.75 ; -0.201V: -2.51 ; -0.000V: -0.24 ; +0.200V: +2.02 ; +0.508V: +5.41 | 10.9@-0.35 ; 11.3@-0.10 ; 11.3@+0.10 ; 11.0@+0.35 | +0.0215 V |
| vicinal step face | Au211 | coll_topspacing-3pct | 48 | 84.7 | 5 | -0.500V: -5.42 ; -0.200V: -2.17 ; -0.007V: -0.00 ; +0.200V: +2.32 ; +0.500V: +5.56 | 10.9@-0.35 ; 11.2@-0.10 ; 11.2@+0.10 ; 10.8@+0.35 | -0.0069 V |
| vicinal step face | Au211 | ideal | 48 | 84.7 | 5 | -0.507V: -5.44 ; -0.194V: -1.89 ; -0.000V: +0.35 ; +0.201V: +2.67 ; +0.500V: +6.00 | 11.3@-0.35 ; 11.6@-0.10 ; 11.5@+0.10 ; 11.1@+0.35 | -0.0306 V |
| vicinal step face | Au211 | pert05 | 48 | 84.7 | 5 | -0.509V: -5.47 ; -0.200V: -2.09 ; +0.000V: +0.16 ; +0.200V: +2.40 ; +0.500V: +5.64 | 10.9@-0.35 ; 11.3@-0.10 ; 11.2@+0.10 ; 10.8@+0.35 | -0.0138 V |
| vicinal step face | Au211 | pert10 | 48 | 84.7 | 5 | -0.500V: -5.60 ; -0.200V: -2.28 ; -0.002V: -0.00 ; +0.200V: +2.31 ; +0.500V: +5.63 | 11.1@-0.35 ; 11.5@-0.10 ; 11.4@+0.10 ; 11.1@+0.35 | -0.0019 V |
| vicinal step face | Au211 | relaxed | 48 | 84.7 | 5 | -0.500V: -5.49 ; -0.200V: -2.22 ; +0.000V: +0.03 ; +0.200V: +2.28 ; +0.500V: +5.54 | 10.9@-0.35 ; 11.3@-0.10 ; 11.2@+0.10 ; 10.9@+0.35 | -0.0030 V |
| vicinal step face | Au221 | coll_strain+1pct | 28 | 52.9 | 5 | -0.507V: -5.56 ; -0.199V: -2.03 ; -0.000V: +0.32 ; +0.200V: +2.66 ; +0.500V: +6.03 | 11.5@-0.35 ; 11.8@-0.10 ; 11.7@+0.10 ; 11.2@+0.35 | -0.0278 V |
| vicinal step face | Au221 | coll_topspacing-3pct | 28 | 51.9 | 5 | -0.506V: -5.32 ; -0.199V: -1.79 ; -0.001V: +0.55 ; +0.200V: +2.88 ; +0.500V: +6.23 | 11.5@-0.35 ; 11.8@-0.10 ; 11.6@+0.10 ; 11.1@+0.35 | -0.0474 V |
| vicinal step face | Au221 | ideal | 28 | 51.9 | 5 | -0.510V: -5.94 ; -0.199V: -2.28 ; -0.000V: +0.15 ; +0.199V: +2.56 ; +0.500V: +6.05 | 11.8@-0.35 ; 12.2@-0.10 ; 12.1@+0.10 ; 11.6@+0.35 | -0.0125 V |
| vicinal step face | Au221 | pert05 | 28 | 51.9 | 5 | -0.506V: -5.37 ; -0.199V: -1.84 ; -0.000V: +0.51 ; +0.200V: +2.84 ; +0.500V: +6.19 | 11.5@-0.35 ; 11.8@-0.10 ; 11.6@+0.10 ; 11.2@+0.35 | -0.0437 V |
| vicinal step face | Au221 | pert10 | 28 | 51.9 | 5 | -0.500V: -6.19 ; -0.200V: -2.77 ; -0.000V: -0.40 ; +0.199V: +1.95 ; +0.507V: +5.47 | 11.4@-0.35 ; 11.9@-0.10 ; 11.8@+0.10 ; 11.4@+0.35 | +0.0334 V |
| vicinal step face | Au221 | relaxed | 28 | 51.9 | 5 | -0.508V: -5.59 ; -0.199V: -2.04 ; +0.000V: +0.33 ; +0.199V: +2.66 ; +0.500V: +6.04 | 11.5@-0.35 ; 11.9@-0.10 ; 11.7@+0.10 ; 11.3@+0.35 | -0.0274 V |
| vicinal step face | Au332 | ideal | 21 | 40.5 | 5 | -0.509V: -6.21 ; -0.199V: -2.60 ; +0.000V: -0.20 ; +0.198V: +2.17 ; +0.505V: +5.71 | 11.6@-0.35 ; 12.1@-0.10 ; 11.9@+0.10 ; 11.5@+0.35 | +0.0168 V |
| vicinal step face | Au554 | ideal | 36 | 70.2 | 5 | -0.499V: -6.47 ; -0.200V: -3.13 ; -0.000V: -0.77 ; +0.200V: +1.56 ; +0.507V: +5.04 | 11.2@-0.35 ; 11.8@-0.10 ; 11.7@+0.10 ; 11.3@+0.35 | +0.0659 V |
| kink / edge rearrangement | Kink-edge1 | ideal | 109 | 179.7 | 5 | -0.499V: -6.25 ; -0.201V: -2.72 ; -0.000V: -0.25 ; +0.201V: +2.22 ; +0.500V: +5.77 | 11.8@-0.35 ; 12.3@-0.10 ; 12.3@+0.10 ; 11.9@+0.35 | +0.0199 V |
| kink / edge rearrangement | Kink-edge1 | path_kinkmove1 | 109 | 179.7 | 4 | -0.499V: -6.35 ; -0.201V: -2.91 ; +0.201V: +1.93 ; +0.500V: +5.39 | 11.5@-0.35 ; 12.0@+0.00 ; 11.6@+0.35 | +0.0406 V |
| kink / edge rearrangement | Kink-edge1 | path_kinkmove2 | 109 | 179.7 | 4 | -0.499V: -6.08 ; -0.201V: -2.64 ; +0.201V: +2.20 ; +0.499V: +5.65 | 11.5@-0.35 ; 12.0@+0.00 ; 11.6@+0.35 | +0.0184 V |
| kink / edge rearrangement | Kink-edge1 | pert05 | 109 | 179.7 | 4 | -0.499V: -6.28 ; -0.200V: -2.83 ; +0.201V: +2.00 ; +0.500V: +5.46 | 11.6@-0.35 ; 12.0@+0.00 ; 11.6@+0.35 | +0.0347 V |
| kink / edge rearrangement | Kink-edge1 | pert10 | 109 | 179.7 | 4 | -0.499V: -6.29 ; -0.201V: -2.85 ; +0.201V: +1.96 ; +0.500V: +5.41 | 11.5@-0.35 ; 12.0@-0.00 ; 11.5@+0.35 | +0.0374 V |
| kink / edge rearrangement | Kink-edge1 | relaxed | 109 | 179.7 | 5 | -0.499V: -6.34 ; -0.201V: -2.91 ; +0.000V: -0.49 ; +0.201V: +1.92 ; +0.500V: +5.38 | 11.5@-0.35 ; 12.1@-0.10 ; 12.0@+0.10 ; 11.6@+0.35 | +0.0409 V |
| kink / edge rearrangement | Kink-edge2 | ideal | 109 | 179.7 | 5 | -0.499V: -6.16 ; -0.201V: -2.61 ; +0.010V: -0.00 ; +0.201V: +2.36 ; +0.499V: +5.90 | 11.9@-0.35 ; 12.4@-0.10 ; 12.3@+0.11 ; 11.9@+0.35 | +0.0096 V |
| kink / edge rearrangement | Kink-edge2 | path_kinkmove1 | 109 | 179.7 | 4 | -0.499V: -5.89 ; -0.201V: -2.43 ; +0.201V: +2.41 ; +0.499V: +5.85 | 11.6@-0.35 ; 12.0@+0.00 ; 11.5@+0.35 | +0.0010 V |
| kink / edge rearrangement | Kink-edge2 | path_kinkmove2 | 109 | 179.7 | 4 | -0.499V: -5.98 ; -0.201V: -2.50 ; +0.201V: +2.36 ; +0.499V: +5.82 | 11.6@-0.35 ; 12.1@+0.00 ; 11.6@+0.35 | +0.0065 V |
| kink / edge rearrangement | Kink-edge2 | pert05 | 109 | 179.7 | 2 | -0.200V: -2.77 ; +0.201V: +2.08 | 12.1@+0.00 | +0.0285 V |
| kink / edge rearrangement | Kink-edge2 | pert10 | 109 | 179.7 | 2 | -0.201V: -2.91 ; +0.201V: +1.90 | 12.0@-0.00 | +0.0420 V |
| kink / edge rearrangement | Kink-edge2 | relaxed | 109 | 179.7 | 3 | -0.201V: -2.75 ; +0.000V: -0.32 ; +0.201V: +2.10 | 12.1@-0.10 ; 12.0@+0.10 | +0.0265 V |
| kink / edge rearrangement | Step-8x2_edge-vacancy_plus_foot-adatom | ideal | 72 | 119.8 | 5 | -0.499V: -5.81 ; -0.200V: -2.21 ; -0.000V: +0.29 ; +0.201V: +2.79 ; +0.499V: +6.37 | 12.0@-0.35 ; 12.5@-0.10 ; 12.4@+0.10 ; 12.0@+0.35 | -0.0236 V |
| single-layer island | Island-19-8x8 | ideal | 275 | 479.1 | 2 | -0.007V: -0.23 ; +0.193V: +2.24 | 12.3@+0.09 | +0.0117 V |
| single-layer island | Island-7-8x8 | ideal | 263 | 479.1 | 1 | +0.008V: -0.33 |  | not bracketed (nearest +0.008 V, -0.33) |
| single-layer island | Island-7-compact | ideal | 151 | 269.5 | 3 | -0.196V: -1.59 ; -0.001V: +0.85 ; +0.193V: +3.23 | 12.5@-0.10 ; 12.3@+0.10 | -0.0688 V |
| single-layer island | Island-7-compact | path_detach1 | 151 | 269.5 | 2 | -0.198V: -1.91 ; +0.201V: +2.86 | 12.0@+0.00 | -0.0384 V |
| single-layer island | Island-7-compact | path_detach2 | 151 | 269.5 | 2 | -0.197V: -1.72 ; +0.199V: +3.06 | 12.1@+0.00 | -0.0547 V |
| single-layer island | Island-7-compact | pert05 | 151 | 269.5 | 2 | -0.199V: -2.05 ; +0.203V: +2.71 | 11.8@+0.00 | -0.0262 V |
| single-layer island | Island-7-compact | pert10 | 151 | 269.5 | 2 | -0.199V: -1.94 ; +0.203V: +2.82 | 11.8@+0.00 | -0.0352 V |
| single-layer island | Island-7-compact | relaxed | 151 | 269.5 | 3 | -0.199V: -1.98 ; -0.000V: +0.40 ; +0.201V: +2.75 | 12.0@-0.10 ; 11.7@+0.10 | -0.0333 V |
| single-layer island | Island-7-elongated | ideal | 151 | 269.5 | 3 | -0.194V: -1.77 ; -0.006V: +0.59 ; +0.192V: +3.05 | 12.6@-0.10 ; 12.4@+0.09 | -0.0537 V |
| single-layer island | Island-7-elongated | pert05 | 151 | 269.5 | 2 | -0.197V: -1.79 ; +0.199V: +2.99 | 12.1@+0.00 | -0.0492 V |
| single-layer island | Island-7-elongated | pert10 | 151 | 269.5 | 2 | -0.196V: -1.90 ; +0.197V: +2.90 | 12.2@+0.00 | -0.0406 V |
| single-layer island | Island-7-elongated | relaxed | 151 | 269.5 | 3 | -0.197V: -1.88 ; -0.000V: +0.52 ; +0.199V: +2.90 | 12.2@-0.10 ; 11.9@+0.10 | -0.0428 V |
| single-layer pit | Pit-19-8x8 | ideal | 237 | 479.1 | 2 | +0.006V: -0.62 ; +0.196V: +1.67 | 12.1@+0.10 | +0.0576 V |
| single-layer pit | Pit-7-8x8 | ideal | 249 | 479.1 | 1 | +0.006V: -1.04 |  | not bracketed (nearest +0.006 V, -1.04) |
| single-layer pit | Pit-7-compact | ideal | 137 | 269.5 | 3 | -0.198V: -2.38 ; -0.000V: -0.00 ; +0.201V: +2.39 | 12.1@-0.10 ; 11.9@+0.10 | -0.0004 V |
| single-layer pit | Pit-7-compact | path_rimin1 | 137 | 269.5 | 2 | -0.198V: -2.21 ; +0.201V: +2.56 | 12.0@+0.00 | -0.0134 V |
| single-layer pit | Pit-7-compact | path_rimin2 | 137 | 269.5 | 2 | -0.197V: -1.98 ; +0.199V: +2.81 | 12.1@+0.00 | -0.0332 V |
| single-layer pit | Pit-7-compact | pert05 | 137 | 269.5 | 2 | -0.198V: -2.20 ; +0.202V: +2.56 | 11.9@+0.00 | -0.0133 V |
| single-layer pit | Pit-7-compact | pert10 | 137 | 269.5 | 2 | -0.198V: -2.80 ; +0.201V: +1.96 | 11.9@+0.00 | +0.0365 V |
| single-layer pit | Pit-7-compact | relaxed | 137 | 269.5 | 3 | -0.198V: -2.18 ; -0.000V: +0.20 ; +0.202V: +2.58 | 12.0@-0.10 ; 11.8@+0.10 | -0.0171 V |
| single-layer pit | Pit-7-trench | ideal | 137 | 269.5 | 3 | -0.196V: -2.52 ; -0.001V: -0.14 ; +0.199V: +2.26 | 12.2@-0.10 ; 12.0@+0.10 | +0.0108 V |
| single-layer pit | Pit-7-trench | pert05 | 137 | 269.5 | 2 | -0.197V: -2.36 ; +0.200V: +2.41 | 12.0@+0.00 | -0.0007 V |
| single-layer pit | Pit-7-trench | pert10 | 137 | 269.5 | 2 | -0.199V: -2.46 ; +0.202V: +2.30 | 11.9@+0.00 | +0.0083 V |
| single-layer pit | Pit-7-trench | relaxed | 137 | 269.5 | 3 | -0.198V: -2.30 ; +0.000V: +0.09 ; +0.200V: +2.47 | 12.1@-0.10 ; 11.9@+0.10 | -0.0071 V |
| composite | C1-island-near-step | ideal | 151 | 239.6 | 2 | -0.000V: +0.39 ; +0.202V: +2.89 | 12.3@+0.10 | not bracketed (nearest -0.000 V, +0.39) |
| composite | C2-island+pit | ideal | 256 | 479.1 | 2 | -0.005V: -0.26 ; +0.193V: +2.18 | 12.3@+0.09 | +0.0165 V |

- PZC bracketed by the sampled range: 104 of 107 geometries; not bracketed: 3.

## Within one cell: morphology only

Cross-cell PZC comparisons carry a numerical offset: audit K.12 records that two *flat* cells differ by 72 meV in neutral mu_e purely from cell shape, k-mesh and PREC. Inside one cell all of that is common, so the columns below are morphology. `ref` is the flat member of that cell where one exists.

| A_proj (A^2) | n | structures | reference | secant C med (range) | U_pzc med (range, mV) | spread (mV) |
|---|---|---|---|---|---|---|
| 119.8 | 28 | A1-fcc, A1-hcp, A3, R1-hcp-terminated, T-4x4, V1, V2, V3 | T-4x4__ideal__8fa1fcc9706c | 11.31 (10.66-12.86) | -14 (-127..+25) | 152 |
| 269.5 | 20 | Island-7-compact, Island-7-elongated, Pit-7-compact, Pit-7-t | - | 12.01 (11.71-12.60) | -30 (-69..+37) | 105 |
| 179.7 | 12 | Kink-edge1, Kink-edge2 | - | 11.89 (11.51-12.39) | +28 (+1..+42) | 41 |
| 119.8 | 11 | Flat-8x2, Step-8x2, Step-8x2_edge-vacancy_plus_foot-adatom | Flat-8x2__ideal__8eddbc7a2602 | 11.95 (10.86-12.50) | +25 (-24..+73) | 97 |
| 119.8 | 6 | Flat-16x1, R2-stripe-wall, Step-16x1 | Flat-16x1__ideal__234bab346496 | 11.04 (10.68-11.83) | +92 (+44..+127) | 83 |
| 51.9 | 5 | Au221 | - | 11.64 (11.15-12.21) | -27 (-47..+33) | 81 |
| 84.7 | 5 | Au211 | - | 11.18 (10.84-11.61) | -7 (-31..-2) | 29 |
| 479.1 | 5 | C2-island+pit, Island-19-8x8, Island-7-8x8, Pit-19-8x8, Pit- | - | 12.32 (12.08-12.32) | +17 (+12..+58) | 46 |
| 122.2 | 4 | A1-fcc, A1-hcp, T-4x4, V1 | - | 11.33 (10.86-11.65) | -41 (-96..+15) | 111 |
| 239.6 | 4 | Step-16x2 | - | 11.73 (11.56-11.94) | +54 (+43..+60) | 18 |
| 52.9 | 1 | Au221 | - | 11.58 (11.23-11.82) | -28 (-28..-28) | - |
| 86.4 | 1 | Au211 | - | 11.15 (10.85-11.30) | +21 (+21..+21) | - |
| 40.5 | 1 | Au332 | - | 11.77 (11.52-12.08) | +17 (+17..+17) | - |
| 70.2 | 1 | Au554 | - | 11.49 (11.20-11.78) | +66 (+66..+66) | - |
| 59.9 | 1 | Step-8x1 | - | 12.02 (11.75-12.29) | +32 (+32..+32) | - |
| 179.7 | 1 | Step-24x1 | - | 11.46 (11.22-11.80) | +66 (+66..+66) | - |
| 239.6 | 1 | C1-island-near-step | - | 12.34 (12.34-12.34) | - | - |

## Which matters more at a given potential: the PZC shift or the capacitance?

Writing sigma(U) = C (U - U_pzc), the spread of sigma across geometries at a fixed U splits into `<C> x spread(U_pzc)` and `spread(C) x |U - <U_pzc>|`.

| set | n | spread U_pzc (mV) | spread C (%) | sigma spread from PZC | from C | ratio |
|---|---|---|---|---|---|---|
| all geometries | 104 | 254 | 15.3 | 2.95 uC/cm2 | 0.35 uC/cm2 | 8.5x |
| within cell A=120 A^2 (A1-fcc, A1-hcp, A3...) | 28 | 152 | 14.2 | 1.74 uC/cm2 | 0.35 uC/cm2 | 5.0x |
| within cell A=270 A^2 (Island-7-compact, Island-7-elongated, Pit-7-compact...) | 20 | 105 | 5.5 | 1.26 uC/cm2 | 0.15 uC/cm2 | 8.4x |
| within cell A=180 A^2 (Kink-edge1, Kink-edge2) | 12 | 41 | 4.6 | 0.48 uC/cm2 | 0.09 uC/cm2 | 5.1x |
| within cell A=120 A^2 (Flat-8x2, Step-8x2, Step-8x2_edge-vacancy_plus_foot-adatom) | 11 | 97 | 9.5 | 1.15 uC/cm2 | 0.20 uC/cm2 | 5.8x |
| within cell A=120 A^2 (Flat-16x1, R2-stripe-wall, Step-16x1) | 6 | 83 | 6.7 | 0.91 uC/cm2 | 0.08 uC/cm2 | 11.5x |
| within cell A=52 A^2 (Au221) | 5 | 81 | 3.1 | 0.94 uC/cm2 | 0.08 uC/cm2 | 11.5x |
| within cell A=85 A^2 (Au211) | 5 | 29 | 3.6 | 0.32 uC/cm2 | 0.08 uC/cm2 | 3.8x |
| within cell A=479 A^2 (C2-island+pit, Island-19-8x8, Island-7-8x8...) | 3 | 46 | 1.9 | 0.57 uC/cm2 | 0.04 uC/cm2 | 12.9x |
| within cell A=122 A^2 (A1-fcc, A1-hcp, T-4x4...) | 4 | 111 | 2.1 | 1.26 uC/cm2 | 0.06 uC/cm2 | 21.5x |
| within cell A=240 A^2 (Step-16x2) | 4 | 18 | 1.6 | 0.21 uC/cm2 | 0.03 uC/cm2 | 7.5x |

(evaluated at U = +0.2 V; the ratio grows as U approaches the median PZC and shrinks far from it.)

## Potential-induced relative stabilisation, same composition and cell (fixed U window $\pm$0.2 V)

d(Omega_A - Omega_B) = -int [N_A(mu) - N_B(mu)] dmu over the SAME window for every pair, so the values are comparable. Negative means A is relatively stabilised as U becomes more positive. This is the potential-INDUCED change only: it does not rank stability at the reference potential, and it is never taken across different Au counts (that would need a reservoir term). N_A - N_B is nearly constant over this window, so the value tracks the PZC offset, shown alongside.

242 pairs cover the window; 96 are between DIFFERENT structures. Those, largest first:

| A | B | N | d(Omega_A-Omega_B) (eV) | mean N_A-N_B (e) | dU_pzc (mV) |
|---|---|---|---|---|---|
| Step-8x2__pert10__a5a115f445a9 | Step-8x2_edge-vacancy_plus_foot-adatom__ideal__cc6ba4ff8a63 | 72 | -0.0291 | +0.0727 | +80.0 |
| A1-fcc__pert10__10fb3f2d2638 | A1-hcp__pert10__6c834986bfe0 | 65 | +0.0233 | -0.0583 | -68.0 |
| Step-8x2__relaxed__b3d8731e142a | Step-8x2_edge-vacancy_plus_foot-adatom__ideal__cc6ba4ff8a63 | 72 | -0.0223 | +0.0558 | +61.3 |
| Step-8x2__coll_edgebend0.15__85b0ce5989f7 | Step-8x2_edge-vacancy_plus_foot-adatom__ideal__cc6ba4ff8a63 | 72 | -0.0219 | +0.0548 | +60.0 |
| A1-fcc__coll_topspacing-3pct__cff7df59fe39 | A1-hcp__pert10__6c834986bfe0 | 65 | +0.0219 | -0.0547 | -64.0 |
| Kink-edge1__relaxed__e391de3e71c6 | Kink-edge2__path_kinkmove1__8dca1d700e9e | 109 | -0.0216 | +0.0540 | +39.9 |
| Kink-edge1__path_kinkmove1__4c97421862c4 | Kink-edge2__path_kinkmove1__8dca1d700e9e | 109 | -0.0214 | +0.0534 | +39.6 |
| Step-8x2__ideal__25c159395cb3 | Step-8x2_edge-vacancy_plus_foot-adatom__ideal__cc6ba4ff8a63 | 72 | -0.0205 | +0.0512 | +55.5 |
| A1-fcc__pert05__138462b661be | A1-hcp__pert10__6c834986bfe0 | 65 | +0.0204 | -0.0511 | -59.6 |
| A1-fcc__pert10__10fb3f2d2638 | A1-hcp__ideal__f0af6d51e549 | 65 | +0.0202 | -0.0505 | -60.7 |
| Kink-edge1__pert10__5faa30511248 | Kink-edge2__path_kinkmove1__8dca1d700e9e | 109 | -0.0196 | +0.0490 | +36.4 |
| A1-fcc__relaxed__e7fadb33a876 | A1-hcp__pert10__6c834986bfe0 | 65 | +0.0189 | -0.0473 | -55.0 |
| A1-fcc__path_bridge__1857661e9126 | A1-hcp__pert10__6c834986bfe0 | 65 | +0.0188 | -0.0469 | -54.7 |
| A1-fcc__coll_topspacing-3pct__cff7df59fe39 | A1-hcp__ideal__f0af6d51e549 | 65 | +0.0187 | -0.0469 | -56.8 |
| Kink-edge1__relaxed__e391de3e71c6 | Kink-edge2__path_kinkmove2__81e14678d5f0 | 109 | -0.0186 | +0.0466 | +34.4 |
| A1-fcc__ideal__dae654c15784 | A1-hcp__pert05__349e8880c105 | 65 | -0.0186 | +0.0465 | +56.1 |
| Kink-edge1__path_kinkmove1__4c97421862c4 | Kink-edge2__path_kinkmove2__81e14678d5f0 | 109 | -0.0184 | +0.0460 | +34.1 |
| A1-fcc__ideal__dae654c15784 | A1-hcp__coll_topspacing-3pct__6df671924f54 | 65 | -0.0183 | +0.0458 | +55.1 |
| Kink-edge1__pert05__d3baabf3cfd9 | Kink-edge2__path_kinkmove1__8dca1d700e9e | 109 | -0.0182 | +0.0455 | +33.7 |
| Step-8x2__pert05__9aa197e38f03 | Step-8x2_edge-vacancy_plus_foot-adatom__ideal__cc6ba4ff8a63 | 72 | -0.0179 | +0.0448 | +48.8 |
| Step-8x2__path_detach1__c6d4a8e29e02 | Step-8x2_edge-vacancy_plus_foot-adatom__ideal__cc6ba4ff8a63 | 72 | -0.0178 | +0.0446 | +48.5 |
| A1-fcc__pert05__138462b661be | A1-hcp__ideal__f0af6d51e549 | 65 | +0.0173 | -0.0433 | -52.4 |
| Kink-edge1__pert10__5faa30511248 | Kink-edge2__path_kinkmove2__81e14678d5f0 | 109 | -0.0166 | +0.0415 | +31.0 |
| Kink-edge1__relaxed__e391de3e71c6 | Kink-edge2__ideal__6ffdcdbfea14 | 109 | -0.0166 | +0.0415 | +31.3 |
| Kink-edge1__path_kinkmove1__4c97421862c4 | Kink-edge2__ideal__6ffdcdbfea14 | 109 | -0.0164 | +0.0409 | +31.0 |
| A1-fcc__relaxed__e7fadb33a876 | A1-hcp__ideal__f0af6d51e549 | 65 | +0.0158 | -0.0395 | -47.8 |
| A1-fcc__path_bridge__1857661e9126 | A1-hcp__ideal__f0af6d51e549 | 65 | +0.0156 | -0.0391 | -47.4 |
| Kink-edge1__pert05__d3baabf3cfd9 | Kink-edge2__path_kinkmove2__81e14678d5f0 | 109 | -0.0152 | +0.0381 | +28.3 |
| A1-fcc__ideal__dae654c15784 | A1-hcp__relaxed__723e5d5fa1ac | 65 | -0.0151 | +0.0377 | +45.4 |
| Kink-edge1__pert10__5faa30511248 | Kink-edge2__ideal__6ffdcdbfea14 | 109 | -0.0146 | +0.0365 | +27.8 |
| Step-8x2__coll_topspacing-3pct__b08ca5bfbfcf | Step-8x2_edge-vacancy_plus_foot-adatom__ideal__cc6ba4ff8a63 | 72 | -0.0133 | +0.0333 | +35.9 |
| Kink-edge1__pert05__d3baabf3cfd9 | Kink-edge2__ideal__6ffdcdbfea14 | 109 | -0.0132 | +0.0330 | +25.1 |
| Kink-edge1__path_kinkmove2__a8157b517224 | Kink-edge2__pert10__678a8f2ea128 | 109 | +0.0126 | -0.0316 | -23.6 |
| R1-hcp-terminated__pert10__9e9d78f478ca | T-4x4__coll_topspacing-3pct__d9010df4b5dc | 64 | -0.0117 | +0.0293 | +35.0 |
| Kink-edge1__ideal__0a4b4cbb19ec | Kink-edge2__pert10__678a8f2ea128 | 109 | +0.0115 | -0.0287 | -22.1 |
| Step-8x2__path_detach2__dd1e5ef5486e | Step-8x2_edge-vacancy_plus_foot-adatom__ideal__cc6ba4ff8a63 | 72 | -0.0109 | +0.0272 | +29.2 |
| Kink-edge1__ideal__0a4b4cbb19ec | Kink-edge2__path_kinkmove1__8dca1d700e9e | 109 | -0.0106 | +0.0264 | +18.9 |
| R1-hcp-terminated__ideal__7618cb82f683 | T-4x4__pert10__59bc520c9795 | 64 | +0.0094 | -0.0236 | -28.6 |
| Kink-edge1__path_kinkmove2__a8157b517224 | Kink-edge2__path_kinkmove1__8dca1d700e9e | 109 | -0.0094 | +0.0235 | +17.4 |
| R1-hcp-terminated__pert05__ccebe5455555 | T-4x4__pert10__59bc520c9795 | 64 | +0.0085 | -0.0212 | -25.7 |

- Same-structure pairs (ideal / relaxed / perturbed / collective / path images of one parent, 146 pairs): |dOmega| median 0.0066 eV, max 0.0272 eV. That is the scale on which sampling configurations of ONE morphology already differ, i.e. the bar a cross-morphology difference must clear to be meaningful.

