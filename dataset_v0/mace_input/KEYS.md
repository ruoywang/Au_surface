# Keys in *.mace.extxyz (FermiMACE fork, env cp-mace-0524)

- `forces` (arrays, eV/A): pass `--forces_key=forces`. Raw forces on all atoms (fixed atoms not zeroed).
- `potential` (info, eV): actual converged mu_e (CPM-ion line). Literal key required by the fork (data/utils.py:145).
  It is a TARGET (the fork predicts mu_e from structure + N_e).
- `electron` (info): actual converged N_e (CPM-ion line). Literal key required (data/utils.py:146). It is a model INPUT.
- energy candidates (info, eV), choose ONE with `--energy_key=<name>`:
  `E_free_TOTEN_eV` (free energy; the forces are its derivatives at fixed N_e), `E_without_entropy_eV`
  (the group's existing convention in surface_charge/6-Au converters), `E_sigma0_eV`, `GCE_code_eV`
  (= E_sigma0 - mu_e*(N_e - N0)), `TOTEN_minus_mu_dN_eV` (= F - mu_e*(N_e - N0)).
  No generic `energy` key is written on purpose: the mapping to the trainer's energy target is UNCONFIRMED.
- `config_type=Default` (for --config_type_weights), `state_id`, `geometry_id`.
