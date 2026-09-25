# Au(111) surface defects — constant-potential DFT / electric double layer

Constant-potential DFT (CP-DFT) study of Au(111) surface defects using a
single-sided implicit-solvent CEP-DIP code (VASP + VASPsol++ + CP-VASP,
locally modified), aimed at a training dataset for non-specific anion
enrichment near defects as a function of electrode potential.

**Live report**: https://claude.ai/artifact/C3tDzBtSYwmLaiMbxap2Ft

## What's in this repo

This is the lightweight subset of the project directory — scripts, input
files, analysis, and results. It deliberately excludes:

- **Raw VASP binary outputs** (`WAVECAR`, `CHG`/`CHGCAR`, `vaspout.h5`,
  `PROCAR`, `OUTCAR`, `vasprun.xml`) and the CEP-DIP 3D field dumps
  (`POT`, `PHI`, `RHOB`, `RHOION`, `SDIEL`, `SION`, etc.) — these are large
  (tens of GB in aggregate) and are scratch/intermediate data, not
  something meant to be version-controlled.
- **`POTCAR` files** — VASP's licensed pseudopotential data, not
  redistributable.

The full run directories (including the excluded files above) live on
Purdue Anvil at `/anvil/scratch/x-rywang/Au_Cl`.

## Layout

- `00_audit/` — interface audit of the local CEP-DIP build (environment,
  model conventions, parameter map / decision log)
- `02_bulk_eos/` — bulk Au lattice-constant fit (7-point PBE scan)
- `01_min_validation/`, `03_pilot/`, `04_speedtest/` — input decks
  (`INCAR`/`KPOINTS`/job scripts/`POSCAR`) and lightweight logs for the
  validation run, the 9-point pilot (T/V1/A1-fcc × 3 potentials), and
  parallelization/speed testing
- `03_pilot/all_defect_structures/`, `03_pilot/report_assets/` — generated
  defect-structure geometries and rendered diagrams
- `scripts/` — structure builders, geometry-check utilities, renderers
- `defect_plan.md` — current defect/morphology sampling plan (steps,
  kinks, islands, pits, reconstruction placeholders) with geometry QC
  fields and corrections log
- `reports/pilot_report.html` — source of the published report page above

## Status

Work in progress. See `defect_plan.md` for the current plan and open
items, and `00_audit/parameter_map.md` for the full decision log.
