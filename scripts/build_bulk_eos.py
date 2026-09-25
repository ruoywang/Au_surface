#!/usr/bin/env python3
"""Build a 7-point fcc Au bulk EOS scan (v2 doc Sec.3.1) to replace the
a0=4.180 A literature placeholder used in 01_min_validation/ with a locally
fitted PBE value."""
import os
from ase.build import bulk
from ase.io import write

A_VALUES = [4.00, 4.05, 4.10, 4.15, 4.20, 4.25, 4.30]  # Angstrom, brackets typical PBE Au (~4.15-4.18)
OUTBASE = "02_bulk_eos"

for a in A_VALUES:
    tag = f"a{a:.2f}"
    d = f"{OUTBASE}/{tag}"
    os.makedirs(d, exist_ok=True)
    atoms = bulk("Au", "fcc", a=a)
    write(f"{d}/POSCAR", atoms, format="vasp", direct=True, sort=True)
    print(tag, atoms.get_volume(), "A^3/atom")
