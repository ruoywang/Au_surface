#!/usr/bin/env python3
"""Deliberate small-sample split for the INTERFACE test (not a generalisation test):
train = the four Flat16x1 / Step-16x1 states (2 geometries x 2 charge states),
valid = Step-8x1_muref_fastcfg (1 state).  Same-geometry states are kept on the same
side of the split on purpose; nothing here is a held-out test of transferability.
Run with the trainer env's python (ASE 3.25) or scripts/pyrun.sh.
"""
from ase.io import read, write

src = "../v0_main_production.mace.extxyz"
cfgs = read(src, ":")
train = [a for a in cfgs if a.info["state_id"] in ("Flat16x1_muref", "Flat16x1_dUp02", "Step-16x1_muref_fastcfg", "Step-16x1_dUp02_fastcfg")]
valid = [a for a in cfgs if a.info["state_id"] == "Step-8x1_muref_fastcfg"]
assert len(train) == 4 and len(valid) == 1
write("train.extxyz", train, format="extxyz")
write("valid.extxyz", valid, format="extxyz")
for name, lst in (("train", train), ("valid", valid)):
    for a in lst:
        print(f"{name}: {a.info['state_id']:26s} N={len(a):3d} electron={a.info['electron']:.6f} potential={a.info['potential']:.6f} "
              f"E_woS={a.info['E_without_entropy_eV']:.5f} |F|max={abs(a.arrays['forces']).max():.3f}")
