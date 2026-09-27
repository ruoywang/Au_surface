#!/usr/bin/env python3
"""DIAGNOSTIC variant of the interface test (2026-09-27): identical to interface_test/ except that
the fork's hard-coded `electron` input is written as the EXCESS electron count
    electron = N_e_final - N_neutral   (|.| <= 0.14 here)
instead of the absolute N_e (703.9 ... 792.0, i.e. +-44 around the train-set mean when the two
geometries have 64 and 72 atoms). The absolute count is kept as `electron_Ne_actual` for provenance.
Purpose: test whether the immediate divergence of interface_test (RMSE_E 5.2 eV/atom before any
update) is caused by the scale of the `electron` node attribute. This is a data-convention
question for the user, not a model change; nothing in the fork is modified.
"""
from ase.io import read, write
cfgs = read("../v0_main_production.mace.extxyz", ":")
for a in cfgs:
    a.info["electron_Ne_actual"] = a.info["electron"]
    a.info["N_neutral"] = 11.0 * len(a)
    a.info["electron"] = a.info["electron"] - 11.0 * len(a)
train = [a for a in cfgs if a.info["state_id"] in ("Flat16x1_muref", "Flat16x1_dUp02", "Step-16x1_muref_fastcfg", "Step-16x1_dUp02_fastcfg")]
valid = [a for a in cfgs if a.info["state_id"] == "Step-8x1_muref_fastcfg"]
assert len(train) == 4 and len(valid) == 1
write("train.extxyz", train, format="extxyz"); write("valid.extxyz", valid, format="extxyz")
for name, lst in (("train", train), ("valid", valid)):
    for a in lst:
        print(f"{name}: {a.info['state_id']:26s} N={len(a):3d} electron(excess)={a.info['electron']:+.6f} Ne={a.info['electron_Ne_actual']:.6f} potential={a.info['potential']:.6f}")
