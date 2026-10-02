#!/usr/bin/env python3
"""Regression cases for fieldio.check on small synthetic CHGCAR-like files (2 ions, 2x2x2 grid). Every case states
the verdict it must get; the script exits non-zero if any verdict differs. Includes the two counterexamples from the
review of 2026-10-02 (truncated last augmentation block; NaN inside a block) and the plain cases.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/test_fieldio.py
"""
import os
import sys
import tempfile

sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1")
import fieldio  # noqa: E402

HEAD = """Au test
   1.00000000000000
     5.0  0.0  0.0
     0.0  5.0  0.0
     0.0  0.0 10.0
   Au
     2
Direct
  0.0 0.0 0.1
  0.5 0.5 0.3

    2    2    2
"""
GRID_OK = " 0.1 0.2 0.3 0.4 0.5\n 0.6 0.7 0.8\n"
GRID_SHORT = " 0.1 0.2 0.3 0.4 0.5\n 0.6\n"
GRID_NAN = " 0.1 0.2 nan 0.4 0.5\n 0.6 0.7 0.8\n"
AUG_OK = "augmentation occupancies   1   5\n 0.1 0.2 0.3 0.4 0.5\naugmentation occupancies   2   5\n 0.1 0.2 0.3 0.4 0.5\n"
AUG_TRUNC = "augmentation occupancies   1   5\n 0.1 0.2 0.3 0.4 0.5\naugmentation occupancies   2   5\n 0.1 0.2\n"
AUG_NAN = "augmentation occupancies   1   5\n 0.1 0.2 0.3 0.4 0.5\naugmentation occupancies   2   5\n 0.1 nan 0.3 0.4 0.5\n"
AUG_MISSING_ION = "augmentation occupancies   1   5\n 0.1 0.2 0.3 0.4 0.5\n"

CASES = [  # (name as given to check, body, expected ok, label)
    ("CHGCAR", HEAD + GRID_OK + AUG_OK, True, "complete CHGCAR"),
    ("CHGCAR", HEAD + GRID_OK + AUG_TRUNC, False, "last augmentation block declares 5 values, has 2"),
    ("CHGCAR", HEAD + GRID_OK + AUG_NAN, False, "NaN inside an augmentation block"),
    ("CHGCAR", HEAD + GRID_OK + AUG_MISSING_ION, False, "one augmentation block for two ions"),
    ("CHGCAR", HEAD + GRID_OK, False, "augmented field without augmentation blocks"),
    ("CHGCAR", HEAD + GRID_SHORT + AUG_OK, False, "main grid truncated (6 of 8 values)"),
    ("CHGCAR", HEAD + GRID_NAN + AUG_OK, False, "NaN in the main grid"),
    ("RHOB", HEAD + GRID_OK, True, "complete non-augmented field"),
    ("RHOB", HEAD + GRID_SHORT, False, "non-augmented field truncated"),
    ("RHOB", HEAD + GRID_OK + AUG_OK, False, "augmentation block in a non-augmented field"),
]


def main():
    bad = 0
    with tempfile.TemporaryDirectory() as d:
        for k, (name, body, expect, label) in enumerate(CASES):
            p = os.path.join(d, f"{name}_{k}"); open(p, "w").write(body)
            ok, note, _ = fieldio.check(p, name)
            verdict = "PASS" if ok == expect else "FAIL"
            if ok != expect: bad += 1
            print(f"[{verdict}] {label:55s} -> check={ok!s:5s} expected={expect!s:5s} ({note})")
        p = os.path.join(d, "missing"); ok, note, _ = fieldio.check(p, "PHI"); print(f"[{'PASS' if not ok else 'FAIL'}] {'missing file':55s} -> ({note})"); bad += int(ok)
    print(f"{len(CASES) + 1 - bad} of {len(CASES) + 1} cases as expected")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
