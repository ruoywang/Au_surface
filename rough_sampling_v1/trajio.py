"""Trajectory access for the rough-sampling scripts: frames are addressed by their LAMMPS TIMESTEP, never by
array index, and every frame carries the thermostat stage and the measured temperature of the mobile group.

The protocol (md_driver.SEGMENTS) gives the stage of every time; the dump gives the timesteps actually written;
log.lammps gives the temperature at every thermo step. A frame is accepted for sampling only if the requested
time exists in the dump, lies in a 300 K hold, and the logged temperature is within T_TOL of 300 K.
"""
import os
import re

import numpy as np
from ase.io import read

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
import sys  # noqa: E402
sys.path.insert(0, R)
from md_driver import SEGMENTS, DT_PS  # noqa: E402

TOTAL_PS = sum(s[3] for s in SEGMENTS)
TOTAL_STEPS = int(round(TOTAL_PS / DT_PS))
T_TOL = 50.0
MD_DIR = f"{R}/md"          # overridable by test drivers (size_complex_test_v1 keeps its own md/)


class TrajError(RuntimeError):
    pass


def stage_of(t_ps):
    """(name, T0, T1) of the protocol segment that contains time t_ps."""
    t = 0.0
    for name, T0, T1, ps in SEGMENTS:
        if t <= t_ps <= t + ps + 1e-9: return name, T0, T1
        t += ps
    raise TrajError(f"time {t_ps} ps is outside the {TOTAL_PS:.0f} ps protocol")


def dump_steps(path):
    """All TIMESTEP values in a LAMMPS text dump, in file order (a cheap line scan; no coordinates read)."""
    steps = []
    with open(path) as f:
        for line in f:
            if line.startswith("ITEM: TIMESTEP"):
                steps.append(int(next(f).split()[0]))
    return steps


def thermo_table(log_path):
    """step -> (time_ps, temp_K) from every thermo block of log.lammps (custom step time temp ...)."""
    out = {}
    if not os.path.exists(log_path): return out
    cols = None
    for line in open(log_path):
        p = line.split()
        if not p: continue
        if p[0] == "Step": cols = p; continue
        if cols and p[0].lstrip("-").isdigit() and len(p) == len(cols):
            try: out[int(p[0])] = (float(p[cols.index("Time")]), float(p[cols.index("Temp")]))
            except (ValueError, IndexError): pass
            continue
        if cols and p[0] == "Loop": cols = None
    return out


def run_complete(pid):
    """(ok, why): the run printed its DONE line for the full protocol and the dump holds the last step."""
    d = f"{MD_DIR}/{pid}"
    so = f"{d}/stdout.txt"; tr = f"{d}/traj.lammpstrj"
    if not os.path.exists(tr): return False, "no trajectory"
    if not (os.path.exists(so) and re.search(rf"^DONE {re.escape(pid)} total_ps {TOTAL_PS:.0f}\b", open(so).read(), re.M)):
        return False, "no DONE line for the full protocol in stdout.txt"
    steps = dump_steps(tr)
    if not steps or steps[-1] != TOTAL_STEPS: return False, f"dump ends at step {steps[-1] if steps else None}, expected {TOTAL_STEPS}"
    return True, f"{len(steps)} frames, last step {steps[-1]}"


def frames_at(pid, times_ps, require_300K=True):
    """Frames of parent pid at the requested times, each as (step, time_ps, T_K, Atoms). Raises TrajError naming
    the parent and the missing or off-stage time; never falls back."""
    d = f"{MD_DIR}/{pid}"; tr = f"{d}/traj.lammpstrj"
    ok, why = run_complete(pid)
    if not ok: raise TrajError(f"{pid}: run incomplete ({why})")
    steps = dump_steps(tr); pos = {s: k for k, s in enumerate(steps)}
    th = thermo_table(f"{d}/log.lammps")
    out = []
    for t in times_ps:
        step = int(round(t / DT_PS))
        if step not in pos: raise TrajError(f"{pid}: no frame at {t} ps (step {step}) in the dump")
        name, T0, T1 = stage_of(t)
        if require_300K and not (T0 == 300.0 and T1 == 300.0): raise TrajError(f"{pid}: {t} ps is in segment '{name}', not a 300 K hold")
        if step not in th: raise TrajError(f"{pid}: no thermo row at step {step} in log.lammps")
        T = th[step][1]
        if require_300K and abs(T - 300.0) > T_TOL: raise TrajError(f"{pid}: T = {T:.0f} K at {t} ps, outside 300 +- {T_TOL:.0f} K")
        at = read(tr, index=pos[step], format="lammps-dump-text")
        at.set_chemical_symbols(["Au"] * len(at)); at.set_pbc((True, True, False))
        out.append((step, float(t), float(T), at))
    return out


def frame_by_step(pid, step):
    """One frame by TIMESTEP (for the cut and the figures)."""
    tr = f"{MD_DIR}/{pid}/traj.lammpstrj"
    steps = dump_steps(tr)
    if step not in steps: raise TrajError(f"{pid}: step {step} not in dump")
    at = read(tr, index=steps.index(step), format="lammps-dump-text")
    at.set_chemical_symbols(["Au"] * len(at)); at.set_pbc((True, True, False)); return at
