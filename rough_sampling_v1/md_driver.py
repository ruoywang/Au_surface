#!/usr/bin/env python3
"""Generate (and optionally run) the LAMMPS inputs that evolve each parent surface with the public Au FLARE
potential, following the fixed 200 ps protocol.

Protocol per parent (an engineering exploration budget, not a claim of equilibrated roughness):
    pre-relaxation   : a short minimisation of the movable atoms
    300 K            20 ps
    300 -> 600 K     20 ps   (linear ramp)
    600 K           100 ps
    600 -> 300 K     20 ps
    300 K            40 ps
                    200 ps total, 2 fs step, NVT (Nose-Hoover, 100 fs damping) on the MOVABLE atoms only
Fixed atoms (bottom two layers) are not integrated and have zero velocity; the thermostat and the reported
temperature use the movable group only. One frame is written every 1 ps (200 frames per parent) plus the final
data file. If nothing merges or reconstructs, that is recorded as found: no temperature is raised and no run
is extended to chase it.

The potential is called exactly as its authors specify:  newton on / pair_style flare / pair_coeff * * <file>.
Its energies and forces are candidate-generation tools only and never become CP-DFT labels.

Usage (from Au_Cl/):
    rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/md_driver.py --write            # inputs for all 32 parents
    rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/md_driver.py --write --only PA1_s11 --steps 200 --outdir <dir>  # timing run
    rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/md_driver.py --write --slurm --np 32   # also one SLURM script per parent
Outputs: rough_sampling_v1/md/<parent_id>/{in.lammps, data.lammps, run.sh[, slurm.sh]}, md/md_manifest.json
Run:     NP=<ranks> md/<parent_id>/run.sh        (the binary is MPI; run.sh uses the recorded mpirun)
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np
from ase.io import read

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
sys.path.insert(0, R)
from lmpio import write_data, group_lines  # noqa: E402
POT = f"{R}/potential/Au_training/lmp_t0.0001_no_bulk_vac_fix3.header2025.flare"
LMP = f"{R}/env/src/lammps-22Jul2025/build_mpi/lmp"
MPIRUN = "/apps/spack/anvil/apps/openmpi/4.0.6-gcc-11.2.0-3navcwb/bin/mpirun"
MPILIB = "/apps/spack/anvil/apps/openmpi/4.0.6-gcc-11.2.0-3navcwb/lib"
DT_PS = 0.002
SEGMENTS = [("300 K", 300.0, 300.0, 20.0), ("ramp 300->600 K", 300.0, 600.0, 20.0), ("600 K", 600.0, 600.0, 100.0),
            ("ramp 600->300 K", 600.0, 300.0, 20.0), ("300 K", 300.0, 300.0, 40.0)]
DUMP_EVERY_PS = 1.0
TDAMP_PS = 0.1


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()


def write_inputs(pid, seed, steps_override=None, outdir=None):
    at = read(f"{R}/parents/{pid}.extxyz")
    fixed = at.get_array("fixed").astype(bool)
    d = outdir or f"{R}/md/{pid}"; os.makedirs(d, exist_ok=True)
    # one atom type (pair_style flare maps types to species); the frozen bottom two layers are an id group
    write_data(at, f"{d}/data.lammps", comment=f"{pid}: 32x32x4 Au(111) parent, {int(fixed.sum())} frozen atoms selected by id in in.lammps")

    total_ps = sum(s[3] for s in SEGMENTS); dump_every = int(round(DUMP_EVERY_PS / DT_PS))
    L = [f"# {pid}: Au FLARE candidate generation, 200 ps protocol (2 fs), NVT on mobile atoms only",
         f"# potential sha256 {sha256(POT)}", "units           metal", "atom_style      atomic", "boundary        p p f",
         "newton          on", f"read_data       data.lammps", "", "pair_style      flare",
         f"pair_coeff      * * {POT}", "", "neighbor        2.0 bin", "neigh_modify    every 1 delay 0 check yes", "",
         group_lines("frozen", fixed), "group           mobile subtract all frozen",
         "velocity        frozen set 0.0 0.0 0.0", "fix             hold frozen setforce 0.0 0.0 0.0",
         "compute         Tm mobile temp", "",
         "thermo_style    custom step time temp pe ke etotal press", "thermo_modify   temp Tm", "thermo          500", "",
         "# low-cost pre-relaxation of the movable atoms", "min_style       cg", "minimize        1e-6 1e-4 300 3000", "",
         f"reset_timestep  0", f"timestep        {DT_PS}",
         f"velocity        mobile create 300.0 {seed} dist gaussian mom yes rot yes",
         f"dump            d1 all custom {dump_every} traj.lammpstrj id type x y z", "dump_modify     d1 sort id",
         ""]
    t = 0.0
    for k, (name, T0, T1, ps) in enumerate(SEGMENTS):
        n = int(round(ps / DT_PS))
        if steps_override is not None: n = min(n, steps_override)
        L += [f"# segment {k+1}: {name}, {ps:.0f} ps", f"fix             nvt mobile nvt temp {T0} {T1} {TDAMP_PS}",
              "fix_modify      nvt temp Tm", f"run             {n}", "unfix           nvt", ""]
        t += ps
        if steps_override is not None: break
    L += ["write_data      final.data", f"print           \"DONE {pid} total_ps {t:.0f}\""]
    open(f"{d}/in.lammps", "w").write("\n".join(L) + "\n")
    open(f"{d}/run.sh", "w").write(f"#!/bin/bash\n# NP=<ranks> {d}/run.sh\nexport LD_LIBRARY_PATH={MPILIB}:$LD_LIBRARY_PATH\n"
                                   f"cd {d}\n{MPIRUN} -np ${{NP:-1}} {LMP} -in in.lammps -log log.lammps > stdout.txt 2>&1\n")
    os.chmod(f"{d}/run.sh", 0o755)
    return dict(parent_id=pid, dir=d, n_atoms=len(at), n_mobile=int((~fixed).sum()), n_fixed=int(fixed.sum()),
                seed=seed, total_ps=t, steps=int(round(t / DT_PS)) if steps_override is None else steps_override)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true"); ap.add_argument("--only", default=None)
    ap.add_argument("--steps", type=int, default=None, help="cap each segment at this many steps (timing runs)")
    ap.add_argument("--slurm", action="store_true"); ap.add_argument("--outdir", default=None)
    ap.add_argument("--np", type=int, default=32, help="MPI ranks per run in the SLURM scripts")
    ap.add_argument("--walltime", default="12:00:00")
    a = ap.parse_args()
    M = json.load(open(f"{R}/parents/parents_manifest.json"))
    rows = []
    for m in M:
        if a.only and m["parent_id"] != a.only: continue
        rows.append(write_inputs(m["parent_id"], 1000 + m["seed"], a.steps, a.outdir))
    info = dict(potential=POT, potential_sha256=sha256(POT), lammps_binary=LMP, mpirun=MPIRUN, dt_ps=DT_PS, segments=SEGMENTS,
                dump_every_ps=DUMP_EVERY_PS, tdamp_ps=TDAMP_PS, runs=rows)
    if a.slurm:
        info["slurm"] = dict(partition="shared", account="CHE190065", ranks=a.np, walltime=a.walltime)
        for r in rows:
            open(f"{r['dir']}/slurm.sh", "w").write(
                "#!/bin/bash\n#SBATCH -J md_" + r["parent_id"] + "\n#SBATCH --account=CHE190065\n#SBATCH --partition=shared\n"
                f"#SBATCH --nodes=1\n#SBATCH --ntasks={a.np}\n#SBATCH --cpus-per-task=1\n#SBATCH --time={a.walltime}\n"
                f"#SBATCH -o {r['dir']}/slurm.out\n#SBATCH -e {r['dir']}/slurm.err\n"
                f"export LD_LIBRARY_PATH={MPILIB}:$LD_LIBRARY_PATH\ncd {r['dir']}\n"
                f"{MPIRUN} -np {a.np} {LMP} -in in.lammps -log log.lammps > stdout.txt 2>&1\n")
    if a.outdir is None:
        os.makedirs(f"{R}/md", exist_ok=True)
        json.dump(info, open(f"{R}/md/md_manifest.json", "w"), indent=1)
    print(f"{len(rows)} input sets written under {R}/md/  (protocol {sum(s[3] for s in SEGMENTS):.0f} ps, "
          f"{int(sum(s[3] for s in SEGMENTS)/DT_PS)} steps each)")


if __name__ == "__main__":
    main()
