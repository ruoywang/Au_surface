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
    scripts/pyrun.sh rough_sampling_v1/md_driver.py --write            # inputs for all 32 parents
    scripts/pyrun.sh rough_sampling_v1/md_driver.py --write --only PA1_s11 --steps 200   # a short timing run
    scripts/pyrun.sh rough_sampling_v1/md_driver.py --slurm            # also write one SLURM script per parent
Outputs: rough_sampling_v1/md/<parent_id>/{in.lammps, data.lammps, run.sh}, md/README.md
"""
import argparse
import hashlib
import json
import os
import subprocess

import numpy as np
from ase.io import read, write

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
POT = f"{R}/potential/Au_training/lmp_t0.0001_no_bulk_vac_fix3.flare"
LMP = f"{R}/env/src/lammps-22Jul2025/build/lmp"
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
    # LAMMPS atom types carry the fixed/mobile split so the groups are unambiguous in the data file itself:
    # type 1 = mobile Au, type 2 = fixed Au (both mass 196.967, both the same species to the potential)
    at.set_tags(np.where(fixed, 2, 1))
    write(f"{d}/data.lammps", at, format="lammps-data", atom_style="atomic", masses=True,
          specorder=None, units="metal")
    # ase writes one type per species; rewrite the Atoms section to use the tag as the type and declare 2 types
    txt = open(f"{d}/data.lammps").read().splitlines()
    out = []; in_atoms = False; tags = at.get_tags()
    for line in txt:
        if line.strip().endswith("atom types"): line = "2 atom types"
        if line.startswith("Masses"):
            out += ["Masses", "", "1 196.967   # Au, mobile", "2 196.967   # Au, fixed (bottom two layers)", ""]
            skip = True; continue
        if line.startswith("Atoms"): in_atoms = True; out.append(line); continue
        if in_atoms and line.strip():
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                parts[1] = str(int(tags[int(parts[0]) - 1])); line = " ".join(parts)
        if line.startswith("Velocities"): in_atoms = False
        out.append(line)
    # drop the masses block ase wrote (we re-emitted it); crude but the format is simple
    cleaned = []; i = 0
    while i < len(out):
        if out[i] == "Masses" and i + 1 < len(out) and out[i + 1] == "" and i + 2 < len(out) and out[i + 2].startswith("1 196.967   # Au, mobile"):
            cleaned += out[i:i + 5]; i += 5
            # skip any following original masses lines ("1 196.967")
            while i < len(out) and (out[i].strip() == "" or out[i].strip().startswith("1 196.96")): i += 1
            continue
        cleaned.append(out[i]); i += 1
    open(f"{d}/data.lammps", "w").write("\n".join(cleaned) + "\n")

    total_ps = sum(s[3] for s in SEGMENTS); dump_every = int(round(DUMP_EVERY_PS / DT_PS))
    L = [f"# {pid}: Au FLARE candidate generation, 200 ps protocol (2 fs), NVT on mobile atoms only",
         f"# potential sha256 {sha256(POT)}", "units           metal", "atom_style      atomic", "boundary        p p f",
         "newton          on", f"read_data       data.lammps", "", "pair_style      flare",
         f"pair_coeff      * * {POT}", "", "neighbor        2.0 bin", "neigh_modify    every 1 delay 0 check yes", "",
         "group           mobile type 1", "group           frozen type 2",
         "velocity        frozen set 0.0 0.0 0.0", "fix             hold frozen setforce 0.0 0.0 0.0",
         "compute         Tm mobile temp", "thermo_modify   temp Tm", "",
         "thermo_style    custom step time temp pe ke etotal press", "thermo          500", "",
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
    open(f"{d}/run.sh", "w").write(f"#!/bin/bash\ncd {d}\n{LMP} -in in.lammps -log log.lammps > stdout.txt 2>&1\n")
    os.chmod(f"{d}/run.sh", 0o755)
    return dict(parent_id=pid, dir=d, n_atoms=len(at), n_mobile=int((~fixed).sum()), n_fixed=int(fixed.sum()),
                seed=seed, total_ps=t, steps=int(round(t / DT_PS)) if steps_override is None else steps_override)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true"); ap.add_argument("--only", default=None)
    ap.add_argument("--steps", type=int, default=None, help="cap each segment at this many steps (timing runs)")
    ap.add_argument("--slurm", action="store_true"); ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    M = json.load(open(f"{R}/parents/parents_manifest.json"))
    rows = []
    for m in M:
        if a.only and m["parent_id"] != a.only: continue
        rows.append(write_inputs(m["parent_id"], 1000 + m["seed"], a.steps, a.outdir))
    info = dict(potential=POT, potential_sha256=sha256(POT), lammps_binary=LMP, dt_ps=DT_PS, segments=SEGMENTS,
                dump_every_ps=DUMP_EVERY_PS, tdamp_ps=TDAMP_PS, runs=rows)
    if a.slurm:
        for r in rows:
            open(f"{r['dir']}/slurm.sh", "w").write(
                "#!/bin/bash\n#SBATCH -J md_" + r["parent_id"] + "\n#SBATCH -A che190065\n#SBATCH -p shared\n"
                "#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH -t 48:00:00\n"
                f"#SBATCH -o {r['dir']}/slurm.out\n#SBATCH -e {r['dir']}/slurm.err\n"
                f"cd {r['dir']}\n{LMP} -in in.lammps -log log.lammps > stdout.txt 2>&1\n")
    os.makedirs(f"{R}/md", exist_ok=True)
    json.dump(info, open(f"{R}/md/md_manifest.json", "w"), indent=1)
    print(f"{len(rows)} input sets written under {R}/md/  (protocol {sum(s[3] for s in SEGMENTS):.0f} ps, "
          f"{int(sum(s[3] for s in SEGMENTS)/DT_PS)} steps each)")


if __name__ == "__main__":
    main()
