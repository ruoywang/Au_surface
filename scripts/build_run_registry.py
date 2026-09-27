#!/usr/bin/env python3
"""Provenance registry for every CP-DFT run directory under 03_pilot/ (and 04_speedtest/):
numerical configuration + final CPM state + timing, so that results from different PREC /
k-mesh / parallel configurations are never merged silently (2026-09-27 decision: keep old
and new PREC data, record the configuration version, do not apply a blanket meV/atom
correction when energy labels are combined).

Output: 03_pilot/run_registry.json and 03_pilot/run_registry.md.
Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_run_registry.py
"""
import glob
import json
import os
import re
import subprocess

from ase.io import read

TAGS = ["PREC", "ENCUT", "ALGO", "EDIFF", "ISMEAR", "SIGMA", "NPAR", "NCORE", "KPAR", "LREAL", "ISYM",
        "TARGETMU", "FERMICONVERGE", "NESCHEME", "C_MOLAR", "R_ION", "SOL_Z0", "SOL_Z1", "SOL_SIGMA", "DIPOL",
        "AMIX", "BMIX", "AMIN", "NELECT", "IBRION", "NSW", "LVAC"]


def incar_tags(path):
    tags = {}
    for line in open(path):
        line = line.split("#")[0].split("!")[0].strip()
        if "=" in line:
            k, v = [x.strip() for x in line.split("=", 1)]
            if k.upper() in TAGS:
                tags[k.upper()] = v
    return tags


def kpoints(path):
    if not os.path.exists(path):
        return None
    lines = open(path).read().splitlines()
    return " ".join(lines[3].split()[:3]) if len(lines) > 3 else None


def parallel_layout(jobrun):
    if not os.path.exists(jobrun):
        return None
    txt = open(jobrun).read()
    nt = re.search(r"--ntasks=(\d+)", txt); nc = re.search(r"--cpus-per-task=(\d+)", txt)
    nn = re.search(r"--nodes=(\d+)", txt); part = re.search(r"--partition=(\S+)", txt)
    return dict(nodes=int(nn.group(1)) if nn else None, mpi_ranks=int(nt.group(1)) if nt else None,
                threads_per_rank=int(nc.group(1)) if nc else 1, partition=part.group(1) if part else None)


def final_state(d):
    rec = dict(cpm_ion_lines=0)
    log = f"{d}/log.out"
    if os.path.exists(log):
        txt = open(log).read()
        m = re.findall(r"CPM-ion:.*?N_ele=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)\s+TARGETMU=\s*([-\d.]+)", txt)
        rec["cpm_ion_lines"] = len(m)
        if m:
            rec.update(N_e_final=float(m[-1][0]), mu_e_final=float(m[-1][1]), targetmu_log=float(m[-1][2]))
    osz = f"{d}/OSZICAR"
    if os.path.exists(osz):
        rec["scf_steps"] = sum(1 for l in open(osz) if l.startswith(("DAV", "RMM")))
    out = f"{d}/OUTCAR"
    if os.path.exists(out):
        toten = None; elapsed = None
        for line in open(out):
            if "free  energy   TOTEN" in line:
                toten = float(line.split("=")[1].split()[0])
            elif "Elapsed time (sec)" in line:
                elapsed = float(line.split(":")[1])
        rec.update(TOTEN_eV=toten, elapsed_s=elapsed)
    return rec


def slurm_state(d):
    ids = sorted(int(re.search(r"myjob\.o(\d+)", f).group(1)) for f in glob.glob(f"{d}/myjob.o*"))
    if not ids:
        return dict(job_ids=[], last_state=None)
    try:
        st = subprocess.run(["sacct", "-j", str(ids[-1]), "-n", "-X", "-o", "State,Elapsed"], capture_output=True, text=True, timeout=30).stdout.split()
        return dict(job_ids=ids, last_state=st[0] if st else None, last_elapsed=st[1] if len(st) > 1 else None)
    except Exception:
        return dict(job_ids=ids, last_state=None)


reg = {}
for d in sorted(glob.glob("03_pilot/*/") + glob.glob("04_speedtest/*/")):
    if not (os.path.exists(f"{d}/INCAR") and os.path.exists(f"{d}/POSCAR")):
        continue
    if os.path.basename(d.rstrip("/")) in ("structures", "all_defect_structures", "report_assets"):
        continue
    name = d.rstrip("/").split("/", 1)[1]
    atoms = read(f"{d}/POSCAR")
    tags = incar_tags(f"{d}/INCAR")
    rec = dict(dir=d.rstrip("/"), n_atoms=len(atoms), cell_in_plane=[[round(x, 4) for x in atoms.get_cell()[i][:2]] for i in range(2)],
               metal_top_z=round(float(atoms.get_positions()[:, 2].max()), 3), kpoints=kpoints(f"{d}/KPOINTS"),
               incar=tags, parallel=parallel_layout(f"{d}/job-run"))
    rec.update(final_state(d)); rec.update(slurm_state(d))
    rec["config_version"] = ("production-2026-09-26 (PREC=Normal, ALGO=Fast, NPAR=16 hybrid)" if tags.get("PREC", "").lower().startswith("n") and tags.get("ALGO", "").lower() == "fast"
                             else "pilot-2026-09-22 (PREC=Accurate, ALGO=Normal, NCORE=8)" if tags.get("PREC", "").lower().startswith("a")
                             else "other/test")
    reg[name] = rec

json.dump(reg, open("03_pilot/run_registry.json", "w"), indent=1)

cols = ["run", "atoms", "k-mesh", "PREC", "ALGO", "TARGETMU", "FERMICONV", "ranks×thr", "SCF", "N_e final", "μ_e final", "TOTEN (eV)", "wall", "state"]
rows = []
for name, r in reg.items():
    p = r.get("parallel") or {}
    rows.append([name, r["n_atoms"], r.get("kpoints") or "-", r["incar"].get("PREC", "-"), r["incar"].get("ALGO", "-"),
                 r["incar"].get("TARGETMU", "-"), r["incar"].get("FERMICONVERGE", "-"),
                 f"{p.get('mpi_ranks','-')}×{p.get('threads_per_rank','-')}", r.get("scf_steps", "-"),
                 f"{r['N_e_final']:.6f}" if "N_e_final" in r else "-", f"{r['mu_e_final']:.6f}" if "mu_e_final" in r else "-",
                 f"{r['TOTEN_eV']:.5f}" if r.get("TOTEN_eV") is not None else "-", r.get("last_elapsed") or "-", r.get("last_state") or "-"])
md = ["# Run registry (auto-generated by scripts/build_run_registry.py)", "",
      "One row per run directory with INCAR+POSCAR. `config_version` in the JSON separates the pilot configuration "
      "(PREC=Accurate, ALGO=Normal, 128 MPI, NCORE=8) from the production configuration adopted 2026-09-26 "
      "(PREC=Normal, ALGO=Fast, NPAR=16, 16 MPI × 8 OpenMP). TOTEN values from different PREC settings carry a "
      "grid-dependent offset (~0.3 meV/atom measured on Step-8x1/16x1) and are NOT to be mixed as energy labels "
      "without recording this field; no blanket correction is applied. TARGETMU is the internal μ₀ reference, not V vs RHE.", "",
      "**What this table is and is not (2026-09-27).** It is an index for locating raw runs, their actual converged state "
      "(N_e, μ_e from the CPM-ion lines) and their numerical configuration. It is NOT a list of training samples: the rows "
      "include completed, timed-out, cancelled and speed-test runs, and repeated calculations of the same geometry in "
      "different directories are not additional structural coverage. Force, 3D electron-density or solvent-field labels must "
      "come from the raw outputs of the corresponding run, not from the integrated/averaged analysis JSONs in 03_pilot/. "
      "The first step of the next stage is curation of these runs (which to keep, which label standard, which training "
      "targets the existing outputs support), not new DFT.", "",
      "| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
md += ["| " + " | ".join(str(x) for x in row) + " |" for row in rows]
open("03_pilot/run_registry.md", "w").write("\n".join(md) + "\n")
print(f"{len(reg)} runs registered -> 03_pilot/run_registry.json, 03_pilot/run_registry.md")
for row in rows:
    print(" | ".join(str(x) for x in row))
