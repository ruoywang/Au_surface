#!/bin/bash
#SBATCH -J md_PB2_s51
#SBATCH --account=CHE190065
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH -o /anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1/md/PB2_s51/slurm.out
#SBATCH -e /anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1/md/PB2_s51/slurm.err
export LD_LIBRARY_PATH=/apps/spack/anvil/apps/openmpi/4.0.6-gcc-11.2.0-3navcwb/lib:$LD_LIBRARY_PATH
cd /anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1/md/PB2_s51
/apps/spack/anvil/apps/openmpi/4.0.6-gcc-11.2.0-3navcwb/bin/mpirun -np 16 /anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1/env/src/lammps-22Jul2025/build_mpi/lmp -in in.lammps -log log.lammps > stdout.txt 2>&1
