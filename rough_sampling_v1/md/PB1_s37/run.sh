#!/bin/bash
# NP=<ranks> /anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1/md/PB1_s37/run.sh
export LD_LIBRARY_PATH=/apps/spack/anvil/apps/openmpi/4.0.6-gcc-11.2.0-3navcwb/lib:$LD_LIBRARY_PATH
cd /anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1/md/PB1_s37
/apps/spack/anvil/apps/openmpi/4.0.6-gcc-11.2.0-3navcwb/bin/mpirun -np ${NP:-1} /anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1/env/src/lammps-22Jul2025/build_mpi/lmp -in in.lammps -log log.lammps > stdout.txt 2>&1
