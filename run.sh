#!/bin/bash

#SBATCH --job-name=iso_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --threads-per-core=1
#SBATCH --time=WALLTIME
#SBATCH --partition=devel
#SBATCH --output basilisk.output
#SBATCH --exclusive

LEVEL=7

module purge
module load mpi/openmpi-x86_64
module load libGLU/9.0.1-GCCcore-10.3.0
module load Mesa/21.1.1-GCCcore-10.3.0
module load intel/2021a

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so

mpicc -Wall -std=c99 -O2 _isotropic.c -o isotropic -I$BASILISK -L$BASILISK/gl -L$/home/reub0024/local/lib -lglutils -lfb_osmesa -lOSMesa -lGLU -lm
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS \
     ./isotropic -m WALLTIME $LEVEL \
     2> log-$LEVEL-$SLURM_NTASKS > out-$LEVEL-$SLURM_NTASKS
