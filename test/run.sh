#!/bin/bash
#SBATCH -J CTSEG
#SBATCH -p ccq
#SBATCH --constraint=rome
#SBATCH -N4
#SBATCH --ntasks=448
#SBATCH --ntasks-per-core=1
##SBATCH --exclusive
#SBATCH --time=0-00:10:00
##SBATCH --mail-type=all


module load triqs/3_unst_llvm_ompi

export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=TRUE
export OMP_NESTED=FALSE
export OMP_WAIT_POLICY=ACTIVE
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

mpirun -np 448 --map-by core python test_J.py