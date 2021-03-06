#!/bin/bash
####################################
#
# Performance profiling on the cluster 
#
# with one argument: the maximal number of processes
#
####################################

# request resources
#SBATCH --job-name=submission
#SBATCH --output=result_processes_${1}.txt
#
#SBATCH --ntasks=$1
#SBATCH --ntasks-per-node=4
#SBATCH --time=10:00

# load modules
module use /usr/local.nfs/sgs/modulefiles
module load gcc/10.2
module load openmpi/3.1.6-gcc-10.2
module load vtk/9.0.1
module load cmake/3.18.2

# enable custom build openmpi 3.1 that works with slurm
export CPATH=/scratch-nfs/maierbn/openmpi/install-3.1/include
export PATH=/scratch-nfs/maierbn/openmpi/install-3.1/bin:$PATH

cd build
srun -n $1 ./numsim_parallel ../ini/lid_driven_cavity.txt
