#!/bin/bash
#BSUB -n XXXX
#BSUB -o %MD.out -e %MD.err
#BSUB -J YYYY
#BSUB -q training
#BSUB -W 210
#BSUB -U PROJ_I
#BSUB -R "span[ptile=XXXX]"


source ${MODULESHOME}/init/bash
module purge
module load openmpi
module load gcc/3.3

touch module_loaded.out
echo $(module list) > module_loaded.out

num_procs=$1

mpirun -np $(num_procs) dynamics_mpi input.dat
