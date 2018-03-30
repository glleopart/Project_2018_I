#!/bin/bash
#BSUB -n XXXX
#BSUB -o %MD.out -e %MD.err
#BSUB -J YYYY
#BSUB -q training
#BSUB -W 01:00
#BSUB -U PROJI
#BSUB -R "span[ptile=16]"

num_procs=$1

mpirun -np $(num_procs) dynamics_mpi input.dat
