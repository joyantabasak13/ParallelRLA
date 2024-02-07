#!/bin/bash
#SBATCH --partition=hi-core         # partition for MPI
#SBATCH -N 1                        # Can request up to 16 node
#SBATCH --ntasks=1                  # Can request upto 384 CPU cores
#SBATCH --time=01:30:00             # Can request upto 6 hours
#SBATCH --constraint='epyc128'
#SBATCH --mem=128G

mpicxx -std=c++17 -O3  -o prla_superBlocking_v3 prla_superBlocking_v3.cpp
srun --mpi=openmpi ./prla_superBlocking_v3 ds1_50k
