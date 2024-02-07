#!/bin/bash
#SBATCH --partition=hi-core             # partition for MPI
#SBATCH -N 2                          # Can request up to 16 node
#SBATCH --time=01:30:00                 # Can request upto 6 hours
#SBATCH --constraint='epyc128'
#SBATCH --mem=128G

mpicxx -std=c++17 -O3  -o Demo Demo.cpp
srun --mpi=pmix ./Demo

# ds7_1M_fl --ntasks-per-node=1 --cpus-per-task=6 SBATCH --ntasks=2      Can request upto 384 CPU cores
