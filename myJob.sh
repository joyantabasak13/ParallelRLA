#!/bin/bash
#SBATCH --partition=hi-core             # partition for MPI
#SBATCH -N 2                            # Can request up to 16 node
#SBATCH --ntasks=2                      # Can request upto 384 CPU cores
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=6
#SBATCH --time=01:30:00                 # Can request upto 6 hours
#SBATCH --constraint='epyc128'
#SBATCH --mem=128G

mpicxx -std=c++17 -O3 -g -o prla_dynamic_hybrid prla_dynamic_hybrid.cpp
srun --mpi=pmix ./prla_dynamic_hybrid ds1_50k_fl

# ds7_1M_fl --ntasks-per-node=1 --cpus-per-task=6 SBATCH -N 2 
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -g
