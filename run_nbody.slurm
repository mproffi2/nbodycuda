#!/bin/bash
#SBATCH --job-name=nbody
#SBATCH --output=nbody_%j.out
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=00:10:00

module load cuda/12.4

# Example run:
#   ./nbody <input_or_number> <dt> <steps> <print_every> <blocksize> cpu|gpu

./nbody 10000 0.001 5 5 128 gpu