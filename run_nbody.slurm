#!/bin/bash
#SBATCH --job-name=nbody_sim
#SBATCH --partition=GPU
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --output=nbody_output_%j.log

# Usage: sbatch run_nbody.slurm <num_particles> <dt> <num_steps> <print_interval> <block_size>

NUM_PARTICLES=${1:-1000}
DT=${2:-0.01}
NUM_STEPS=${3:-50}
PRINT_EVERY=${4:-10}
BLOCK_SIZE=${5:-128}

module load cuda/12.4

# Add GPU flag (the 6th argument)
./nbody $NUM_PARTICLES $DT $NUM_STEPS $PRINT_EVERY $BLOCK_SIZE --gpu