1. Make sure you have an NVIDIA GPU and CUDA 12.4 installed.
2. Compile the program by running:
   make
3. Submit the job to Centaurus using SLURM:
   sbatch run_nbody.slurm <num_particles> <dt> <num_steps> <print_interval> <block_size>
   Example:
   sbatch run_nbody.slurm 1000 0.01 50 10 128
4. Check the output log file for simulation results:
   nbody_output_<jobid>.log