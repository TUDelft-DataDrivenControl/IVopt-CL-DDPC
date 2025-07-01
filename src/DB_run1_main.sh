#!/bin/bash
#SBATCH --job-name=seed94
#SBATCH --partition=compute
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900M
#SBATCH --account=research-me-dcsc
#SBATCH --output=./Job1.out
#SBATCH --error=./Job1.err

# Load any necessary modules or set environment variables
#research-me-dcsc
module load matlab

cd /scratch/$USER/IVopt-DDPC/src/

# Your commands or script for each job array
matlab -nosplash -nodesktop -r "seed_val=94; fprintf('seed = %d\n',seed_val); tic; main(N=1e5,seed=seed_val); toc, quit"
# srun matlab -batch "main(seed=10)"
echo "Finished MATLAB calculations"