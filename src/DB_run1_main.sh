#!/bin/bash
#SBATCH --job-name=dRe
#SBATCH --partition=compute
#SBATCH --time=00:03:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --account=research-me-dcsc
#SBATCH --output=./Job1/dRe.%A_%a.out
#SBATCH --error=./Job1/dRe.%A_%a.err

# Load any necessary modules or set environment variables
module load matlab

# Your commands or script for each job array
matlab -nosplash -nodesktop -r "seed_val=10; sprintf('seed = %d\n',seedval); main(seed_val), quit"
echo "Finished MATLAB calculations"