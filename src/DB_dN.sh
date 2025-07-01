#!/bin/bash

# Function to submit SLURM job array
submit_job_array() {
    local index=$1

    # Generate a unique SLURM submission script for each set of parameters
    script_name="./dN/iN${index}/submit_job_iN${index}.sh"

    # Create the SLURM submission script
    cat > "$script_name" <<EOL
#!/bin/bash
#SBATCH --job-name=iN${index}
#SBATCH --partition=compute
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900M
#SBATCH --account=research-me-dcsc
#SBATCH --output=./dN/iN${index}/job.%A_%a.out
#SBATCH --error=./dN/iN${index}/job.%A_%a.err
#SBATCH --array=1-100

# Load any necessary modules or set environment variables
module load matlab

# Use literal $ to defer expansion until job runtime
task_id=\$SLURM_ARRAY_TASK_ID
seed_val=\$(( (${index} - 1)*100 + task_id ))

matlab -nosplash -nodesktop -r "index=${index}; task_id=\${task_id}; seed_val=\${seed_val}; fprintf('seed = %d, index = %d, task_id = %d\n', seed_val, index, task_id); tic; main_dN(10,index,seed=seed_val); toc, quit"
echo "Finished MATLAB calculations"

EOL


    # Make the script executable
    chmod +x "$script_name"

    # Submit the SLURM job
    sbatch "$script_name"

    echo "Submitted job array ${index}"
}

# Specify the range of job arrays
start_index=1
end_index=10

cd ${HOME}/../../scratch/${USER}/IVopt-DDPC/src/

# Iterate over the range of job arrays
for index in $(seq $start_index $end_index); do
    # Call the function to submit the SLURM job array
    mkdir -p "./dN"
    mkdir -p "./dN/iN${index}"
    submit_job_array "$index"
done
echo "Submitted all Slurm job arrays"