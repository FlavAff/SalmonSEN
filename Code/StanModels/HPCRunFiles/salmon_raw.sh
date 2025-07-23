#!/bin/bash
#
# Slurm job script for running the salmon model analysis in R
#
#SBATCH --job-name=salmon_models
#SBATCH --output=salmon_models_%j.out
#SBATCH --error=salmon_models_%j.err
#SBATCH --account=rpp-gonzalez
#SBATCH --time=72:00:00         # Set a maximum runtime
#SBATCH --ntasks=1              # We'll run the R script in a single task
#SBATCH --cpus-per-task=8       # Request number of CPUs for cmdstanr 
#SBATCH --mem-per-cpu=12GB              # Request memory
#SBATCH --mail-type=BEGIN,END,FAIL    # Send email on job end or failure
#SBATCH --mail-user=flavio.affinito@mail.mcgill.ca # Your email address

# Load necessary modules (adjust these based on your cluster's environment)
module load gcc/12.3 r/4.3.1  # Or your preferred R version

# Run the R script
Rscript run_all_models_raw.R