#!/bin/bash
#SBATCH --job-name=covidMRP
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null
#SBATCH --mem-per-cpu=8000
#SBATCH --array=1-7

module load R/4.4.1-gfbf-2023a

R CMD BATCH --no-save --no-restore SBC_posterior.R script_$SLURM_ARRAY_TASK_ID 