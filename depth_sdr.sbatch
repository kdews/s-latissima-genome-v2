#!/bin/bash
#SBATCH -J ygs_sdr
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p oneweek
#SBATCH --time=7-0
#SBATCH --mem=60g
#SBATCH --cpus-per-task=16

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]; then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME" # job_name=sdr_id
fi

