#!/bin/bash
#SBATCH -J depth_sdr
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p oneweek
#SBATCH --time=7-0
#SBATCH --mem=60g
#SBATCH --cpus-per-task=16

# Script variables
# Extract job name from script
if [[ -n ${SLURM_JOB_NAME:-} ]]; then
  job_name="$SLURM_JOB_NAME"
elif [[ -n ${BASH_SOURCE[0]:-} && ${BASH_SOURCE[0]} != "bash" ]]; then
  job_name="$(basename "${BASH_SOURCE[0]}")"
  job_name="${job_name%%.*}"
else
  job_name="depth_sdr"
fi
