#!/bin/bash
#SBATCH -J yaml_to_metadata
#SBATCH -o %x_logs/%x_%j.log
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5g

# Script variables
# Extract job name from script
if [[ -z "${SLURM_JOB_NAME}" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
rscript="$job_name.R" # R script to run
yaml_dir="assemblies_metadata" # input directory containing genome metadata YAMLs
out_tsv="genomes_metadata.tsv" # output TSV
scripts_dir="s-latissima-genome-v2" # scripts directory
if [[ -d "$scripts_dir" ]]; then
  rscript="$scripts_dir/$rscript"
  yaml_dir="$scripts_dir/$yaml_dir"
  out_tsv="$scripts_dir/$out_tsv"
fi

# Help message
show_help () {
  echo "\
Convert one YAML file per assembly into a 10-column metadata TSV.

Usage: sbatch $job_name.sbatch

Input: 
  $yaml_dir/             # input directory containing genome metadata YAMLs
  ├── <label1>.yaml
  └── <label2>.yaml 

Output: $out_tsv

Requires:
 - R v6+ (https://www.r-project.org)"
}
# Catch help flag
[[ $1 == "-h" || $1 == "--help" ]] && { show_help; exit 0; }

# Catch errors in pipes
set -o pipefail

# Print Slurm job information
[[ -n "$SLURM_JOB_ID" ]] && echo "\
==========================================
SLURM_JOB_ID = $SLURM_JOB_ID
SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST
TMPDIR = $TMPDIR
=========================================="

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Load modules for Rscript
module purge
module load r/4.6.0
R --version | head -n1 || exit 1

# Run Rscript
set -x
Rscript "$rscript" "$yaml_dir" "$out_tsv"
set +x

# Print date and time
echo
date +"$date_fmt"
