#!/bin/bash
#SBATCH -J make_guide_tree
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p largemem
#SBATCH --time=05:00:00
#SBATCH --mem=180g
#SBATCH --cpus-per-task=16

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
scripts_dir="s-latissima-genome-v2" # scripts directory
[[ -d "$scripts_dir" ]] && rscript="$scripts_dir/$rscript"
lineage="stramenopiles_odb12.2" # BUSCO lineage name

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"

# Input
# Requires TSV listing assembly ID and FASTA path, i.e.,
# name1  path1
# name2  path2
# ...
# nameN  pathN
if [[ "$#" -eq 1 ]]
then
  # seqFile=s_latissima_prog_align.txt
  seqFile="$1"
  if [[ -f "$seqFile" ]]
  then
    echo "Using seqFile: $seqFile"
  else
    echo "Error: no seqFile found ($seqFile)."
    exit 1
  fi
else
  echo "Error: 1 positional argument expected."
  exit 1
fi

# Load R modules
module purge
# module load apptainer/1.4.5 rstats/4.5.2
module load r/4.6.0
Rscript --version || exit 1
echo

# Run Rscript
set -x
Rscript "$rscript" "$seqFile" "$lineage" "$scripts_dir"
set +x

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
