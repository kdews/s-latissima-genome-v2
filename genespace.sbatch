#!/bin/bash
#SBATCH -J genespace
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p largemem
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=280g

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
rscript="$job_name.R" # R script to run
scripts_dir="s-latissima-genome-v2" # scripts directory
[[ -d "$scripts_dir" ]] && rscript="$scripts_dir/$rscript"

# Help message
show_help () {
  echo "\
Runs GENESPACE on genomes.

Usage: sbatch $job_name.sbatch <in_tsv>

Required input:
  in_tsv      TSV with 10 named columns (see below), with commented (#) header

Output:
  /         # 
  ├──       # 
  └──       # 

Requires:
 - R v6.1+ (https://www.r-project.org)
 - Bioconductor ()
 - 
 - GENESPACE (https://github.com/jtlovell/GENESPACE)
 - MCScanX (https://github.com/wyp1125/MCScanX)"
}
# Catch help flag or missing argument(s)
if [[ $1 == "-h" || $1 == "--help" ]]
then
  show_help
  exit 0
elif [[ "$#" -ne 1 ]]
then
  echo "Error: 1 positional argument required."
  echo
  show_help
  exit 1
fi

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
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"

# Load conda
module purge
cond="$HOME/.conda_for_sbatch.sh"
if [[ -f "$cond" ]]
then
  source "$cond"
else
  echo "Error sourcing $cond"
  exit 1
fi

# Multithread support
nthrd="${SLURM_CPUS_PER_TASK:-1}"

# Input
in_tsv="${1:-s-latissima-genome-v2/genomes_metadata.tsv}"
[[ -r $in_tsv ]] || die "Cannot read input TSV: $in_tsv"

# Load GENESPACE env
conda activate genespace
R --version | head -n1 || exit 1

# Run Rscript
set -x
Rscript "$rscript" "$in_tsv" "$scripts_dir"
set +x

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
