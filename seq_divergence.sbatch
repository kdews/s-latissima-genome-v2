#!/bin/bash
#SBATCH -J seq_divergence
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=15gb
#SBATCH --time=05:00:00

# Script variables
# Extract job name from script
if [[ -z "${SLURM_JOB_NAME}" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
scripts_dir="s-latissima-genome-v2/" # scripts directory
outdir="$job_name" # output directory

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Extracts individual and merged PAF files from chromosome alignments.
Requires output directory from cactus-pangenome (cactusDir) as input.
Usage: sbatch $job_name.sbatch

Required input:
  minlen   minimum length to filter contigs out below

Output:
  └── $outdir
      ├── $(basename "$out_md5")      expected md5 of FASTQs parsed ENA TSV
      └── <ENA_ID>.fastq.gz         locally downloaded FASTQs

Requires:
 - R v4.6.0 (https://www.r-project.org/)"
  exit 0
fi

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"

# Create output directory
mkdir -p "$outdir"

# Print date and time
echo
date +"$date_fmt"

