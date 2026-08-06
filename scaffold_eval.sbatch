#!/bin/bash
#SBATCH -J scaffold_eval
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=8g
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
rscript="$job_name.R" # Rscript
scripts_dir="s-latissima-genome-v2" # scripts directory
[[ -d "$scripts_dir" ]] && rscript="$scripts_dir/$rscript"

# Help message
show_help () {
  echo "\
Given a TSV manifest of genome metadata, parses list of assemblies parsed and 
generates plots comparing assemblies on the basis of scaffolding.

Usage: sbatch $job_name.sbatch <in_tsv>

Required input:
  in_tsv    table of genome metadata, with the required columns:
    # Label  Pretty_Label  ...

Requires:
 - samtools (https://www.htslib.org)
 - R v6.0+ (https://www.r-project.org)"
}
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
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

set -o pipefail

# Print Slurm job information
[[ -n "$SLURM_JOB_ID" ]] && echo "\
==========================================
SLURM_JOB_ID = $SLURM_JOB_ID
SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST
TMPDIR = $TMPDIR
=========================================="

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

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Multithread support
nthrd="${SLURM_CPUS_PER_TASK:-1}"

# Input
in_tsv="${1:-genomes_metadata.tsv}"
if [[ -f "$in_tsv" ]]
then
  echo "Using TSV: $in_tsv"
  # Assembly labels
  mapfile -t labels < <(grep -v "#" "$in_tsv" | cut -f1)
  mapfile -t pretty_labels < <(grep -v "#" "$in_tsv" | cut -f2)
  echo "Found ${#labels[@]} assembly labels."
else
  echo "Error: input TSV ($in_tsv) not found."
  exit 1
fi
echo

# Load samtools conda env
conda activate samtools
samtools --version | head -n2 || exit 1
echo
# Run samtools index on each FASTA
for i in "${!labels[@]}"
do
  label="${labels[i]}"
  pretty_label="${pretty_labels[i]}"
  fasta="assemblies/$label.fa"
  echo "$pretty_label"
  echo "FASTA: $fasta"
  [[ -f "$fasta" ]] || { echo "Error: FASTA not found." ; exit 1 ; }
  if [[ -f "$fasta.fai" ]]
  then
    # Skip indexing if exists
    echo "Found file $fasta.fai, skipping samtools faidx step."
  else
    # Index assembly
    set -x
    samtools faidx -@ "$nthrd" "$fasta"
    set +x
  fi
  echo
done

# Plot and filter assembly indices
# Load R module
module purge
module load r/4.6.0
Rscript --version
# Run Rscript to plot and filter contigs
set -x
Rscript "$rscript" "$in_tsv" "$scripts_dir"
set +x
