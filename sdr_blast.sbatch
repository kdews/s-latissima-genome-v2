#!/bin/bash
#SBATCH -J sdr_blast
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --time=05:00:00
#SBATCH --mem=15gb
#SBATCH --cpus-per-task=12

# Script variables
# Extract job name from script
if [[ -z "${SLURM_JOB_NAME}" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
date_str="$(date +"%m-%d-%y")" # filename date format
outdir="${job_name}_${date_str}_$SLURM_JOB_ID" # output directory
rscript="$job_name.R" # R script
# Scripts directory
scripts_dir="s-latissima-genome-v2/"
[[ -d "$scripts_dir" ]] && rscript="${scripts_dir}${rscript}"

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Builds a reference protein set from the known Ectocarpus SDR genes to use as a 
database for reciprocal BLAST.

Usage: sbatch $job_name.sbatch <seqFile> <cactus_dir>

Output:
  ${job_name}_<date>_<jobid>/
  ├── 
  ├── 
  └── 
Requires:
  - gffread (https://github.com/gpertea/gffread)
  - BLAST ()"
  exit 0
fi

# -e  exit on error
# -o pipefail  catch errors in pipes
set -eo pipefail

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Multithread support
nthrd="${SLURM_CPUS_PER_TASK:-1}"

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

# Load gffread env
conda activate gffread
printf "gffread "
gffread --version || exit 1
echo

# Input
GENOME="Ec32_V5.reference.fasta"
GFF="Ec32_V5.reference.gff"
CONTIG="Ec25_SDR_F"
# Output
OUTDIR="Ec25_SDR_F_extract"
PREFIX="Ec25_SDR_F"

mkdir -p "$OUTDIR"

# 1) Pull out only the SDR contig from the genome
samtools faidx "$GENOME" "$CONTIG" > "$OUTDIR/${PREFIX}.fasta"

# 2) Pull out only GFF3 records on that contig
awk -v c="$CONTIG" '
  /^#/ {print; next}
  $1 == c {print}
' "$GFF" > "$OUTDIR/${PREFIX}.gff"

# 3) Extract CDS, transcripts, and proteins from the contig-specific annotation
gffread "$OUTDIR/${PREFIX}.gff" \
  -g "$OUTDIR/${PREFIX}.fasta" \
  -x "$OUTDIR/${PREFIX}.cds.fa" \
  -y "$OUTDIR/${PREFIX}.proteins.fa" \
  -w "$OUTDIR/${PREFIX}.transcripts.fa"

