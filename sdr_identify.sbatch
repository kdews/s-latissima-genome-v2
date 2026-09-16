#!/bin/bash
#SBATCH -J sdr_identify
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --time=01:00:00
#SBATCH --mem=86g
#SBATCH --cpus-per-task=12

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
ygs_url="https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/YGS.txt"
ygs_pl="YGS.pl"
scripts_dir="s-latissima-genome-v2" # scripts directory
if [[ -d "$scripts_dir" ]]; then
  rscript="$scripts_dir/$rscript"
  ygs_pl="$scripts_dir/$ygs_pl"
fi

# Help message
show_help () {
  echo "\
Runs k-mer based SDR identification from Barrera-Redondo et al. (2025).

Usage: sbatch $job_name.sbatch

Requires:
 - jellyfish v2.3.1 (https://www.htslib.org)"
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

# Download YGS.pl to scripts_dir (if needed)
[[ -f "$ygs_pl" ]] || wget -O "$ygs_pl" "$ygs_url"

# Input
in_fastq="example.fastq.gz"

# Output
outname="output"
out_jf="$outname.jf" # k-mer counts
out_jfa="$outname.fasta.gz" # k-mer FASTA (compressed)

# YGS method from Carvalho and Clark (2013) supplement
# https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/Supplemental_Material_Ykmer_ms_SM.pdf

# Load Jellyfish env
module purge
conda activate jellyfish
jellyfish --version || exit 1

# Generate 15-bp k-mer sequences from Illumina reads
count_cmd=(
  jellyfish count
  -t "$nthrd" # threads
  -m 15 # k-mer size = 15
  -s 10G # initial hash size
  -C # count both strand, canonical representation
  --quality-start=33 # ASCII for quality values
  --min-quality=20 # minimum quality (hard masked)
  -o "$out_jf" # output filename
)
echo "${count_cmd[*]} <($in_fastq)"
"${count_cmd[@]}" <("$in_fastq")

# Convert k-mer FASTQ to FASTA, filter for minimum frequency = 5
set -x
jellyfish dump --lower-count=5 | gzip -c > "$out_jfa"
set +x

