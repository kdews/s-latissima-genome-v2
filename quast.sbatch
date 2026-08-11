#!/bin/bash
#SBATCH -J quast
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p gpu
#SBATCH --time=18:00:00
#SBATCH --mem=180g
#SBATCH --cpus-per-task=12

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]; then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi

# Help message
show_help () {
  echo "\
Runs QUAST on a given assembly 

Usage: sbatch $job_name.sbatch <assembly.fasta> <label>

Requires: 
 - QUAST (https://quast.sourceforge.net)"
}
# Catch help flag or missing argument(s)
if [[ $1 == "-h" || $1 == "--help" ]]; then
  show_help
  exit 0
elif (( $# < 3 )); then
  echo "Error: at least 3 positional arguments required."
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

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Load conda
cond="$HOME/.conda_for_sbatch.sh"
if [[ -f "$cond" ]]; then
  source "$cond"
else
  echo "Error sourcing $cond"
  exit 1
fi

# Load QUAST env
module purge
conda activate quast
quast --version || exit 1
circos --version || exit 1

# Define variables
label="$1" # assembly label
ref="$2" # path to reference FASTA (chromosome-only preferred)
reads_opt="$3" # variable to toggle reads mode on/off
assemblies=("${@:4}") # uncompressed assembly FASTA(s)
# Multithread support
nthrd="${SLURM_CPUS_PER_TASK:-1}"
# Output directory name
out="quast_results/$label"
# Log details
echo
echo "Assembly label: $label"
echo "Assembly FASTA(s): ${assemblies[*]}"
# Check input FASTA(s) exist
for assembly in "${assemblies[@]}"; do
  [[ -f "$assembly" ]] || { echo "Error: Assembly FASTA not found."; exit 1; }
done
[[ -f "$ref" ]] || { echo "Error: Reference FASTA not found."; exit 1; }
echo

# QUAST command
# quast="quast --large -k -b --report-all-metrics --threads $nthrd -o $out"
cmd=(
  quast
  --large # enable large eukaryotic genome settings
  --circos # make Circos plot
  # -s # split assembly scaffolds into contigs for analysis ("<label>_broken")
  -k # use k-mers to determine assembly completion (default k=101bp)
  # -b # run BUSCO on assemblies
  --glimmer # gene finding with GlimmerHMM
  # --rna-finding # ribosomal RNA gene finding
  # -g "$gff" # feature file
  --threads "$nthrd"
  -r "$ref" # reference genome
  -o "$out"
)

# Check for assembly reads
# Find FASTQ(s) named by label stripped of version number
if [[ "$reads_opt" == reads ]]; then
  clean_label="${label%%_v[0-9]*}"
  pe1="reads/${clean_label}_1.fastq.gz"
  pe2="reads/${clean_label}_2.fastq.gz"
  long_pb="reads/${clean_label}_pacbio.fastq.gz"
  long_ont="reads/${clean_label}_nanopore.fastq.gz"
  # Append reads flags to command
  if [[ -f "$pe1" && -f "$pe2" ]]; then
    echo "Using PE reads."
    echo "PE1: $pe1"
    echo "PE2: $pe2"
    cmd+=(--pe1 "$pe1" --pe2 "$pe2")
  elif [[ -f "$long_pb" ]]; then
    echo "Using PacBio reads: $long_pb"
    cmd+=(--pacbio "$long_pb")
  elif [[ -f "$long_ont" ]]; then
    echo "Using Nanopore reads: $long_ont"
    cmd+=(--nanopore "$long_ont")
  fi
fi
# Append assembly path(s) to command
cmd+=("${assemblies[@]}")

# # Set reference for QUAST
# last_idx="$(( ${#assemblies[@]} - 1 ))"
# if [[ "${#assemblies[@]}" -eq 1 ]]; then
#   # Use default reference for single assemblies
#   echo "Reference: $def_ref"
#   cmd+=(-r "$def_ref")
#   cmd+=("${assemblies[0]}")
# else
#   # Use latest version as reference for multiple assembly versions
#   ref="${assemblies[last_idx]}"
#   echo "Reference: $ref"
#   cmd+=(-r "$ref")
#   cmd+=("${assemblies[@]:0:$last_idx}")
# fi

# Run QUAST
echo "${cmd[*]}"
"${cmd[@]}"

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
