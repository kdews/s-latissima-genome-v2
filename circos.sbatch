#!/bin/bash
#SBATCH -J circos
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p oneweek
#SBATCH --time=2-0
#SBATCH --mem=86g
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
Runs Circos using given configuration file(s).

Usage: sbatch $job_name.sbatch <label>

Requires: 
 - Circos (https://circos.ca/)"
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

# Load circos env
module purge
conda activate circos
circos --version || exit 1

# Define variables
# label=corteva_v2
label="$1" # assembly label
# Multithread support
nthrd="${SLURM_CPUS_PER_TASK:-1}"
# Input directory name
in="quast_results/$label/circos"
# Output directory name
out="circos_results/$label"
mkdir -p "$out"

# Copy default QUAST-LG Circos output to new directory
cp -r "$in/data" "$in/circos.conf" "$in/legend.txt" "$out/"

# Set conf file variables
conf="$out/circos.conf"
ticks="$out/data/ticks.conf"
legend="$out/legend.txt"
labels_file="$out/data/labels.txt"
ideogram="$out/data/ideogram.conf"

# Edit circos.conf
# Redirect to new Circos directory
sed -i "s|^dir = .*|dir = $(realpath -e "$out")|" "$conf"

# Edit data/ideogram.conf
# Turn on ideogram labels
sed -i "s/^show_label = no$/show_label = yes/" "$ideogram"
sed -i "s/^label_size = .*$/label_size = 14/" "$ideogram"
sed -i "s/^label_radius = .*$/label_radius = dims(ideogram,radius) + 0.12r/" "$ideogram"
sed -i "s/^label_font = .*$/label_font = condensed/" "$ideogram"

# Edit data/labels.txt with legend.txt
awk -F' - ' '/^[[:space:]]*assembly[0-9]+/ {
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
  print $1"\t"$2
}' "$legend" | while IFS=$'\t' read -r label name; do
  num=${label#assembly}
  track_idx=$((num - 1))

  # escape sed special chars (&, /, \) in the replacement name
  esc_name=$(printf '%s' "$name" | sed -e 's/[&/\]/\\&/g')

  sed -i "s/track${track_idx}=${label}/track${track_idx}=${esc_name}/" "$labels_file"
done

circos -conf "$out/circos.conf"

