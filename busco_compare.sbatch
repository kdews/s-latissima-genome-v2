#!/bin/bash
#SBATCH -J busco_compare
#SBATCH -o %x_logs/%x_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=50g
#SBATCH --cpus-per-task=8

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
sum_dir="busco_summaries" # output directory
# BUSCO lineages
lineages=(
  "eukaryota"
  "stramenopiles"
)

# Help message
show_help () {
  echo "\
Generates plot(s) with conda-installed BUSCO generate_plot.py.

Usage: sbatch $job_name.sbatch <in_tsv>

Requires:
 - BUSCO (https://busco.ezlab.org)"
}
# Catch help flag
[[ $1 == "-h" || $1 == "--help" ]] && { show_help; exit 0; }

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
module purge
cond="$HOME/.conda_for_sbatch.sh"
if [[ -f "$cond" ]]
then
  source "$cond"
else
  echo "Error sourcing $cond"
  exit 1
fi

# BUSCO lineage + ortholog database version (_odb##)
# Activate BUSCO conda env
conda activate busco
odbs=()
for i in "${!lineages[@]}"
do
  odbs+=(
    "$(
      busco --list-datasets \
        | grep "${lineages[i]}" \
        | sed -r "s/.*(${lineages[i]}_odb.+) .*/\1/g"
    )"
  )
done
conda deactivate

# Input
# Assembly metadata
in_tsv="${1:-genomes_metadata.tsv}"
if [[ -f "$in_tsv" ]]
then
  mapfile -t labels < <(grep -v "#" "$in_tsv" | cut -f1)
else
  echo "Error: input TSV ($in_tsv) not found."
fi

# Compile BUSCO summary JSON files (if needed)
mkdir -p "$sum_dir"
for odb in "${odbs[@]}"
do
  wd="$sum_dir/$odb"
  mkdir -p "$wd"
  # Remove any previous JSON files from working directory (if they exist)
  [[ $(ls "$wd"/*.json 2>/dev/null) ]] && rm "$wd"/*.json
  # Copy BUSCO results for each lineage and assembly into summary directory
  for label in "${labels[@]}"
  do
    # busco_base="busco_${label}_$odb"
    busco_base="busco_results/$odb/$label"
    sum_file="$busco_base/short_summary.specific.$odb.$label.json"
    if [[ -f "$sum_file" ]]
    then
      cmd="cp -u $sum_file $wd/"
      echo "$cmd"
      $cmd
    else
      pwd
      echo "Error: $odb BUSCO summary file ($sum_file) not found for $label"
      exit 1
    fi
  done
  # Run BUSCO's native plotting function on each lineage
  # Activate BUSCO conda env
  conda activate busco
  cmd="busco --plot $wd"
  echo "$cmd"
  $cmd
  # Deactivate env
  conda deactivate
  # Generate combined figures of genome analyses
  # Load modules for Rscript
  # module load ver/2506 gcc/14.3.0 r/4.5.2
  module load r/4.6.0
  cmd=("Rscript" "$rscript" "$wd" "$in_tsv" "$scripts_dir")
  echo "${cmd[*]}"
  "${cmd[@]}"
  # Unload modules
  module purge
done

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
