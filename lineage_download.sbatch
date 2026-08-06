#!/bin/bash
#SBATCH -J lineage_download
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=4g

# Script variables
# Extract job name from script
if [[ -z "${SLURM_JOB_NAME}" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi

# Help message
show_help () {
  echo "\
Downloads lastest BUSCO ortholog databases (ODBs) for  Eukaryota an
Stramenopiles lineages.

Usage: sbatch $job_name.sbatch

Output:
  busco_downloads/
  ├── information/
  └── lineages/
      ├── eukaryota_odb<latest>/
      │   ├── hmms/
      │   ├── info/
      │   └── prfl/
      └── stramenopiles_odb<latest>/
          ├── hmms/
          ├── info/
          └── prfl/

Requires: 
- BUSCO (https://busco.ezlab.org)
- GNU parallel (https://www.gnu.org/software/parallel)"
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
cond="$HOME/.conda_for_sbatch.sh"
if [[ -f "$cond" ]]
then
  source "$cond"
else
  echo "Error sourcing $cond"
  exit 1
fi

# Load BUSCO env
module purge
conda activate busco
busco --version || exit 1

# Download required BUSCO lineage datasets
lineages=(
  "eukaryota"
  "stramenopiles"
)
for i in "${!lineages[@]}"
do
  odb="$(
    busco --list-datasets \
      | grep "${lineages[i]}" \
      | sed -r "s/.*(${lineages[i]}_odb.+) .*/\1/g"
  )"
  busco --download "$odb"
done
