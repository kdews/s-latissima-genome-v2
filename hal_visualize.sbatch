#!/bin/bash
#SBATCH -J hal_visualize
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=100gb
#SBATCH --time=10:00:00

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
date_str="$(date +"%m-%d-%y")" # filename date format
res_dir="${job_name}_results" # results directory for all jobs
outdir="$res_dir/${date_str}_$SLURM_JOB_ID" # output directory
workdir="$outdir/toil_workdir" # Toil working directory

# Help message
show_help () {
    echo "\
Summarizes and visualizes input HAL. Requires specification of reference
assembly label.

Usage: sbatch $job_name.sbatch <in_hal> <ref_label>

Requires:
 - virtualenv (https://virtualenv.pypa.io/en/latest)
 - Toil (https://toil.readthedocs.io/en/latest)
 - Cactus (https://github.com/ComparativeGenomicsToolkit/cactus)
   - contains HALtools distribution (http://github.com/glennhickey/hal)"
}
# Catch help flag or missing argument(s)
if [[ $1 == "-h" || $1 == "--help" ]]
then
  show_help
  exit 0
elif (( $# < 2 ))
then
  echo "Error: 2 positional arguments required."
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

# Input
in_hal="${1:-prog_cactus_results/06-10-26_9364670/s_latissima_prog_align.filt_1000000.hal}"
ref_label="${2:-corteva_v2}"

# Output
out_maf="$outdir/$(basename "$in_hal" .hal).maf"

# Load conda
module purge
cond="$HOME/.conda_for_sbatch.sh"
if [[ -e "$cond" ]]
then
  source "$cond"
else
  echo "Error sourcing $cond"
  exit 1
fi

# Load Python module and Cactus virtualenv
conda activate cactus_gpu
env_path="/home1/kdeweese/bin/cactus-bin-v3.2.1/venv-cactus-v3.2.1/bin/activate"
if [[ -e "$env_path" ]]
then
  source "$env_path"
else
  echo "Error on source of $env_path"
  exit 1
fi

# Verify which Cactus install is being sourced
echo "Cactus location: $(which cactus)"
halStats --help 2>&1 | grep "halStats v"

# Get basic information from alignment
set -x
halValidate "$in_hal"
halStats "$in_hal"
halSummarizeMutations "$in_hal"
set +x

# Create output directory
mkdir -p "$outdir"
# Create Toil working directory
mkdir -p "$workdir"

# Test cactus-hal2maf
echo "cactus-hal2maf v$(cactus-hal2maf --version)"
echo
# Convert HAL to MAF
set -x
cactus-hal2maf \
  --refGenome "$ref_label" \
  --workDir "$workdir" \
  "$outdir" \
  "$in_hal" \
  "$out_maf"
set +x

# Deactivate Cactus envs
echo
deactivate
conda deactivate

