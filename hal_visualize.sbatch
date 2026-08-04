#!/bin/bash
#SBATCH -J hal_visualize
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=100gb
#SBATCH --time=10:00:00

# Print Slurm job information
echo "\
==========================================
SLURM_JOB_ID = $SLURM_JOB_ID
SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST
TMPDIR = $TMPDIR
=========================================="

# Input
seqFile="s_latissima_prog_align.txt"
in_hal="prog_cactus_06-10-26_9364670/s_latissima_prog_align.filt_1000000.hal"
ref="corteva_v2"

# Output
workdir="$outdir/toil_workdir" # Toil working directory

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

# Create Toil working directory
mkdir -p "$workdir"

# Test cactus-hal2maf
echo "cactus-hal2maf v$(cactus-hal2maf --version)"
echo
# Convert HAL to MAF
set -x
cactus-hal2maf --refGenome "$ref" \
  "$in_hal"
set +x

# Deactivate Cactus envs
echo
deactivate
conda deactivate

