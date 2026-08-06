#!/bin/bash
#SBATCH -J prog_cactus
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --time=2-0
#SBATCH --mem=10g

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
js="$outdir/jobstore" # Toil jobstore
workdir="$outdir/toil_workdir" # Toil working directory
logdir="$outdir/slurm_logs" # slurm log directory
# Allocate computing resources
partition="gpu" # slurm partition name
gpu_partition="gpu" # slurm GPU-enabled partition name
max_time="10:00:00" # maximum job length on partition
max_jobs=100 # maximum concurrent jobs on partition
max_mem=215 # maximum memory on node of partition (unit = Gi)
# max_mem=120 # maximum memory on node of partition (unit = Gi)
gpu_cores=2 # maximum GPU cores on node of partition
max_cores=32 # maximum cores on node of partition
def_mem="20Gi" # increase default memory from 2Gi
def_cores=2 # increase default cores from 1
def_disk="20Gi" # increase default disk memory from 2Gi

# Help message
show_help () {
  echo "\
Runs Progressive Cactus alignment. Requires as input a Cactus-formatted seqFile.

Usage: sbatch $job_name.sbatch <seqFile>

Required input:
  seqFile   Cactus-formatted seqFile containing a Newick tree and species list
            i.e.:   NEWICK tree
                    name1 path1
                    name2 path2
                    ...
                    nameN pathN

Output:
  $res_dir/
  └── ${date_str}_<job_id>/        # main Cactus output directory
      ├── $(basename "$js")/       # Toil jobstore (left by unsuccessful run)
      ├── $(basename "$workdir")/  # Toil working directory
      └── $(basename "$logdir")/   # extra slurm logs from Cactus and Toil

Requires:
 - virtualenv (https://virtualenv.pypa.io/en/latest)
 - Toil (https://toil.readthedocs.io/en/latest)
 - Cactus (https://github.com/ComparativeGenomicsToolkit/cactus)"
}
# Catch help flag or missing argument(s)
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  show_help
  exit 0
elif [[ "$#" -ne 1 ]]
then
  echo "Error: 1 positional argument expected."
  echo
  show_help
  exit 1
fi


# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"

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
# module load ver/2506 gcc/14.3.0 python/3.11.14
# env_path="/home1/kdeweese/bin/cactus-bin-v3.1.3/venv-cactus-v3.1.3/bin/activate"
if [[ -e "$env_path" ]]
then
  source "$env_path"
else
  echo "Error on source of $env_path"
  exit 1
fi

# Test Cactus and KegAlign
echo "Cactus location:"
which cactus || exit 1
printf "cactus v"
cactus --version
kegalign --version
echo

# Input
# Requires TSV listing assembly ID and FASTA path, i.e.,
# name1  path1
# name2  path2
# ...
# nameN  pathN
seqFile="${1:-s-latissima-genome-v2/s_latissima_prog_align.txt}"
if [[ -f "$seqFile" ]]
then
  echo "Using seqFile: $seqFile"
else
  echo "Error: no seqFile found ($seqFile)."
  exit 1
fi
# Prefix output files by seqFile
outname=$(basename "$seqFile" .txt) # output file prefix

# Output
hal="$outdir/$outname.hal"

# Create output directory
mkdir -p "$outdir"
# Create Toil working directory
mkdir -p "$workdir"

# Run cactus alignment 
# --gpu "$gpu_cores" \ enables GPU alignment acceleration with KegAlign
# --slurmGPUPartition "$gpu_partition" \
# --lastzCores "$gpu_cores" \
set -x
cactus \
  --batchSystem slurm \
  --workDir "$workdir" \
  --batchLogsDir "$logdir" \
  --slurmTime "$max_time" \
  --slurmPartition "$partition" \
  --defaultMemory "$def_mem" \
  --defaultCores "$def_cores" \
  --defaultDisk "$def_disk" \
  --maxJobs "$max_jobs" \
  --maxCores "$max_cores" \
  --maxMemory "${max_mem}Gi" \
  --consCores "$max_cores" \
  --consMemory "${max_mem}Gi" \
  --lastzMemory "$(( max_mem / 40 ))Gi" \
  --doubleMem true \
  "$js" "$seqFile" "$hal"
set +x

# Print date and time
echo
date +"$date_fmt"
