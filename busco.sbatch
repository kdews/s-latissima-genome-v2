#!/bin/bash
#SBATCH -J busco
#SBATCH -o %x_logs/%x_%j.log
#SBATCH --time=10:00:00
#SBATCH --mem=15g
#SBATCH --nodes=2
#SBATCH --cpus-per-task=12
#SBATCH --ntasks-per-node=1

# Script variables
# Extract job name from script
if [[ -z "${SLURM_JOB_NAME}" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi

echo

# Help message
show_help () {
  echo "\
Runs BUSCO on given assembly using the latest versions of the 
Eukaryota and Stramenopiles lineage ortholog databases (ODBs).

Usage: sbatch $job_name.sbatch <assembly.fasta> <label>

Output:
  busco_results/
  ├── eukaryota_odb<latest>/
  |   └── <label>/
  └── stramenopiles_odb<latest>/
      └── <label>/

Requires: 
- BUSCO (https://busco.ezlab.org)
- GNU parallel (https://www.gnu.org/software/parallel)"
}
# Catch help flag or missing argument(s)
if [[ $1 == "-h" || $1 == "--help" ]]
then
  show_help
  exit 0
elif [[ "$#" -ne 2 ]]
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

# Load GNU Parallel
module purge
module load ver/2506 gcc/14.3.0 parallel/20250522
parallel --version || exit 1

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
conda activate busco
busco --version || exit 1

# Define variables
# BUSCO mode
mode="genome"
# BUSCO lineages
lineages=(
  "eukaryota"
  "stramenopiles"
)
# BUSCO lineage + ortholog database version (_odb##)
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
# Define variables
label="$1" # assembly label
assembly="$2" # uncompressed assembly FASTA
# Number of parallel tasks (should = number of ODBs)
ntsk="${SLURM_NTASKS:-2}"
# Multithread support
nthrd="${SLURM_CPUS_PER_TASK:-1}"
# Joblog name
joblog="${job_name}_logs/${job_name}_$label.joblog"
# Check input exists and log details
echo
echo "Assembly label: $label"
echo "Assembly FASTA: $assembly"
[[ -f "$assembly" ]] || { echo "Error - assembly FASTA not found."; exit 1; }
echo

# Increase java heap memory (for assemblies >650Mb)
java_mem="$(( 9 * SLURM_MEM_PER_NODE / 10 ))M"
export _JAVA_OPTIONS="-Xmx$java_mem"
export JAVA_TOOL_OPTIONS="-Xmx$java_mem"

# Invoke parallel with jobs=$SLURM_NTASKS and a joblog file
# srun gives 1 core to the set of threads defined by $SLURM_CPUS_PER_TASK
cmd=(
  parallel
  -v # print command before running
  -j "$ntsk" # number of jobs to run (= number of lineages)
  --joblog "$joblog" # logs status of each job after completion
  srun
  --export=ALL # export environment variables to srun
  --exclusive # each job uses its own node and resources
  -N 1 # number of nodes per job
  -n 1 # number of tasks per job
  -c "$nthrd" # threads per job
  busco
  -f # force overwrite existing output files
  # --metaeuk # use MetaEuk gene predictor
  # --augustus # use AUGUSTUS gene predictor
  --scaffold_composition # writes base content per scaffold to scaffold_composition.txt
  -c "$nthrd" # threads for BUSCO (≤ CPUS_PER_TASK)
  -m "$mode"
  -i "$assembly"
  -l {} # lineage (replaced by input from parallel)
  --out_path "${job_name}_results/{}"
  -o "$label"
  :::
  "${odbs[@]}"
)
echo "${cmd[*]}"
"${cmd[@]}"

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
