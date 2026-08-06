#!/bin/bash
#SBATCH -J cactus_pangenome
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --time=1-0
#SBATCH --mem=10g

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
res_dir="${job_name}_results" # results directory for all jobs
outdir="$res_dir/${date_str}_$SLURM_JOB_ID" # output directory
js="$outdir/jobstore" # Toil jobstore
workdir="$outdir/toil_workdir" # Toil working directory
logdir="$outdir/slurm_logs" # slurm log directory
# Allocate computing resources
partition="main" # slurm partition name
max_jobs=100 # maximum concurrent jobs on partition
max_time="05:00:00" # maximum job length on partition
max_mem=215 # maximum memory on node of partition (unit = G)
max_cores=32 # maximum cores on node of partition

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs Minigraph-Cactus pangenome alignment. Requires as input a Cactus-formatted
sequence file (seqFile) and a reference assembly ID (reference_ID).

Usage: sbatch $job_name.sbatch <seqFile> <reference_ID>

Required input:
  seqFile           TSV listing assembly ID and FASTA path, i.e.,
                      name1 path1
                      name2 path2
                      ...
                      nameN pathN
  reference_ID      assembly ID designated as reference

Output:
  ${job_name}_results/
      └── ${date_str}_<job_id>        # main Cactus output directory
          ├── $(basename "$js")/      # Toil jobstore (left by unsuccessful run)
          ├── $(basename "$workdir")/ # Toil working directory
          └── $(basename "$logdir")/  # extra Slurm logs from Cactus and Toil

Requires:
 - virtualenv (https://virtualenv.pypa.io/en/latest)
 - Toil (https://toil.readthedocs.io/en/latest)
 - Cactus (https://github.com/ComparativeGenomicsToolkit/cactus)"
  exit 0
fi

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"

# Load Python module and Cactus virtualenv
module purge
module load ver/2506 gcc/14.3.0 python/3.11.14
env_path="/home1/kdeweese/bin/cactus-bin-v3.1.3/venv-cactus-v3.1.3/bin/activate"
if [[ -e "$env_path" ]]
then
  source "$env_path"
else
  echo "Error on source of $env_path"
  exit 1
fi

# Test Cactus
which cactus || exit 1
printf "cactus v"
cactus --version
echo

# Input
# Requires TSV listing assembly ID and FASTA path, i.e.,
# name1  path1
# name2  path2
# ...
# nameN  pathN
if [[ "$#" -eq 2 ]]
then
  # seqFile=s_latissima_align.txt
  seqFile="$1"
  reference="$2"
  if [[ -f "$seqFile" ]]
  then
    echo "Using seqFile: $seqFile"
    if grep -q "$reference" "$seqFile"
    then
      echo "Using reference: $reference"
    else
      echo "Error: reference ($reference) not found in seqFile ($seqFile)."
      exit 1
    fi
  else
    echo "Error: no seqFile found ($seqFile)."
    exit 1
  fi
else
  echo "Error: 2 positional arguments expected."
  exit 1
fi
# Prefix output files by seqFile
outname="$(basename "$seqFile" .txt)" # output file prefix
ref_idx="$(grep "$reference" "$seqFile" | cut -f2).fai" # reference FASTA index

# Create output directory
mkdir -p "$outdir"
# Create Toil working directory
mkdir -p "$workdir"

# Run cactus-pangenome
cmd=(
  cactus-pangenome
  "$js" # must be accessible to all worker systems
  "$seqFile"
  --outDir "$outdir" # output directory
  --outName "$outname" # output file prefix
  --reference "$reference" # reference sample name
  --vcf # represents pangenome graph as sites of variation along reference
  --giraffe # make giraffe indexes for the filtered graph
  --gfa # standard text-based graph format
  --gbz # read-only format that scales well with # of paths (req. for giraffe)
  --chrom-og # make chromosomal graphs in odgi format
  --batchSystem slurm # enable slurm
  --workDir "$workdir" # toil workdir
  --batchLogsDir "$logdir" # directory for toil slurm logging
  --slurmTime "$max_time" # same for all slurm jobs, so pick highest needed
  --slurmPartition "$partition" # name of slurm partition to use
  --defaultMemory "20G" # increase default memory from 2Gi
  --defaultCores 2 # increase default cores from 1
  --defaultDisk "20G" # increase default disk memory from 2Gi
  --maxJobs "$max_jobs" # maximum number of jobs to submit to slurm at once
  --maxCores "$max_cores" # maximum cores schedulable on partition
  --maxMemory "${max_mem}G" # maximum memory schedulable on partition
  --mgCores $(( max_cores / 2 )) # cores for minigraph construction
  --mgMemory "${max_mem}G" # memory for minigraph construction
  --mapCores $(( max_cores / 2 )) # cores for each minigraph mapping job
  --consCores "$max_cores" # cores for each cactus_consolidated job (64 good)
  --consMemory "${max_mem}G" # max available when using slurm
  --indexCores $(( max_cores / 2 )) # half available when using slurm
  --indexMemory "${max_mem}G" # max available when using slurm
  --doubleMem true # if OOM, retry asking for double the memory
  --permissiveContigFilter # don't filter out small contigs
)
# Specify file with chr. list for ref assembly (generated by scaffold_eval.R)
ref_contigs_file="$reference-chromosome_ids.txt"
cut -f1 "$ref_idx" | head -n31 > "$ref_contigs_file"
if [[ -f "$ref_contigs_file" ]]
then
  echo "Chromosomes in reference: $ref_contigs_file ($(wc -l < "$ref_contigs_file"))"
  mapfile -t ref_contigs < "$ref_contigs_file"
  cmd+=(
    --refContigs "${ref_contigs[@]}"
    --otherContig chrOther
  )
fi
echo "${cmd[*]}"
"${cmd[@]}"

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
