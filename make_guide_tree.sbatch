#!/bin/bash
#SBATCH -J make_guide_tree
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p largemem
#SBATCH --time=05:00:00
#SBATCH --mem=180g
#SBATCH --cpus-per-task=16

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
rscript="$job_name.R" # R script to run
scripts_dir="s-latissima-genome-v2" # scripts directory
[[ -d "$scripts_dir" ]] && rscript="$scripts_dir/$rscript"
lineage="stramenopiles_odb12.2" # BUSCO lineage name

# Help message
show_help () {
  echo "\
Creates maximum-likelihood (ML) guide tree for Progressive Cactus alignment
by aligning BUSCO-predicted proteins for multiple assemblies using the lineage
$lineage. Requires as input a Cactus-formatted seqFile.

Usage: sbatch $job_name.sbatch <seqFile>

Required input:
  seqFile   Cactus-formatted seqFile containing assembly labels and paths
            i.e.:   label1 path1
                    label2 path2
                    ...
                    labelN pathN

Output:
  <seqFile>    # seqFile updated with Newick tree on first line
  $scripts_dir/
  ├── figures/
  │   ├── <seqFile>_unique_og_dist.png
  │   ├── <seqFile>_NJ_ML_trees.png
  │   └── <seqFile>_msa_ML_tree.png
  ├── <seqFile>_NJ_tree.txt
  └── <seqFile>_ML_tree.txt

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

# Load R modules
module purge
# module load apptainer/1.4.5 rstats/4.5.2
module load r/4.6.0
Rscript --version || exit 1
echo

# Run Rscript
set -x
Rscript "$rscript" "$seqFile" "$lineage" "$scripts_dir"
set +x

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
