#!/bin/bash
#SBATCH -J pangenome_plots
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=5g
#SBATCH --time=05:00:00

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
rscript="$job_name.R" # R script
outdir="$job_name" # output directory
# Scripts directory
scripts_dir="s-latissima-genome-v2"
if [[ -d "$scripts_dir" ]]
then
  rscript="$scripts_dir/$rscript"
  outdir="$scripts_dir/figures/$outdir"
fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Generates alignment summary tables and plots, whole-genome andpairwise dotplots,
and a whole genome ribbon plot.

Usage: sbatch $job_name.sbatch

Output:
  $scripts_dir/
  └── figures/
      └── $(basename "$outdir")/
          ├── alignment_summary.tsv         # alignment summary table
          ├── chr_alignment_summary.tsv     # summary table by chromosome
          ├── chr_whole_genome_dotplot.png  # dotplot of chromosomes only
          ├── coverage_summary.png          # coverage summary plots
          ├── per_chromosome/
          │   ├── chr_<seqid>.png           # chromosome dotplots
          │   └── <sedid>.png               # contig dotplots
          ├── whole_genome_dotplot.png      # whole genome dotplot
          └── whole_genome_ribbon.png       # whole genome ribbon plot

Requires: R (https://www.r-project.org/)"
  exit 0
fi

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Input
# Merged PAF file
paf_file="cactus_pangenome_results/01-28-26_6039802/synteny_PAFs/merged_synteny.paf"
# The --reference assembly prefix (as it appears in PAF tname)
ref="corteva_v2"
# The second assembly name (as in PAF qname)
query="roscoff_v2"
# Cactus pangenome seqFile
seqFile="s-latissima-genome-v2/s_latissima_align1.txt"

# Load R modules
module purge
module load ver/2506 gcc/14.3.0 r/4.5.2
cmd=(
  Rscript
  "$rscript"
  "$paf_file"
  "$ref"
  "$query"
  "$seqFile"
  "$outdir"
)
echo "${cmd[*]}"
"${cmd[@]}"

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo
