#!/bin/bash
#SBATCH -J odgi_viz
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=240g
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=32

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
rscript="$job_name.R" # R script
# Scripts directory
scripts_dir="s-latissima-genome-v2/"
[[ -d "$scripts_dir" ]] && rscript="${scripts_dir}${rscript}"

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Generates whole-pangenome and chromosome-scale graph visualizations from a
Cactus pangenome GFA using ODGI. Produces sorted ODGI graphs, graph statistics,
path lists, and overview PNGs for inspecting structural variation and genome
path relationships across assemblies.

Usage: sbatch $job_name.sbatch <seqFile> <cactus_dir>

Output:
  ${job_name}_results/
      └── ${date_str}_<job_id>/
          ├── graph.raw.og                  # ODGI graph converted from GFA
          ├── graph.sorted.og               # path-guided sorted ODGI graph
          ├── graph.stats.txt               # graph statistics summary
          ├── graph.paths.txt               # list of graph paths/genomes
          ├── graph.overview.compressed.png # compressed whole-pangenome overview
          ├── graph.overview.paths.png      # whole-pangenome path visualization
          └── <seqid>/                      # per-chromosome (ref) visualizations
            ├── <seqid>.full.sorted.og
            ├── <seqid>.full.stats.txt
            ├── <seqid>.full.lay
            ├── <seqid>.full.draw.png
            └── <seqid>.full.viz.png
  
Requires: ODGI (https://github.com/pangenome/odgi)"
  exit 0
fi

# -e  exit on error
# -o pipefail  catch errors in pipes
set -eo pipefail

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Multithread support
nthrd="${SLURM_CPUS_PER_TASK:-1}"

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

# Load odgi env
conda activate odgi
printf "odgi "
odgi version || exit 1
echo

# Input
if [[ "$#" -eq 2 ]]
then
  # seqFile="s_latissima_align.txt"
  seqFile="$1"
  # cactus_dir="cactus_pangenome_05-21-26_8926505"
  cactus_dir="$2"
  if [[ -f "$seqFile" ]]
  then
    echo "Using seqFile: $seqFile"
    echo
    if [[ -d "$cactus_dir" ]]
    then
      echo "Using cactus result directory: $cactus_dir"
    else
      echo "Error: cactus result directory ($cactus_dir) not found."
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
echo
seqFile_no_ext="$(basename "$seqFile" .txt)" # cactus prefix
in_gfa="$cactus_dir/s_latissima_align.gfa.gz" # cactus-pangenome GFA
in_chr_dir="$cactus_dir/s_latissima_align.chroms" # odgi .chroms directory
# Check input files exist
if [[ -f "$in_gfa" ]]
then
  echo "Cactus GFA: $in_gfa"
else
  echo "Error: Cactus GFA ($in_gfa) not found."
  exit 1
fi
echo
if [[ -d "$in_chr_dir" ]]
then
  echo "odgi chroms directory: $in_chr_dir"
else
  echo "Error: odgi chroms directory ($in_chr_dir) not found."
  exit 1
fi
echo

# Create output directory
mkdir -p "$outdir"

# Per-chromosome visualizations

# Input per-chromosome odgi files
mapfile -t og_files < <(find "$in_chr_dir" -type f -name "*.full.og")

# Loop through each odgi file
for og in "${og_files[@]}"
do
  # Check input odgi exists
  if [[ -f "$og" ]]
  then
    echo "Input odgi: $og"
  else
    echo "Error: Input odgi ($og) not found."
    exit 1
  fi

  # Output
  chr="$(basename "$og" .full.og)" # parse chromosome name
  chr_out="$outdir/$chr" # per-chromosome output directory
  mkdir -p "$chr_out" # create output directory
  sorted="$chr_out/${chr}.full.sorted.og" # sorted odgi
  stats="$chr_out/${chr}.full.stats.txt" # graph summary
  paths="$chr_out/${chr}.full.paths.txt" # path names in graph
  colors="$chr_out/${chr}.full.colors.txt" # path colors in graph
  viz="$chr_out/${chr}.full.viz.png" # output graph
  lay="$chr_out/${chr}.full.lay" # layout
  draw="$chr_out/${chr}.full.draw.png" # output graph

  # Sort, summarize, and vizualize each chromosome pangenome graph
  set -x
  odgi sort -t "$nthrd" -p Ygs -i "$og" -o "$sorted"
  odgi stats -S -i "$sorted" > "$stats"
  # Remove contig/unscaffolded alignments
  odgi paths -L -i "$sorted" | grep -vP "contig|#un" > "$paths"
  cat "$paths" | awk '
    BEGIN {
        colors["corteva_v2_polished"]="#0072B2";
        colors["roscoff_v2"]="#CC79A7";
        colors["japonica_li2025"]="#009E73";
    }
    {
        split($0, a, "#");
        asm=a[1];
        print $0 "\t" colors[asm];
    }' > "$colors"
  odgi viz -t "$nthrd" -C -a 20 -X 4 \
    -i "$sorted" -p "$paths" -F "$colors" -o "$viz"
  set +x
  # odgi layout -t "$nthrd" -i "$sorted" -o "$lay"
  # odgi draw -t "$nthrd" -C -i "$sorted" -c "$lay" -p "$draw"
  echo
done

# # Whole-genome visualization

# # Input
# # Handle gzipped input GFAs
# if [[ "$in_gfa" == *.gz ]]
# then
#   in_gfa_d="${in_gfa%.gz}"
#   [[ ! -f "$in_gfa_d" ]] && gunzip -k "$in_gfa"
#   in_gfa="$in_gfa_d"
# fi

# # Output
# raw_og="$outdir/graph.raw.og"
# sorted_og="$outdir/graph.sorted.og"
# stats_txt="$outdir/graph.stats.txt"
# paths_txt="$outdir/graph.paths.txt"
# # Whole-pangenome overview figures
# ov_comp="$outdir/graph.overview.compressed.png"
# ov_paths="$outdir/graph.overview.paths.png"

# # Build an ODGI graph from the GFA.
# # odgi build supports GFAv1 input and writes an .og file.
# set -x
# odgi build -t "$nthrd" -O -g "$in_gfa" -o "$raw_og"
# set +x
# echo

# # Sort the graph for a cleaner 1D layout.
# # gs pipeline = groom, then topological sort.
# set -x
# odgi sort -t "$nthrd" -p gs -i "$raw_og" -o "$sorted_og"
# set +x
# echo

# # Basic graph summary.
# set -x
# odgi stats -S -i "$sorted_og" > "$stats_txt"
# set +x
# echo

# # Write path names to a text file for later subsetting.
# set -x
# odgi paths -L -i "$sorted_og" > "$paths_txt"
# set +x
# echo

# # Whole-chromosome / whole-pangenome overview:
# # compressed mode summarizes coverage across paths in one image.
# set -x
# odgi viz -t "$nthrd" -O -x 6000 -y 1800 -i "$sorted_og" -o "$ov_comp"
# set +x
# echo

# # A second overview with all paths shown as separate rows.
# # For very path-dense graphs, you can change this to a subset file with -p.
# set -x
# odgi viz -t "$nthrd" -k -x 6000 -y 3000 -i "$sorted_og" -o "$ov_paths"
# set +x
# echo

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo
