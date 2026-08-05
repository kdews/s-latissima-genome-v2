#!/bin/bash
#SBATCH -J ntSynt_viz
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
outdir="${job_name}_${date_str}_$SLURM_JOB_ID" # output directory
rscript="$job_name.R" # R script
# Scripts directory
scripts_dir="s-latissima-genome-v2"
[[ -d "$scripts_dir" ]] && rscript="$scripts_dir/$rscript"

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs ntSynt-viz on pangenome GFA result to visualize whole-genome ribbon plots.

Usage: sbatch $job_name.sbatch <seqFile> <cactus_dir>

Output:
  ${job_name}_<date>_<jobid>/
  ├── assembly_names.tsv                     # table of assemblies aligned
  ├── fais/                                  # FASTA indices of assemblies
  │   ├── corteva_v2.fai
  │   └── roscoff_v2.fai
  ├── fais.txt
  ├── name_conversion.tsv                    # table correlating assembly labels
  ├── <seqFile>.chrom-paint-feats.tsv        # coordinates for chrom. painting
  ├── <seqFile>.links.tsv                    # coordinates for links
  ├── <seqFile>.merged-blocks.tsv            # coordinates for merged blocks
  ├── <seqFile>.non-merge-blocks.tsv         # coordinates for non-merged blocks
  ├── <seqFile>.sequence_lengths.sorted.tsv  # sorted table of sequence lengths
  ├── <seqFile>.target_colours.tsv           # table of sequences and colors
  ├── <seqFile>_est-distances.nwk            # distances Newick tree
  ├── <seqFile>_est-distances.order.tsv      # table of ordered distances 
  ├── <seqFile>_est-distances.phylip         # distances Phylip
  ├── <seqFile>_est-distances.tsv            # table of distances
  └── <seqFile>_ribbon-plot.png              # pangenome ribbon plot
Requires: ntSynt-viz (https://github.com/BirolLab/ntSynt-viz/)"
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

# Load ntSynt-viz env
conda activate ntsynt-viz
printf "ntsynt-viz "
ntsynt_viz.py --version || exit 1
echo

# Download conversion Python scripts (if needed)
conv_link="https://raw.githubusercontent.com/BirolLab/ntSynt-viz/main/conversion_scripts/convert_gfa_to_ntsynt_blocks.py"
merge_link="https://raw.githubusercontent.com/BirolLab/ntSynt/main/bin/ntsynt_merge_collinear.py"
conv_py="$(basename "$conv_link")"
merge_py="$(basename "$merge_link")"
if [[ -d "$scripts_dir" ]]
then
  conv_py="${scripts_dir}${conv_py}"
  merge_py="${scripts_dir}${merge_py}"
fi
if [[ ! -f "$conv_py" ]]
then
  if [[ -d "$scripts_dir" ]]
  then
    wget -P "$scripts_dir" "$conv_link"
  else
    wget "$conv_link"
  fi
fi
if [[ ! -f "$merge_py" ]]
then
  if [[ -d "$scripts_dir" ]]
  then
    wget -P "$scripts_dir" "$merge_link"
  else
    wget "$merge_link"
  fi
fi
# Add Python scripts to PATH
chmod +x "$conv_py" "$merge_py"
if [[ -d "$scripts_dir" ]] 
then
  PATH="$PATH:$(realpath -e "$scripts_dir")"
else
  PATH="$PATH:$(pwd)"
fi
export PATH
# Test PATH
"$conv_py" --help > /dev/null || exit 1
"$merge_py" --help > /dev/null || exit 1

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
# Check input GFA exists
if [[ -f "$in_gfa" ]]
then
  # Handle gzipped input GFAs
  if [[ "$in_gfa" == *.gz ]]
  then
    in_gfa_d="${in_gfa%.gz}"
    [[ ! -f "$in_gfa_d" ]] && gunzip -k "$in_gfa"
    in_gfa="$in_gfa_d"
  fi
  echo "Cactus GFA: $in_gfa"
else
  echo "Error: Cactus GFA ($in_gfa) not found."
  exit 1
fi
echo

# Output
mkdir -p "$outdir" # create output directory
ass_names="$outdir/assembly_names.tsv" # assembly names TSV
fais="$outdir/fais.txt" # list of FAI files
n_conv="$outdir/name_conversion.tsv" # TSV to convert assembly names
out_prefix="$outdir/$seqFile_no_ext" # conversion output prefix

# Create "--assembly-names" TSV
cut -f1 "$seqFile" | awk '{print $1 "\t0\t" $1}' > "$ass_names"
# Create "--name_conversion" TSV
{
  cut -f1 "$seqFile" | \
  awk '
    {
      old=$0
      new=$0

      if ($0 ~ /corteva/) {
        new="US_S._latissima"
      } else if ($0 ~ /roscoff/) {
        new="France_S._latissima"
      } else if ($0 ~ /japonica/) {
        new="S._japonica"
      }

      print old "\t" new
    }'
} > "$n_conv"
# Target genome is first listed
target_genome="$(cut -f2 "$n_conv" | head -n1)" # need to use converted name

# Create directory for FAIs specified by "--fais" list file
mkdir -p "$outdir/fais"
# Copy and rename FAIs to outdir to match old prefixes
cat "$n_conv" | \
  while IFS=$'\t' read -r old_name _
  do
    fasta_path="$(awk -v asm="$old_name" '$1 == asm {print $2}' "$seqFile")"

    fai_src="${fasta_path}.fai"
    fai_dst="$outdir/fais/${old_name}.fai"

    cp "$fai_src" "$fai_dst"
    echo "$fai_dst"
  done > "$fais" # save new FAI filenames to "--fai" list file

# Convert cactus-pangenome GFA to ntSynt-format synteny block TSV
set -x
"$conv_py" \
  --merge 250000 \
  --gfa "$in_gfa" \
  --assembly-names "$ass_names" \
  --prefix "$out_prefix"
set +x

# Run ntSynt-viz
set -x
ntsynt_viz.py --normalize --format png --ribbon_adjust 0.14 --scale 1e9 \
  --target-genome "$target_genome" \
  --fais "$fais" \
  --name_conversion "$n_conv" \
  --blocks "$out_prefix.merged-blocks.tsv" \
  --prefix "$out_prefix"
set +x

# Run Rscript to improve ribbon plot
set -x
Rscript "$rscript" "$seqFile_no_ext" "$outdir" \
  --format png \
  --height 20 \
  --width 50 \
  --ratio 0.2 \
  --scale 100000000.0 \
  --no-arrow \
  -s "$outdir/$seqFile_no_ext.sequence_lengths.sorted.tsv" \
  -l "$outdir/$seqFile_no_ext.links.tsv" \
  -p "$outdir/${seqFile_no_ext}_ribbon-plot" \
  -c "$outdir/$seqFile_no_ext.chrom-paint-feats.tsv" \
  --colour_indices "$outdir/$seqFile_no_ext.target_colours.tsv"
set +x
# "--tree "$outdir/${seqFile_no_ext}_est-distances.nwk"

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo
