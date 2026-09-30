#!/bin/bash
#SBATCH -J sdr_id
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p oneweek
#SBATCH --time=01:00:00
#SBATCH --mem=10g
#SBATCH --cpus-per-task=16

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]; then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME" # job_name=sdr_id
fi
# URL to YGS.pl from Carvalho and Clark (2013)
ygs_url="https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/YGS.txt"
ygs_pl="YGS.pl" # output name of YGS script once downloaded
sdr_parse="sdr_parse_multiqc.R" # R script parses reads for YGS from MultiQC data
ygs_job="ygs_sdr.sbatch" # YGS SDR ID job submission script
depth_job="depth_sdr.sbatch" # read depth SDR ID job submission script
scripts_dir="s-latissima-genome-v2" # scripts directory
if [[ -d "$scripts_dir" ]]; then
  sdr_parse="$scripts_dir/$sdr_parse"
  ygs_pl="$scripts_dir/$ygs_pl"
  ygs_job="$scripts_dir/$ygs_job"
  depth_job="$scripts_dir/$depth_job"
fi
ygs_pl="$(realpath -e "$ygs_pl")"
outdir="${job_name}_results" # main output directory
window=500000 # window size = 500 kb
win_str="$((window / 1000))kb" # get string from window size
k=15 # k-mer size = 15
target_cov=25 # target read coverage = 25x
os_reads_file="os_reads.txt" # list of opposite sex read paths
ss_reads_file="ss_reads.txt" # list of same sex read paths

# Function to check output file exists before next step
check_output () {
  outfile="$1"
  if [[ -f "$outfile" ]]; then
    echo "Output file: $outfile"
    return 0
  else
    echo "Error: Output file ($outfile) not found."
    exit 1
  fi
}

# Help message
show_help () {
  echo "\
Runs k-mer-based SDR identification from Barrera-Redondo et al. (2025) and 
Carvalho and Clark (2013) using k=$k. Reads must be from opposite sex (os)
with respect to the assembly, with a set of same sex (ss) reads used to
validate assembly k-mers.

Usage: sbatch $job_name.sbatch <genome.fasta(.gz)> <asm_id>

Input:
  <asm_id>                assembly sample ID
  <genome.fasta(.gz)>     genome assembly FASTA

Output:
  ${job_name}_results/
  ├── <genome>_$win_str.contig_result  # log for YGS contig mode
  ├── <genome>_$win_str.ctg_rep.gz     # YGS bit-array of genomic repeat ${k}-mers found in > 1 contig
  ├── <genome>_$win_str.fasta.gz       # genome FASTA split into $win_str sliding windows
  ├── <genome>_$win_str.gen_rep.gz     # YGS bit-array of genomic repeat ${k}-mers
  ├── <genome>_${win_str}_${os_reads}_k${k}_ss_reads_k${k}.final_result
  ├── os_ids_table.tsv             # table of sample IDs and metadata (os)
  ├── os_reads.txt                 # list of read file paths (os)
  ├── os_reads_k$k.fasta.gz        # ${k}-mers extracted into FASTA (os)
  ├── os_reads_k$k.jf              # ${k}-mer Jellyfish binary hash (os)
  ├── os_reads_k$k.trace.gz        # YGS ${k}-mer bit-array (ss)
  ├── os_reads_k$k.trace_result    # log for YGS trace mode (os)
  ├── ss_ids_table.tsv             # table of sample IDs and metadata (ss)
  ├── ss_reads.txt                 # list of same sex read file paths (ss)
  ├── ss_reads_k$k.fasta.gz        # ${k}-mers extracted into FASTA (ss)
  ├── ss_reads_k$k.jf              # ${k}-mer Jellyfish binary hash (ss)
  ├── ss_reads_k$k.trace.gz        # YGS ${k}-mer bit-array (ss)
  └── ss_reads_k$k.trace_result    # log for YGS trace mode (ss)

Requires:
 - seqkit v2.3.1+ (https://bioinf.shenwei.me/seqkit/)
 - R v6+ (https://www.r-project.org)
 - Jellyfish v2.3.1+ (https://github.com/gmarcais/jellyfish)
 - pigz v2.8+ (https://zlib.net/pigz/)
 - Perl v5.23.1 (https://www.perl.org/)"
}
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]; then
  show_help
  exit 0
# elif [[ "$#" -ne 1 ]]; then
#   echo "Error: 1 positional argument required."
#   echo
#   show_help
#   exit 1
fi

# Catch status codes in pipes
set -o pipefail

# Print Slurm job information
[[ -n "$SLURM_JOB_ID" ]] && echo "\
==========================================
SLURM_JOB_ID = $SLURM_JOB_ID
SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST
TMPDIR = $TMPDIR
=========================================="

# Load conda
module purge
cond="$HOME/.conda_for_sbatch.sh"
if [[ -f "$cond" ]]; then
  source "$cond"
else
  echo "Error sourcing $cond"
  exit 1
fi

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Multithread support
nthrd="${SLURM_CPUS_PER_TASK:-1}"

# Download YGS.pl to scripts_dir (if needed)
[[ -f "$ygs_pl" ]] || wget -O "$ygs_pl" "$ygs_url"

# Input
genome="${1:-"assemblies/corteva_v2_ha.fa"}" # genome FASTA of one sex (e.g., female, U)
asm_id="${2:-"SL-CT1-FG-3"}" # genome assembly ID

# Output
mkdir -p "$outdir" # create main output directory
genome_basename="$(basename "$genome")" # genome filename basename
genome_base="${genome_basename%%.fa*}" # genome filename base

genomes=("$genome")
asm_ids=("$asm_id")

for i in "${!genomes[@]}"; do
  sbatch "$ygs_job" "${genomes[i]}" "${asm_ids[i]}" "$outdir"
  sbatch "$depth_job" "${genomes[i]}" "${asm_ids[i]}" "$outdir"
done

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
