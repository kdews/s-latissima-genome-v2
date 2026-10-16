#!/bin/bash
#SBATCH -J ygs_plot
#SBATCH -o %x_logs/%x_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=15g
#SBATCH --cpus-per-task=16

# Script variables
# Extract job name from script
if [[ -n ${SLURM_JOB_NAME:-} ]]; then
  job_name="$SLURM_JOB_NAME"
elif [[ -n ${BASH_SOURCE[0]:-} && ${BASH_SOURCE[0]} != "bash" ]]; then
  job_name="$(basename "${BASH_SOURCE[0]}")"
  job_name="${job_name%%.*}"
else
  job_name="ygs_plot"
fi
rscript="$job_name.R" # R script to run
scripts_dir="s-latissima-genome-v2" # scripts directory
[[ -d "$scripts_dir" ]] && rscript="$scripts_dir/$rscript"

# Help message
show_help () {
  echo "\
Plots result of k-mer-based SDR identification (Carvalho and Clark 2013,
Barrera-Redondo et al. 2025).

Usage: sbatch $job_name.sbatch <ygs_out> [outdir]

Input:
  ygs_out     path to YGS output directory for specific genome
  outdir     output directory for plot images (optional)

Output:
  <ygs_out>/
  ├── ${win_str}_windows.contig_result                  # log for YGS contig mode
  ├── ${win_str}_windows.ctg_rep.gz                     # YGS bit-array of genomic repeat ${k}-mers in > 1 contig
  ├── ${win_str}_windows.fasta.gz                       # genome FASTA split into $win_str sliding windows
  ├── ${win_str}_windows.gen_rep.gz                     # YGS bit-array of genomic repeat ${k}-mers
  ├── ${win_str}_windows_os_k${k}_ss_k${k}.final_result # speaks for itself
  ├── os_table.tsv                                      # table of sample IDs and metadata (os)
  ├── os.txt                                            # list of read file paths (os)
  ├── os_k$k.fasta.gz                                   # ${k}-mers extracted into FASTA (os)
  ├── os_k$k.jf                                         # ${k}-mer Jellyfish binary hash (os)
  ├── os_k$k.trace.gz                                   # YGS ${k}-mer bit-array (ss)
  ├── os_k$k.trace_result                               # log for YGS trace mode (os)
  ├── ss_table.tsv                                      # table of sample IDs and metadata (ss)
  ├── ss.txt                                            # list of same sex read file paths (ss)
  ├── ss_k$k.fasta.gz                                   # ${k}-mers extracted into FASTA (ss)
  ├── ss_k$k.jf                                         # ${k}-mer Jellyfish binary hash (ss)
  ├── ss_k$k.trace.gz                                   # YGS ${k}-mer bit-array (ss)
  └── ss_k$k.trace_result                               # log for YGS trace mode (ss)

Requires:
 - R v6+ (https://www.r-project.org)"
}
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]; then
  show_help
  exit 0
elif [[ "$#" -ne 1 ]]; then
  echo "Error: 1 positional argument required."
  echo
  show_help
  exit 1
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

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Input
ygs_out="$1" # YGS output directory
outdir="${2:-$scripts_dir}" # output directory for plot images (optional)

# Analyze outputs
# Load R module
module --latest load r
Rscript --version || exit 1
# Run Rscript
set -x
Rscript "$rscript" "$ygs_out" "$outdir"
{ set +x; } 2>/dev/null

