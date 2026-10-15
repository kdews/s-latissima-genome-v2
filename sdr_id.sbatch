#!/bin/bash
#SBATCH -J sdr_id
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p htcondor
#SBATCH --time=01:00:00
#SBATCH --mem=40g
#SBATCH --cpus-per-task=16

# Script variables
# Extract job name from script
if [[ -n ${SLURM_JOB_NAME:-} ]]; then
  job_name="$SLURM_JOB_NAME"
elif [[ -n ${BASH_SOURCE[0]:-} && ${BASH_SOURCE[0]} != "bash" ]]; then
  job_name="$(basename "${BASH_SOURCE[0]}")"
  job_name="${job_name%%.*}"
else
  job_name="sdr_id"
fi
main_outdir="${job_name}_results" # main SDR output directory
# YGS SDR job variables
ygs_job_name="ygs_sdr"
ygs_job="$ygs_job_name.sbatch" # Slurm job submission script
ygs_outdir="${ygs_job_name}_results" # output subdirectory
# URL to YGS.pl from Carvalho and Clark (2013)
ygs_url="https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/YGS.txt"
ygs_pl="YGS.pl" # local name of YGS script
sdr_parse="sdr_parse_multiqc.R" # R script parses reads for YGS from MultiQC data
# Read depth SDR job variables
depth_job_name="depth_sdr"
depth_job="$depth_job_name.sbatch" # Slurm job submission script
depth_outdir="${depth_job_name}_results" # output subdirectory
scripts_dir="s-latissima-genome-v2" # scripts directory
if [[ -d "$scripts_dir" ]]; then
  ygs_job="$scripts_dir/$ygs_job"
  ygs_pl="$scripts_dir/$ygs_pl"
  sdr_parse="$scripts_dir/$sdr_parse"
  depth_job="$scripts_dir/$depth_job"
fi
ygs_pl="$(realpath -e "$ygs_pl")"
window=500000 # window size = 500 kb
win_str="$((window / 1000))kb" # get string from window size
k=15 # k-mer size = 15
target_cov=25 # target read coverage = 25x
os_reads_file="os.txt" # list of paths to opposite sex reads
ss_reads_file="ss.txt" # list of paths to same sex reads
depth_reads_file="os_ss.txt" # paths to same and opposite sex reads
depth_tab_file="os_ss_table.tsv" # table of sample IDs

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
  asm_id                assembly sample ID
  genome.fasta(.gz)     genome assembly FASTA

Output:
  $main_outdir/
  ├── $ygs_outdir/
  │   └── <genome>/
  │       ├── ${win_str}_windows.contig_result                  # log for YGS contig mode
  │       ├── ${win_str}_windows.ctg_rep.gz                     # YGS bit-array of genomic repeat ${k}-mers in > 1 contig
  │       ├── ${win_str}_windows.fasta.gz                       # genome FASTA split into $win_str sliding windows
  │       ├── ${win_str}_windows.gen_rep.gz                     # YGS bit-array of genomic repeat ${k}-mers
  │       ├── ${win_str}_windows_os_k${k}_ss_k${k}.final_result # speaks for itself
  │       ├── os_table.tsv                                      # table of sample IDs and metadata (os)
  │       ├── os.txt                                            # list of read file paths (os)
  │       ├── os_k$k.fasta.gz                                   # ${k}-mers extracted into FASTA (os)
  │       ├── os_k$k.jf                                         # ${k}-mer Jellyfish binary hash (os)
  │       ├── os_k$k.trace.gz                                   # YGS ${k}-mer bit-array (ss)
  │       ├── os_k$k.trace_result                               # log for YGS trace mode (os)
  │       ├── ss_table.tsv                                      # table of sample IDs and metadata (ss)
  │       ├── ss.txt                                            # list of same sex read file paths (ss)
  │       ├── ss_k$k.fasta.gz                                   # ${k}-mers extracted into FASTA (ss)
  │       ├── ss_k$k.jf                                         # ${k}-mer Jellyfish binary hash (ss)
  │       ├── ss_k$k.trace.gz                                   # YGS ${k}-mer bit-array (ss)
  │       └── ss_k$k.trace_result                               # log for YGS trace mode (ss)
  └── $depth_outdir/
      ├── 
      ├── 
      ├── 
      ├── 
      └── 

Requires:
 - seqkit v2.3.1+ (https://bioinf.shenwei.me/seqkit/)
 - R v6+ (https://www.r-project.org)
 - Jellyfish v2.3.1+ (https://github.com/gmarcais/jellyfish)
 - pigz v2.8+ (https://zlib.net/pigz/)
 - Perl v5.23.1 (https://www.perl.org/)
 - bwa-mem2 v2.2.1 (https://bio-bwa.sourceforge.net/)
 - samtools v1.24 (https://www.htslib.org/)
 - mosdepth v0.3.14 (https://github.com/brentp/mosdepth)"
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
genomes=("$genome")
asm_ids=("$asm_id")

# Iterate through genomes
for i in "${!genomes[@]}"; do
  # Output
  genome_basename="$(basename "$genome")" # genome filename basename
  genome_base="${genome_basename%%.fa*}" # genome filename base
  genome_dict="${genome%%.fa*}.dict" # genome GATK4 dictionary (same path as FASTA)
  ygs_out="$main_outdir/$ygs_outdir/$genome_base" # YGS job output subdirectory
  depth_out="$main_outdir/$depth_outdir/$genome_base" # depth job output subdirectory
  mkdir -p "$ygs_out" "$depth_out" # create output directories

  echo "Filtering for reads of opposite sex and same location as assembly sample..."
  # Get assembly size
  # Load seqkit env
  conda activate seqkit
  seqkit --help > /dev/null 2>&1 || exit 1
  echo "seqkit $(seqkit --help | grep "^Version" || exit 1)"
  asm_size="$(seqkit stats -j "$nthrd" -T "$genome" | tail -n+2 | cut -f5)"
  conda deactivate
  # Load R module
  module --latest load r
  Rscript --version || exit 1
  # Run Rscript
  set -x
  Rscript "$sdr_parse" "$target_cov" "$asm_size" "$asm_id" "$ygs_out" "$depth_out"
  { set +x; } 2>/dev/null
  check_output "$ygs_out/$os_reads_file"
  check_output "$ygs_out/$ss_reads_file"
  check_output "$depth_out/$depth_reads_file"
  check_output "$depth_out/$depth_tab_file"
  module purge
  echo; echo
  
  # Prep genome for depth analysis
  # Index genome for bwa-mem2 alignment
  conda activate bwa-mem2
  set -x
  bwa-mem2 index -p "$depth_out/$genome_base" "$genome"
  { set +x; } 2>/dev/null
  conda deactivate
  # Load GATK4
  module purge
  module load gcc/13.3.0 gatk/4.5.0.0 python/3.13.2 openjdk/21.0.0_35
  gatk --version || exit 1
  # Build GATK4 sequence dictionary from genome FASTA
  set -x
  gatk CreateSequenceDictionary -R "$genome" -O "$genome_dict"
  { set +x; } 2>/dev/null
  module purge

  # Submit YGS and depth jobs
  sbatch "$ygs_job" "${genomes[i]}" "$main_outdir"
  sbatch "$depth_job" "${genomes[i]}" "$main_outdir"
done

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
