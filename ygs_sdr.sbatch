#!/bin/bash
#SBATCH -J ygs_sdr
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p htcondor
#SBATCH --time=2-0
#SBATCH --mem=59g
#SBATCH --cpus-per-task=16

# Script variables
# Extract job name from script
if [[ -n ${SLURM_JOB_NAME:-} ]]; then
  job_name="$SLURM_JOB_NAME"
elif [[ -n ${BASH_SOURCE[0]:-} && ${BASH_SOURCE[0]} != "bash" ]]; then
  job_name="$(basename "${BASH_SOURCE[0]}")"
  job_name="${job_name%%.*}"
else
  job_name="ygs_sdr"
fi
# URL to YGS.pl from Carvalho and Clark (2013)
ygs_url="https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/YGS.txt"
ygs_pl="YGS.pl" # local name of YGS script
scripts_dir="s-latissima-genome-v2" # scripts directory
main_outdir="sdr_id_results" # main output directory
ygs_outdir="$main_outdir/${job_name}_results" # job type output directory
window=500000 # window size = 500 kb
win_str="$((window / 1000))kb" # get string from window size
k=15 # k-mer size = 15
os_reads_file="os.txt" # list of paths to opposite sex reads
ss_reads_file="ss.txt" # list of paths to same sex reads

# Function to check output file exists before next step
check_output () {
  outfile="$1"
  if [[ -s "$outfile" ]]; then
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

Usage: sbatch $job_name.sbatch <genome>.fa <target> [scripts_dir]

Input:
  genome.fa             # genome assembly FASTA (.fa suffix REQUIRED)
  target                # sex chromosome to identify (default: U)
                          - U (female-specific)
                          - V (male-specific)
                          - UV run separately with each target
  scripts_dir           # optional: scripts directory (default: s-latissima-genome-v2)

Output:
  $ygs_outdir/<genome>/target_<target>/
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
 - seqkit v2.3.1+ (https://bioinf.shenwei.me/seqkit/)
 - Jellyfish v2.3.1+ (https://github.com/gmarcais/jellyfish)
 - pigz v2.8+ (https://zlib.net/pigz/)
 - Perl v5.23.1 (https://www.perl.org/)
 - R v6+ (https://www.r-project.org)"
}
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]; then
  show_help
  exit 0
elif [[ $# -lt 2 || $# -gt 3 ]]; then
  echo 'Expected genome, target and (optional) scripts_dir' >&2
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

# Input
genome="${1:-"assemblies/corteva_v2_ha.fa"}" # genome FASTA of one sex (e.g., female, U)
target="${2:-"U"}" # target sex
# $target specifies assembly sex as U or V
[[ "$target" == U || "$target" == V ]] || { echo 'target = U or V' >&2; exit 1; }
scripts_dir="${3:-$scripts_dir}" # optional user-set scripts directory name

# Download YGS.pl to scripts_dir (if needed), then resolve its absolute path
if [[ ! -s "$scripts_dir/$ygs_pl" ]]; then
  wget -O "$scripts_dir/$ygs_pl" "$ygs_url" || exit 1
fi
ygs_pl="$(realpath -e "$scripts_dir/$ygs_pl")" || exit 1

# Output
genome_base="$(basename "${genome%.fa}")" # genome filename base
ygs_out="$ygs_outdir/$genome_base/target_$target" # job output subdirectory
mkdir -p "$ygs_out" # create output directory
g_win_base="${win_str}_windows" # basename for windowized genome
os_base="os_k$k" # basename for opposite sex outputs
ss_base="ss_k$k" # basename for same sex outputs
out_final="${g_win_base}_${os_base}_${ss_base}.final_result" # final run output

# Verify input files: os/ss.txt lists of reads
if ! [[ -s "$ygs_out/$os_reads_file" && -s "$ygs_out/$ss_reads_file" ]]; then
  echo "Error: missing or empty YGS read lists in $ygs_out" >&2
  exit 1
fi

echo "Step 1: Slice genome into windows for YGS (Barrera-Redondo et al. 2025, p. 2138)"
# Load seqkit env
conda activate seqkit
seqkit --help > /dev/null 2>&1 || exit 1
echo "seqkit $(seqkit --help | grep "^Version" || exit 1)"
# Section genome FASTA into non-overlapping windows (default = 500kb)
# --greedy = don't toss end sequences shorter than window size
set -x
seqkit sliding -j "$nthrd" -t dna --greedy -s "$window" -W "$window" \
  -o "$ygs_out/$g_win_base.fasta.gz" "$genome" || exit 1
{ set +x; } 2>/dev/null
check_output "$ygs_out/$g_win_base.fasta.gz"
conda deactivate
echo; echo

echo "Step 2: YGS method from Carvalho and Clark (2013) supplement"
# https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/Supplemental_Material_Ykmer_ms_SM.pdf
# Load Jellyfish env
conda activate jellyfish
jellyfish --version || exit 1
echo "Preparing $k-mers for..."; echo
n_g="$((nthrd > 4 ? 4 : 1))" # concurrent generators (=4 if threads >4)
count_threads="$((nthrd > n_g ? nthrd - n_g : 1))" # calculate remaining threads
for reads_type in os ss; do
  # Set filter for minimum frequency = 5 (opposite sex) or 1 (same sex)
  if [[ "$reads_type" == os ]]; then
    type_str="opposite"
    min_freq=5
  elif [[ "$reads_type" == ss ]]; then
    type_str="same"
    min_freq=1
  fi
  echo "...$type_str sex reads (w.r.t. assembly sample)"
  base="${reads_type}_k$k"
  # Import read list into array
  mapfile -t reads < "$ygs_out/$reads_type.txt"
  # Write list of decompression commands to file
  printf 'pigz -dc %q\n' "${reads[@]}" > "$ygs_out/$reads_type.generators.txt"
  # Generate k-mer sequences from reads
  if [[ ! -s "$ygs_out/$base.jf" ]]; then # skip if output exists
    set -x
    jellyfish count -t "$count_threads" -m "$k" -s 10G -C --quality-start=33 \
      --min-quality=20 -g "$ygs_out/$reads_type.generators.txt" -G "$n_g" \
      -o "$ygs_out/$base.jf" || exit 1
    { set +x; } 2>/dev/null
  fi
  check_output "$ygs_out/$base.jf"
  # Convert k-mer FASTQ to FASTA
  if [[ ! -s "$ygs_out/$base.fasta.gz" ]]; then # skip if output exists
    set -x
    jellyfish dump --lower-count="$min_freq" "$ygs_out/$base.jf" | gzip -c \
      > "$ygs_out/$base.fasta.gz" || exit 1
    { set +x; } 2>/dev/null
  fi
  check_output "$ygs_out/$base.fasta.gz"
done
conda deactivate
echo; echo

echo "Step 3: Run YGS.pl"
# Activate YGS Perl env
# mamba create -n ygs conda-forge::perl conda-forge::perl-bit-vector \
#   conda-forge::perl-carp-clan bioconda::perl-perlio-gzip
conda activate ygs
perl --version | grep "version" || exit 1
perl -MBit::Vector -MCarp::Clan -MPerlIO::gzip \
  -e 'print "All required modules found\n"' || exit 1
# Change to output directory to control output location
cd "$ygs_out" || { echo "Erorr: Output directory ($ygs_out) not found."; exit 1; }
echo "YGS output will be in: $ygs_out"

# Produce bit-array from opposite sex k-mers, filtered by quality and frequency
echo "Producing bit-array from..."; echo
for reads_type in os ss; do
  if [[ "$reads_type" == os ]]; then
    type_str="opposite"
  elif [[ "$reads_type" == ss ]]; then
    type_str="same"
  fi
  echo "...${k}-mers of $type_str sex reads (w.r.t. assembly sample)"
  base="${reads_type}_k$k"
  if [[ ! -s "$base.trace.gz" ]]; then # skip if output exists
    set -x
    perl "$ygs_pl" kmer_size="$k" mode=trace trace="$base.fasta.gz" || exit 1
    { set +x; } 2>/dev/null
    check_output "$base.trace.gz"
  fi
  # Print log for trace run
  echo "Log file: $base.trace_result"
  grep -v '^$' "$base.trace_result" | sed 's/^/  /'
done

# Produce genomic repeat k-mers (copy number > 1) bit-array
echo "...genomic repeats"
if [[ ! -s "$g_win_base.gen_rep.gz" ]]; then # skip if output exists
  set -x
  perl "$ygs_pl" kmer_size="$k" mode=contig contig="$g_win_base.fasta.gz" \
    || exit 1
  { set +x; } 2>/dev/null
fi
check_output "$g_win_base.gen_rep.gz"
# Print log for genomic repeat run
echo "Log file: $g_win_base.contig_result"
grep -v '^$' "$g_win_base.contig_result" | grep -vP '^>|^defline' | sed 's/^/  /'
echo

# Predict SDR for each chromosome 
echo "Starting YGS final run..."
if [[ ! -s "$out_final" ]]; then # skip if output exists
  set -x
  perl "$ygs_pl" kmer_size="$k" mode=final_run contig="$g_win_base.fasta.gz" \
    trace="$os_base.trace.gz" gen_rep="$g_win_base.gen_rep.gz" \
    male_trace="$ss_base.trace.gz" || exit 1
  { set +x; } 2>/dev/null
fi
check_output "$out_final"
grep -v '^$' "$out_final" | grep -vP '^>|^GI' | sed 's/^/  /'

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
