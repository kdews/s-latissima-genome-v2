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
sdr_parse="sdr_parse_multiqc.R" # R script parses reads for YGS from MultiQC data
scripts_dir="s-latissima-genome-v2" # scripts directory
if [[ -d "$scripts_dir" ]]; then
  sdr_parse="$scripts_dir/$sdr_parse"
  ygs_pl="$scripts_dir/$ygs_pl"
fi
ygs_pl="$(realpath -e "$ygs_pl")"
main_outdir="${2:-sdr_id_results}" # main output directory
ygs_outdir="$main_outdir/${job_name}_results" # job type output directory
window=500000 # window size = 500 kb
win_str="$((window / 1000))kb" # get string from window size
k=15 # k-mer size = 15
target_cov=25 # target read coverage = 25x
os_reads_file="os.txt" # list of paths to opposite sex reads
ss_reads_file="ss.txt" # list of paths to same sex reads

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

Usage: sbatch $job_name.sbatch <genome.fasta(.gz)> [main_outdir]

Input:
  genome.fasta(.gz)     # genome assembly FASTA
  main_outdir           # optional: main output directory (default: sdr_id_results)

Output:
  $ygs_outdir/<genome>/
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
elif [[ "$#" -lt 1 ]]; then
  echo "Error: at least 1 positional arguments required."
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

# Download YGS.pl to scripts_dir (if needed)
[[ -f "$ygs_pl" ]] || wget -O "$ygs_pl" "$ygs_url"

# Input
genome="${1:-"assemblies/corteva_v2_ha.fa"}" # genome FASTA of one sex (e.g., female, U)

# Output
genome_basename="$(basename "$genome")" # genome filename basename
genome_base="${genome_basename%%.fa*}" # genome filename base
ygs_out="$ygs_outdir/$genome_base" # job output subdirectory
mkdir -p "$ygs_out" # create output directories
genome_windows="${win_str}_windows.fasta.gz" # genome split into windows
os="os_k$k" # basename for opposite sex outputs
os_jf="$os.jf" # opposite sex k-mer counts
os_jfa="$os.fasta.gz" # opposite sex k-mer FASTA (compressed)
os_trc="$os.trace.gz" # opposite sex reads k-mer trace
os_trc_log="$os.trace_result" # log for opposite sex trace run
ss="ss_k$k" # basename for same sex outputs
ss_jf="$ss.jf" # same sex k-mer counts
ss_jfa="$ss.fasta.gz" # same sex k-mer FASTA (compressed)
ss_trc="$ss.trace.gz" # same sex reads k-mer trace
ss_trc_log="$ss.trace_result" # log for same sex trace run
ctg="${win_str}_windows.gen_rep.gz" # genomic repeat k-mers
ctg_log="${win_str}_windows.contig_result" # log for contig run
out_final="${win_str}_windows_${os}_${ss}.final_result" # final run output

echo "Step 1: Slice genome into windows for YGS (Barrera-Redondo et al. 2025, p. 2138)"
# Load seqkit env
conda activate seqkit
seqkit --help > /dev/null 2>&1 || exit 1
echo "seqkit $(seqkit --help | grep "^Version" || exit 1)"
# Section genome FASTA into non-overlapping windows (default = 500kb)
seq_cmd=(
  seqkit sliding
  -j "$nthrd" # threads
  -t dna # sequence type
  --greedy # don't toss end sequences shorter than window size
  -s "$window" # step size
  -W "$window" # window size
  -o "$ygs_out/$genome_windows" # output binned FASTA
  "$genome" # input genome FASTA
)
set -x
"${seq_cmd[@]}"
{ set +x; } 2>/dev/null
check_output "$ygs_out/$genome_windows"
echo; echo
conda deactivate

echo "Step 2: YGS method from Carvalho and Clark (2013) supplement"
# https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/Supplemental_Material_Ykmer_ms_SM.pdf
# Load Jellyfish env
conda activate jellyfish
jellyfish --version || exit 1
echo "Preparing $k-mers for..."
echo
# Generate k-mer sequences from opposite sex Illumina reads
echo "1. Opposite sex Illumina reads (w.r.t. assembly sample)"
if ! [[ -f "$ygs_out/$os_jf" ]]; then # skip if output exists
  # Import list into array
  mapfile -t os_reads < "$ygs_out/$os_reads_file"
  # Write list of decompression commands to file
  os_gen=os_generators.txt # generator filename
  printf 'pigz -dc %q\n' "${os_reads[@]}" > "$ygs_out/$os_gen"
  n_g=4 # number of concurrent generators
  count_cmd=(
    jellyfish count
    -t "$((nthrd - n_g))" # threads
    -m "$k" # k-mer size = 15 (script default)
    -s 10G # initial hash size
    -C # count both strand, canonical representation
    --quality-start=33 # ASCII for quality values (Illumina 1.9)
    --min-quality=20 # minimum quality (hard masked)
    -g "$ygs_out/$os_gen" # list of decompression commands
    -G "$n_g" # number of concurrent generators to run at once
    -o "$ygs_out/$os_jf" # output filename
  )
  set -x
  "${count_cmd[@]}"
  { set +x; } 2>/dev/null
fi
check_output "$ygs_out/$os_jf" && rm "$ygs_out/$os_gen"
# Convert k-mer FASTQ to FASTA, filter for minimum frequency = 5
if ! [[ -f "$ygs_out/$os_jfa" ]]; then # skip if output exists
  set -x
  jellyfish dump --lower-count=5 "$ygs_out/$os_jf" | gzip -c > "$ygs_out/$os_jfa"
  { set +x; } 2>/dev/null
fi
check_output "$ygs_out/$os_jfa"
echo

# Generate k-mer sequences from same sex Illumina reads
echo "2. Same sex Illumina reads (w.r.t. assembly sample)"
if ! [[ -f "$ygs_out/$ss_jf" ]]; then # skip if output exists
  # Import list into array
  mapfile -t ss_reads < "$ygs_out/$ss_reads_file"
  # Write list of decompression commands to file
  ss_gen=ss_generators.txt # generator filename
  printf 'pigz -dc %q\n' "${ss_reads[@]}" > "$ygs_out/$ss_gen"
  n_g=4 # number of concurrent generators
  count_cmd=(
    jellyfish count
    -t "$((nthrd - n_g))" # threads
    -m "$k" # k-mer size = 15 (script default)
    -s 10G # initial hash size
    -C # count both strand, canonical representation
    --quality-start=33 # ASCII for quality values (Illumina 1.9)
    --min-quality=20 # minimum quality (hard masked)
    -g "$ygs_out/$ss_gen" # list of decompression commands
    -G "$n_g" # number of concurrent generators to run at once
    -o "$ygs_out/$ss_jf" # output filename
  )
  set -x
  "${count_cmd[@]}"
  { set +x; } 2>/dev/null
fi
check_output "$ygs_out/$ss_jf" && rm "$ygs_out/$ss_gen"
# Convert k-mer FASTQ to FASTA, filter for minimum frequency = 1
if ! [[ -f "$ygs_out/$ss_jfa" ]]; then # skip if output exists
  set -x
  jellyfish dump --lower-count=1 "$ygs_out/$ss_jf" | gzip -c > "$ygs_out/$ss_jfa"
  { set +x; } 2>/dev/null
fi
check_output "$ygs_out/$ss_jfa"
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
echo "Producing bit-array from..."
echo
echo "1. ${k}-mers of opposite sex reads (w.r.t. assembly sample)"
if ! [[ -f "$os_trc" ]]; then # skip if output exists
  set -x
  perl "$ygs_pl" kmer_size="$k" mode=trace trace="$os_jfa"
  { set +x; } 2>/dev/null
fi
check_output "$os_trc"
echo "Log file: $os_trc_log"
grep -v '^$' "$os_trc_log" | sed 's/^/  /'
echo

# Produce validating bit-array from same sex k-mers, filtered by quality and frequency
echo "2. ${k}-mers of same sex reads (w.r.t. assembly sample)"
if ! [[ -f "$ss_trc" ]]; then # skip if output exists
  set -x
  perl "$ygs_pl" kmer_size="$k" mode=trace trace="$ss_jfa"
  { set +x; } 2>/dev/null
fi
check_output "$ss_trc"
echo "Log file: $ss_trc_log"
grep -v '^$' "$ss_trc_log" | sed 's/^/  /'
echo

# Produce genomic repeat k-mers (copy number > 1) bit-array
echo "3. Genomic repeats"
if ! [[ -f "$ctg" ]]; then # skip if output exists
  set -x
  perl "$ygs_pl" kmer_size="$k" mode=contig contig="$genome_windows"
  { set +x; } 2>/dev/null
fi
check_output "$ctg"
echo "Log file: $ctg_log"
grep -v '^$' "$ctg_log" | grep -vP '^>|^defline' | sed 's/^/  /'
echo

# Predict SDR for each chromosome 
echo "Starting YGS final run..."
if ! [[ -f "$out_final" ]]; then # skip if output exists
  set -x
  perl "$ygs_pl" kmer_size="$k" mode=final_run contig="$genome_windows" \
    trace="$os_trc" gen_rep="$ctg" male_trace="$ss_trc"
  { set +x; } 2>/dev/null
fi
check_output "$out_final"
grep -v '^$' "$out_final" | grep -vP '^>|^GI' | sed 's/^/  /'

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
