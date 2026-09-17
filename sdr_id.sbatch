#!/bin/bash
#SBATCH -J sdr_id
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p oneweek
#SBATCH --time=1-0
#SBATCH --mem=180g
#SBATCH --cpus-per-task=12

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
rscript="$job_name.R" # Rscript
ygs_url="https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/YGS.txt"
ygs_pl="YGS.pl"
scripts_dir="s-latissima-genome-v2" # scripts directory
if [[ -d "$scripts_dir" ]]; then
  rscript="$scripts_dir/$rscript"
  ygs_pl="$scripts_dir/$ygs_pl"
fi
ygs_pl="$(realpath -e "$ygs_pl")"
k=15 # k-mer size = 15
outdir="${job_name}_results" # main output directory

# Function to check output file exists before next step
check_output () {
  outfile="$1"
  if [[ -f "$outfile" ]]
  then
    echo "Output file: $outfile"
  else
    echo "Error: Output file ($outfile) not found."
    exit 1
  fi
}

# Help message
show_help () {
  echo "\
Runs k-mer-based SDR identification from Barrera-Redondo et al. (2025) using
k=$k. Input genome and reads must be from different sexes.

Usage: sbatch $job_name.sbatch <genome.fasta(.gz)> <reads.fastq(.gz)>

Requires:
 - seqkit v2.3.1+ (https://bioinf.shenwei.me/seqkit/)
 - Jellyfish v2.3.1+ (https://github.com/gmarcais/jellyfish)
 - Perl v5.23.1 (https://www.perl.org/)"
}
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  show_help
  exit 0
# elif [[ "$#" -ne 1 ]]
# then
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
if [[ -f "$cond" ]]
then
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
genome="${1:-assemblies/corteva_v2.fa}" # full female (U) genome FASTA
in_fastq="${2:-wgs/example.fastq.gz}" # input male (V) reads

# Output
mkdir -p "$outdir" # create output directory
genome_basename="$(basename "$genome")" # genome filename basename
genome_base="${genome_basename%%.fa*}" # genome filename base
genome_windows="$genome_base.500kbwin.fasta.gz" # genome split into windows
outname="output" # reads filename base
out_jf="$outname.jf" # k-mer counts
out_jfa="$outname.fasta.gz" # k-mer FASTA (compressed)
out_m_trc="$outname.trace.gz" # male k-mer trace
out_gen_rep="$outname.gen_rep.gz" # genomic repeat k-mers
out_final="$genome_base.500kbwin_${outname}.final_run" # final run output

# Convert genome into windows for YGS with Barrera-Redondo et al. (2025, p. 2138)

# Load seqkit env
module purge
conda activate seqkit
seqkit --help > /dev/null 2>&1 || exit 1
echo "seqkit $(seqkit --help | grep "^Version" || exit 1)"

# Section genome FASTA into 500-kb windows
seq_cmd=(
  seqkit sliding
  -j "$nthrd" # threads
  -t dna # sequence type
  --greedy # don't toss end sequences shorter than window size
  -s 500000 # 500-kb step size
  -W 500000 # 500-kb window size
  -o "$outdir/$genome_windows"
  "$genome"
)
echo "${seq_cmd[*]}"
"${seq_cmd[@]}"
check_output "$outdir/$genome_windows"

conda deactivate
echo

# YGS method from Carvalho and Clark (2013) supplement
# https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/Supplemental_Material_Ykmer_ms_SM.pdf

# Load Jellyfish env
conda activate jellyfish
jellyfish --version || exit 1

# Generate k-mer sequences from Illumina reads
count_cmd=(
  jellyfish count
  -t "$nthrd" # threads
  -m "$k" # k-mer size = 15 (script default)
  -s 10G # initial hash size
  -C # count both strand, canonical representation
  --quality-start=33 # ASCII for quality values (CHECK??)
  --min-quality=20 # minimum quality (hard masked)
  -o "$outdir/$out_jf" # output filename
)
echo "${count_cmd[*]} <($in_fastq)"
"${count_cmd[@]}" <("$in_fastq")
check_output "$outdir/$out_jf"

# Convert k-mer FASTQ to FASTA, filter for minimum frequency = 5
set -x
jellyfish dump --lower-count=5 | gzip -c > "$outdir/$out_jfa"
set +x
check_output "$outdir/$out_jfa"

conda deactivate
echo

# # Run YGS.pl

# # Activate YGS env
# # mamba create -n ygs conda-forge::perl conda-forge::perl-bit-vector \
# #   conda-forge::perl-carp-clan bioconda::perl-perlio-gzip
# conda activate ygs
# perl --version | grep "version" || exit 1
# perl -MBit::Vector -MCarp::Clan -MPerlIO::gzip \
#   -e 'print "All required modules found\n"' || exit 1

# # Change to output directory to control output location
# cd "$outdir" || { echo "Erorr: Output directory ($outdir) not found."; exit 1; }

# # Produce male k‐mers bit‐array, filtered by quality and frequency
# set -x
# perl "$ygs_pl" kmer_size="$k" mode=trace trace="$out_jfa"
# set +x
# check_output "$out_m_trc"

# # Produce genomic repeat k‐mers (copy number > 1) bit‐array
# set -x
# perl "$ygs_pl" kmer_size="$k" mode=contig contig="$genome_windows"
# set +x
# check_output "$out_gen_rep"

# # Predict SDR for each chromosome 
# set -x
# perl "$ygs_pl" kmer_size="$k" mode=final_run contig="$genome_windows" \
#   trace="$out_m_trc" gen_rep="$out_gen_rep"
# set +x
# check_output "$out_final"
