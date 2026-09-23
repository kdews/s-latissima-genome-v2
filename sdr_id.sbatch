#!/bin/bash
#SBATCH -J sdr_id
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p largemem
#SBATCH --time=1-0
#SBATCH --mem=100g
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
scripts_dir="s-latissima-genome-v2" # scripts directory
if [[ -d "$scripts_dir" ]]; then
  sdr_parse="$scripts_dir/$sdr_parse"
  ygs_pl="$scripts_dir/$ygs_pl"
fi
ygs_pl="$(realpath -e "$ygs_pl")"
k=15 # k-mer size = 15
target_cov=60 # target coverage with reads = 60x
outdir="${job_name}_results" # main output directory
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
Carvalho and Clark (2013) using k=$k. Input genome and reads must be from
opposite sexes.

Usage: sbatch $job_name.sbatch <asm_id> <genome.fasta(.gz)>

Input:
  <asm_id>                assembly sample ID
  <genome.fasta(.gz)>     genome assembly FASTA

Requires:
 - seqkit v2.3.1+ (https://bioinf.shenwei.me/seqkit/)
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
asm_id="${1:-"SL-CT1-FG-3"}" # genome assembly ID
genome="${2:-"assemblies/corteva_v2.fa"}" # genome FASTA of one sex (e.g., female, U)

# Output
mkdir -p "$outdir" # create output directory
genome_basename="$(basename "$genome")" # genome filename basename
genome_base="${genome_basename%%.fa*}" # genome filename base
genome_windows="${genome_base}_500kb.fasta.gz" # genome split into windows
os="os_reads_k$k" # basename for opposite sex outputs
os_jf="$os.jf" # opposite sex k-mer counts
os_jfa="$os.fasta.gz" # opposite sex k-mer FASTA (compressed)
os_trc="$os.trace.gz" # opposite sex reads k-mer trace
os_trc_log="$os.trace_result" # log for opposite sex trace run
ss="ss_reads_k$k" # basename for same sex outputs
ss_jf="$ss.jf" # same sex k-mer counts
ss_jfa="$ss.fasta.gz" # same sex k-mer FASTA (compressed)
ss_trc="$ss.trace.gz" # same sex reads k-mer trace
ss_trc_log="$ss.trace_result" # log for same sex trace run
ctg="${genome_base}_500kb.gen_rep.gz" # genomic repeat k-mers
ctg_log="${genome_base}_500kb.contig_result" # log for contig run
out_final="${genome_base}_500kb_${os}_${ss}.final_result" # final run output

echo "Step 1: Slice genome into windows for YGS (Barrera-Redondo et al. 2025, p. 2138)"
# Load seqkit env
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
  -o "$outdir/$genome_windows" # output binned FASTA
  "$genome" # input genome FASTA
)
set -x
"${seq_cmd[@]}"
{ set +x; } 2>/dev/null
check_output "$outdir/$genome_windows"
echo; echo

echo "Step 2: Filter for reads of opposite sex and same location as assembly sample"
# Get assembly size
asm_size="$(seqkit stats -j "$nthrd" -T "$genome" | tail -n+2 | cut -f5)"
conda deactivate
# Load R module
module --latest load r
Rscript --version || exit 1
# Run Rscript
set -x
Rscript "$sdr_parse" "$target_cov" "$asm_size" "$asm_id" "$outdir"
{ set +x; } 2>/dev/null
check_output "$outdir/$os_reads_file"
check_output "$outdir/$ss_reads_file"
module purge
echo; echo

echo "Step 3: YGS method from Carvalho and Clark (2013) supplement"
# https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/Supplemental_Material_Ykmer_ms_SM.pdf
# Load Jellyfish env
conda activate jellyfish
jellyfish --version || exit 1
echo "Preparing $k-mers for..."
echo
# Generate k-mer sequences from opposite sex Illumina reads
echo "1. Opposite sex Illumina reads (w.r.t. assembly sample)"
if ! [[ -f "$outdir/$os_jf" ]]; then # skip if output exists
  # Import list into array
  # mapfile -t os_reads < "$outdir/$os_reads_file"
  mapfile -t os_reads < <(head -n2 "$outdir/$os_reads_file") # DELETE
  # Write list of decompression commands to file
  os_gen=os_generators.txt # generator filename
  printf 'pigz -dc %q\n' "${os_reads[@]}" > "$outdir/$os_gen"
  n_g=4 # number of concurrent generators
  count_cmd=(
    jellyfish count
    -t "$((nthrd - n_g))" # threads
    -m "$k" # k-mer size = 15 (script default)
    -s 10G # initial hash size
    -C # count both strand, canonical representation
    --quality-start=33 # ASCII for quality values (Illumina 1.9)
    --min-quality=20 # minimum quality (hard masked)
    -g "$outdir/$os_gen" # list of decompression commands
    -G "$n_g" # number of concurrent generators to run at once
    -o "$outdir/$os_jf" # output filename
  )
  set -x
  "${count_cmd[@]}"
  { set +x; } 2>/dev/null
fi
check_output "$outdir/$os_jf" && rm "$outdir/$os_gen"
# Convert k-mer FASTQ to FASTA, filter for minimum frequency = 5
if ! [[ -f "$outdir/$os_jfa" ]]; then # skip if output exists
  set -x
  jellyfish dump --lower-count=5 "$outdir/$os_jf" | gzip -c > "$outdir/$os_jfa"
  { set +x; } 2>/dev/null
fi
check_output "$outdir/$os_jfa"
echo

# Generate k-mer sequences from same sex Illumina reads
echo "2. Same sex Illumina reads (w.r.t. assembly sample)"
if ! [[ -f "$outdir/$ss_jf" ]]; then # skip if output exists
  # Import list into array
  # mapfile -t ss_reads < "$outdir/$ss_reads_file"
  mapfile -t ss_reads < <(head -n2 "$outdir/$ss_reads_file") # DELETE
  # Write list of decompression commands to file
  ss_gen=ss_generators.txt # generator filename
  printf 'pigz -dc %q\n' "${ss_reads[@]}" > "$outdir/$ss_gen"
  n_g=4 # number of concurrent generators
  count_cmd=(
    jellyfish count
    -t "$((nthrd - n_g))" # threads
    -m "$k" # k-mer size = 15 (script default)
    -s 10G # initial hash size
    -C # count both strand, canonical representation
    --quality-start=33 # ASCII for quality values (Illumina 1.9)
    --min-quality=20 # minimum quality (hard masked)
    -g "$outdir/$ss_gen" # list of decompression commands
    -G "$n_g" # number of concurrent generators to run at once
    -o "$outdir/$ss_jf" # output filename
  )
  set -x
  "${count_cmd[@]}"
  { set +x; } 2>/dev/null
fi
check_output "$outdir/$ss_jf" && rm "$outdir/$ss_gen"
# Convert k-mer FASTQ to FASTA, filter for minimum frequency = 1
if ! [[ -f "$outdir/$ss_jfa" ]]; then # skip if output exists
  set -x
  jellyfish dump --lower-count=1 "$outdir/$ss_jf" | gzip -c > "$outdir/$ss_jfa"
  { set +x; } 2>/dev/null
fi
check_output "$outdir/$ss_jfa"
conda deactivate
echo; echo

echo "Step 4: Run YGS.pl"
# Activate YGS env
# mamba create -n ygs conda-forge::perl conda-forge::perl-bit-vector \
#   conda-forge::perl-carp-clan bioconda::perl-perlio-gzip
conda activate ygs
perl --version | grep "version" || exit 1
perl -MBit::Vector -MCarp::Clan -MPerlIO::gzip \
  -e 'print "All required modules found\n"' || exit 1
# Change to output directory to control output location
cd "$outdir" || { echo "Erorr: Output directory ($outdir) not found."; exit 1; }
echo "YGS output will be in: $outdir"

# Produce k-mers bit-array from opposite sex reads, filtered by quality and frequency
echo "Producing bit-array from..."
echo
echo "1. Opposite sex Illumina reads (w.r.t. assembly sample)"
if ! [[ -f "$os_trc" ]]; then # skip if output exists
  set -x
  perl "$ygs_pl" kmer_size="$k" mode=trace trace="$os_jfa"
  { set +x; } 2>/dev/null
fi
check_output "$os_trc"
echo "Log file: $os_trc_log"
grep -v '^$' "$os_trc_log" | sed 's/^/  /'
echo

# Produce validating bit-array from same sex reads, filtered by quality and frequency
echo "2. Same sex Illumina reads (w.r.t. assembly sample)"
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
