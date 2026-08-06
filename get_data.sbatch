#!/bin/bash
#SBATCH -J get_data
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p largemem
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=280g

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
# Max Plank's Edmond Dataverse API
DATAVERSE_API_BASE='https://edmond.mpg.de'

# Help message
show_help () {
  echo "\
Usign a TSV manifest of genome metadata, copies from local or downloads from
remote all genome FASTAs and reads. Uses curl or NCBI CLI to download genome
FASTAs and SRA Toolkit to download WGS reads.

Usage: sbatch $job_name.sbatch <in_tsv>

Required input:
  in_tsv      TSV with 10 named columns (see below), with commented (#) header

Input TSV columns:
  <Column_Name>           <Example>
  1. Label                corteva_v2
  2. Pretty_Name          Saccharina latissima (US v2)
  3. Genome_Source        local
  4. Genome_Method        file
  5. Genome_Locator       /project2/noujdine_61/kdeweese/latissima/all_s_latissima_assemblies/corteva_v2/Slattisima_v0.2_Genome_Polished_111925.fasta
  6. Genome_Selector      [blank]*
  7. Reads_Source         local
  8. Reads_Method         glob
  9. Reads_Type           pacbio
  10. Reads_Locator        /project2/noujdine_61/kdeweese/latissima/all_s_latissima_assemblies/corteva_v2/Slattisima_v0.2_ONT_raw_reads_Part[12].fastq

* Blank genome or read fields skip that data type for the row.

Supported genome sources:
  file        local FASTA (plain or gzip-compressed)
  url         direct FASTA URL (plain or gzip-compressed)
  url_tar     tar archive URL; Genome_Selector is the member basename
  figshare    Figshare article ID; selector is a FASTA filename or regex
  dataverse   DOI; selector is a FASTA filename or regex
  ncbi        NCBI assembly accession

Supported read sources:
  glob        local shell glob; matching FASTQ/FASTQ.GZ files are concatenated
  sra         one or more comma-separated SRA run accessions

Output:
  .
  ├── assemblies/             # directory containing all assemblies
  │   ├── <label>.fa          # renamed assembly FASTA (uncompressed)
  │   └── <label>.fa.fai      # FASTA index
  └── reads/                  # directory containing all reads
  │   └── <label>.fastq.gz    # renamed reads FASTQ (gzip-compressed)
  └── <out_tsv>               # output TSV of local paths

Requires:
 - curl v8+ (https://curl.se)
 - jq v1.7+ (https://jqlang.org)
 - pigz v2.8+ (https://zlib.net/pigz)
 - NCBI Datasets command-line tools (CLI) v18.31.0+ (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/)
 - UnZip v6.00+ (https://infozip.sourceforge.net/UnZip.html)
 - samtools v1.23+ (https://www.htslib.org)
 - SRA Toolkit v3.4.1+ (https://github.com/ncbi/sra-tools)"
}
# Catch help flag or missing argument(s)
if [[ $1 == "-h" || $1 == "--help" ]]
then
  show_help
  exit 0
elif [[ "$#" -ne 1 ]]
then
  echo "Error: 1 positional argument required."
  echo
  show_help
  exit 1
fi

set -o pipefail
shopt -s nullglob

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

# Multithread support
nthrd="${SLURM_CPUS_PER_TASK:-1}"

# Input
in_tsv="${1:-genomes_metadata.tsv}"
[[ -r $in_tsv ]] || die "Cannot read input TSV: $in_tsv"

# Output
work_dir="${job_name}_tmp" # temporary working directory
mkdir -p "$work_dir" assemblies reads # create output directories
parsed="$work_dir/parsed.records"

# Test required packages
echo "Checking required packages and versions..."
# Load samtools env
conda activate samtools
samtools --version | head -n1 || exit 1
conda deactivate
# Load curl and jq modules
module load usc curl jq
jq --version || exit 1
curl --version | head -n1 || exit 1
module purge
# Load unzip module
module load unzip
unzip -v | head -n1 || exit 1
module purge
# Load pigz env
conda activate pigz
pigz --version || exit 1
conda deactivate
# Load NCBI CLI env
conda activate ncbi
datasets --version | head -n1 || exit 1
# Load SRA tools
conda activate sratools
prefetch --version | head -n1 || exit 1
fasterq-dump --version | head -n1 || exit 1
conda deactivate

# Kill function
die () {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

# Log function
log () {
  printf '[%(%F %T)T] %s\n' -1 "$*"
}

# Copy a plain/gzipped genome FASTA to a normalized, uncompressed output path
norm_genome () {
  local input="$1" output="$2"
  log "Input: $input"
  log "Output: $output"
  # Load pigz env
  conda activate pigz
  if pigz -t -p "$nthrd" "$input" 2>/dev/null
  then
    pigz -cd -p "$nthrd" -- "$input" > "$output"
  else
    cp -- "$input" "$output"
  fi
  conda deactivate
  # Verify resulting file is not empty
  if [[ -s "$output" ]]
  then
    log "Successfully retrived genome: $output"
  else
    die "Created empty output: $output"
  fi
}

# Print the selected FASTA file record from a Figshare or Dataverse JSON response
select_fasta () {
  local selector="$1"
  jq -er --arg selector "$selector" '
    map(select(.name | test("\\.(fa|fna|fasta|fas)(\\.(gz|bgz))?$"; "i")))
    | (if $selector == "" then .
       else map(select((.name == $selector) or (.name | test($selector))))
       end)
    | if length == 1 then .[0]
      elif length == 0 then error("no matching FASTA")
      else error("multiple matching FASTAs; set Genome_Selector")
      end
  '
}

get_genome () {
  local label="$1" source="$2" method="$3" locator="$4" selector="$5"
  local output="assemblies/$label.fa"
  local download="$work_dir/$label.genome"
  local json record url member

  [[ "$source" == local || "$source" == remote ]] \
    || die "$label: Genome_Source must be either 'local' or 'remote'"
  log "$label: preparing genome ($source/$method)"

  case "$source:$method" in
    local:file)
      [[ -r "$locator" ]] || die "$label: cannot read genome: $locator"
      norm_genome "$locator" "$output"
      ;;
    remote:url)
      # Load curl module
      module load usc curl
      # Download from URL with curl
      curl -fL --retry 3 --output "$download" "$locator"
      module purge
      norm_genome "$download" "$output"
      ;;
    remote:url_tar)
      [[ -n "$selector" ]] || die "$label: url_tar requires Genome_Selector"
      # Load curl module
      module load usc curl
      # Download tarball from URL with curl
      curl -fL --retry 3 --output "$download.tar" "$locator"
      module purge
      tar -tf "$download.tar" \
        | awk -v name="$selector" 'BEGIN { FS="/" } $NF == name' \
        > "$work_dir/members"
      mapfile -t members < "$work_dir/members"
      (( ${#members[@]} == 1 )) \
        || die "$label: expected one FASTA named $selector; found ${#members[@]}"
      member="${members[0]}"
      tar -xOf "$download.tar" "$member" > "$download"
      norm_genome "$download" "$output"
      ;;
    remote:figshare)
      # Load curl and jq modules
      module load usc curl jq
      # Parse download URL from Figshare API
      json="$(curl -fsSL --retry 3 "https://api.figshare.com/v2/articles/$locator/files")"
      record="$(jq -ce '[.[] | {name: .name, url: .download_url}]' <<<"$json" \
        | select_fasta "$selector")" \
        || die "$label: could not select single Figshare FASTA"
      url="$(jq -r .url <<<"$record")"
      # Download from URL with curl
      curl -fL --retry 3 --output "$download" "$url"
      module purge
      norm_genome "$download" "$output"
      ;;
    remote:dataverse)
      # Load curl and jq modules
      module load usc curl jq
      # Parse download URL from Dataverse API
      json="$(curl -fsSL --retry 3 \
        "$DATAVERSE_API_BASE/api/datasets/:persistentId/?persistentId=doi:$locator")"
      record="$(jq -ce --arg base "$DATAVERSE_API_BASE" '
        [.data.latestVersion.files[].dataFile
          | {name: .filename, url: ($base + "/api/access/datafile/" + (.id | tostring))}]
      ' <<<"$json" | select_fasta "$selector")" \
        || die "$label: could not select one Dataverse FASTA"
      url="$(jq -r .url <<<"$record")"
      # Download from URL with curl
      curl -fL --retry 3 --output "$download" "$url"
      module purge
      norm_genome "$download" "$output"
      ;;
    remote:ncbi)
      # Load NCBI CLI env
      conda activate ncbi
      # Download genome accession with NCBI CLI
      datasets download genome accession "$locator" --no-progressbar \
        --include genome --filename "$download.zip"
      conda deactivate
      # Load unzip module
      module load unzip
      # Unzip NCBI dataset
      unzip -Z1 "$download.zip" \
        | awk -v acc="$locator" \
          '$0 ~ "^ncbi_dataset/data/" acc "/.*\\.(fna|fa|fasta)(\\.gz)?$"' \
        > "$work_dir/members"
      mapfile -t members < "$work_dir/members"
      (( ${#members[@]} == 1 )) \
        || die "$label: expected one NCBI FASTA for $locator, found ${#members[@]}"
      unzip -p "$download.zip" "${members[0]}" > "$download"
      module purge
      norm_genome "$download" "$output"
      ;;
    *)
      die "$label: unsupported genome source/method: $source/$method"
      ;;
  esac

  # Load samtools env
  conda activate samtools
  # Index genome FASTA
  samtools faidx "$output"
  conda deactivate
}

# Compress FASTQ input (if needed) and append to output file
append_fastq () {
  local input="$1" output="$2"
  log "Input: $input"
  log "Output: $output"
  log "Checking input FASTQ compression..."
  # Load pigz env
  conda activate pigz
  if pigz -t -p "$nthrd" "$input" 2>/dev/null
  then
    log "Concatenating compressed reads..."
    cat -- "$input" >> "$output"
  else
    log "Compressing then concatenating reads..."
    pigz -c -p "$nthrd" -- "$input" >> "$output"
  fi
  conda deactivate
}

get_reads () {
  local clean_label="$1" source="$2" method="$3" type="$4" locator="$5"
  local longread_output="reads/${clean_label}_$type.fastq.gz"
  local pe1_output="reads/${clean_label}_1.fastq.gz"
  local pe2_output="reads/${clean_label}_2.fastq.gz"
  local accession fastq
  local -a inputs accessions outputs

  [[ "$source" == local || "$source" == remote ]] \
    || die "$clean_label: Reads_Source must be either 'local' or 'remote'"
  log "$clean_label: preparing reads ($source/$method)"

  case "$source:$method" in
    local:glob)
      # Look for glob
      compgen -G "$locator" | sort > "$work_dir/read_inputs" || true
      mapfile -t inputs < "$work_dir/read_inputs"
      (( ${#inputs[@]} > 0 )) || die "$clean_label: glob matched no files: $locator"
      # Initialize empty output FASTQ file
      : > "$longread_output"
      # Append all FASTQs into single file
      for fastq in "${inputs[@]}"
      do
        append_fastq "$fastq" "$longread_output"
      done
      ;;
    remote:sra)
      IFS=',' read -r -a accessions <<<"$locator"
      for accession in "${accessions[@]}"
      do
        accession="${accession//[[:space:]]/}" # remove any whitespaces
        [[ "$accession" =~ ^[SED]RR[0-9]+$ ]] \
          || die "$clean_label: invalid SRA run accession: $accession" # check SRA
      # Load sra-tools env
      conda activate sratools
        # Fetch SRA
        prefetch --progress --max-size 60000000 \
          --output-directory "$TMPDIR" "$accession"
        # Make temporary directory for fasterq-dump
        dump_tmp="$TMPDIR/dump_tmp_$accession"
        mkdir -p "$dump_tmp"
        # Dump FASTQ(s)
        fasterq-dump "$TMPDIR/$accession" \
          --progress --threads "$nthrd" --temp "$dump_tmp" --mem 15GB \
          --outdir "$TMPDIR"
        conda deactivate
        # Find FASTQ(s) named by accession
        longread="$TMPDIR/${accession}.fastq"
        pe1="$TMPDIR/${accession}_1.fastq"
        pe2="$TMPDIR/${accession}_2.fastq"
        # Handle paired-end reads
        if [[ -f "$pe1" && -f "$pe2" ]]
        then
          # Initialize empty output PE FASTQ files
          : > "$pe1_output"
          : > "$pe2_output"
          append_fastq "$pe1" "$pe1_output"
          append_fastq "$pe2" "$pe2_output"
          outputs=("$pe1_output" "$pe2_output")
        # Handle long reads
        elif [[ -f "$longread" ]]
        then
          # Initialize empty output FASTQ file (if needed)
          [[ -f "$longread_output" ]] || : > "$longread_output"
          append_fastq "$longread" "$longread_output"
          outputs=("$longread_output")
        else
          die "$clean_label: fasterq-dump produced no FASTQ(s) for $accession"
        fi
        # Clean up temp and working directories
        rm -rf -- "${TMPDIR:?}/"*
      done
      ;;
    *)
      die "$clean_label: unsupported read source/method: $source/$method"
      ;;
  esac

  for output in "${outputs[@]}"
  do
    # Load pigz env
    conda activate pigz
    # Verify final FASTQ(s) are compressed
    pigz -t -p "$nthrd" "$output" || die "$clean_label: invalid gzip output: $output"
    conda deactivate
    # Verify final FASTQ(s) are not empty
    if [[ -s "$output" ]]
    then
      log "Successfully retrived reads: $output"
    else
      die "Created empty output: $output"
    fi
  done
}

# Check TSV column names
exp_header=$'# Label\tPretty_Label\tGenome_Source\tGenome_Method\tGenome_Locator\tGenome_Selector\tReads_Source\tReads_Method\tReads_Type\tReads_Locator'
IFS= read -r header < "$in_tsv" # read just first line of TSV
header="${header%$'\r'}" # strip carriage return
[[ "$header" == "$exp_header" ]] || die "Expected TSV header: $exp_header"
# Verify TSV format
awk -F '\t' '
  NR > 1 && $0 !~ /^#/ && NF > 0 && NF != 10 {
    printf "ERROR: line %d has %d TSV columns; expected 10\n", NR, NF > "/dev/stderr"
    bad = 1
  }
  END { exit bad }
' "$in_tsv" || exit 1
# Parse TSV using file separator IFS for "read"
awk -F '\t' -v OFS='\034' '
  NR > 1 && $0 !~ /^#/ && NF > 0 {
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
  }
' "$in_tsv" > "$parsed"

# Iterate through parsed TSV
declare -A seen_labels=()
while IFS=$'\034' read -r label pretty_label genome_source genome_method \
  genome_locator genome_selector reads_source reads_method reads_type \
  reads_locator extra
do
  label="${label%$'\r'}"
  reads_locator="${reads_locator%$'\r'}"
  [[ -n "$label" ]] || continue
  [[ "$label" != \#* ]] || continue
  [[ "$label" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || die "Invalid label: $label"
  [[ -z "${seen_labels[$label]+x}" ]] || die "Duplicate label: $label"
  seen_labels[$label]=1
  [[ -z "${extra:-}" ]] || die "$label: too many TSV columns"
  echo
  log "Starting log for $label [$pretty_label]"
  # Fetch genome assemblies
  if [[ -n "$genome_source" || -n "$genome_method" || -n "$genome_locator" ]]
  then
    log "Getting genome for $label..."
    [[ -n "$genome_source" && -n "$genome_method" && -n "$genome_locator" ]] \
      || die "$label: incomplete genome fields"
    # get_genome "$label" "$genome_source" "$genome_method" "$genome_locator" "$genome_selector"
  fi

  # Strip version from label (same reads for all versions)
  clean_label="${label%%_v[0-9]*}"
  # Skip reads download if single/paired-end read FASTQ(s) exist
  if [[ -f "reads/${clean_label}_$reads_type.fastq.gz" || \
    -f "reads/${clean_label}_1.fastq.gz" && -f "reads/${clean_label}_2.fastq.gz" ]]
  then
    log "Already retrieved $reads_type reads for $clean_label. Skipping..."
  else
    # Fetch assembly reads
    if [[ -n "$reads_source" || -n "$reads_method" || -n "$reads_type" \
      || -n "$reads_locator" ]]
    then
      log "Getting $reads_type reads for $clean_label..."
      [[ -n "$reads_source" && -n "$reads_method" && -n "$reads_type" \
        && -n "$reads_locator" ]] || die "$label: incomplete read fields"
      # get_reads "$clean_label" "$reads_source" "$reads_method" "$reads_type" "$reads_locator"
    fi
  fi
done < "$parsed"

# Clean up working directory
rm -rf -- "$work_dir"

# Log output
echo
echo "Finished. Output:"
tree -A assemblies reads

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo


