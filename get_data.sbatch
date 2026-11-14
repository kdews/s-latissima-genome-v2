#!/bin/bash
#SBATCH -J get_data
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p oneweek
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=185g

# Script variables
# Extract job name from script
if [[ -n ${SLURM_JOB_NAME:-} ]]; then
  job_name="$SLURM_JOB_NAME"
elif [[ -n ${BASH_SOURCE[0]:-} && ${BASH_SOURCE[0]} != "bash" ]]; then
  job_name="$(basename "${BASH_SOURCE[0]}")"
  job_name="${job_name%%.*}"
else
  job_name="get_data"
fi
# Max Plank's Edmond Dataverse API
DATAVERSE_API_BASE='https://edmond.mpg.de'

# Help message
show_help () {
  echo "\
Using a TSV manifest of genome metadata, copies from local or downloads remote
genome FASTAs, optional GFF annotations, and reads used to construct each
assembly. Separately prepares sex-specific Illumina read lists for downstream
SDR and depth analyses.

Usage: sbatch $job_name.sbatch <in_tsv>

Required input:
  in_tsv      TSV with 20 named columns (see below), with commented (#) header

Input TSV columns:
  <Column_Name>           <Example>
  1. Label                corteva_v2
  2. Pretty_Name          Saccharina latissima (US v2)
  3. Genome_Source        local
  4. Genome_Method        file
  5. Genome_Locator       /project2/noujdine_61/kdeweese/latissima/all_s_latissima_assemblies/corteva_v2/Slattisima_v0.2_Genome_Polished_111925.fasta
  6. Genome_Selector      [blank]*
  7. Assembly_Reads_Source          local
  8. Assembly_Reads_Method          glob
  9. Assembly_Reads_Platform        nanopore
  10. Assembly_Reads_Locator        /path/to/assembly_reads/*.fastq
  11. Analysis_Reads_Platform       illumina
  12. Analysis_Reads_Method         script
  13. Analysis_Reads_Locator        saccharina_reads_by_sex.R
  14. Analysis_Reads_Female_Locator [blank for method=script]
  15. Analysis_Reads_Male_Locator   [blank for method=script]
  16. Annotation_GFF                annotation.gff3.gz**
  17. Sample_ID                     SL-CT1-FG3
  18. Region                        US
  19. Sex                           V
  20. Ploidy                        1

* Assembly and analysis read fields are independently optional. Blank assembly
  read fields do not prevent genome or analysis-read retrieval.

** Optional. For local assemblies, provide the full GFF path. For Figshare,
Dataverse, or url_tar sources, provide the GFF filename. For direct URLs,
provide either the full GFF URL or a filename beside the genome URL.

Supported genome sources:
  file        local FASTA (plain or gzip-compressed)
  url         direct FASTA URL (plain or gzip-compressed)
  url_tar     tar archive URL; Genome_Selector is the member basename
  figshare    Figshare article ID; selector is a FASTA filename or regex
  dataverse   DOI; selector is a FASTA filename or regex
  ncbi        NCBI assembly accession

Supported assembly-read sources:
  glob        local shell glob; matching FASTQ/FASTQ.GZ files are concatenated
  sra         one or more comma-separated SRA run accessions

Supported analysis-read methods:
  script      run Analysis_Reads_Locator to generate <label>_female.list and
              <label>_male.list; FASTQs referenced by the lists are not copied
  sra         download paired runs from Analysis_Reads_Female_Locator and
              Analysis_Reads_Male_Locator; write paths to sex-specific lists

Output:
  assemblies/             # directory containing all assemblies
  ├── <label>.fa          # renamed assembly FASTA (uncompressed)
  └── <label>.fa.fai      # FASTA index
  reads/                              # reads used to construct assemblies
  ├── <assembly-series>_<platform>.fastq.gz
  ├── <assembly-series>_illumina_1.fastq.gz
  └── <assembly-series>_illumina_2.fastq.gz
  read_lists/                         # analysis-read paths
  ├── <label>_female.list
  └── <label>_male.list
  annotations/            # directory containing available annotations
  └── <label>.gff3        # renamed annotation GFF (uncompressed)

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

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

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
# Allow SRA downloads up to 90% of requested memory
max_sra_mem=$(( SLURM_MEM_PER_NODE * 90 / 100 ))

# Input
in_tsv="${1:-s-latissima-genome-v2/genomes_metadata.tsv}"
[[ -r $in_tsv ]] || { echo "Cannot read input TSV: $in_tsv"; exit 1; }
scripts_dir="s-latissima-genome-v2"
[[ -d "$scripts_dir" ]] || scripts_dir="."

# Output
work_dir="${job_name}_tmp" # temporary working directory
mkdir -p "$work_dir" assemblies reads read_lists annotations # create temp/output dirs
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

# Select one GFF/GFF3 record from a Figshare or Dataverse JSON response
select_gff () {
  local selector="$1"
  jq -er --arg selector "$selector" '
    map(select(.name | test("\\.(gff|gff3)(\\.(gz|bgz))?$"; "i")))
    | map(select((.name == $selector) or (.name | test($selector))))
    | if length == 1 then .[0]
      elif length == 0 then error("no matching GFF")
      else error("multiple matching GFFs; use an exact Annotation_GFF filename")
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

# Copy/decompress an annotation to a normalized, uncompressed output path
norm_annotation () {
  local input="$1" output="$2"
  log "Input: $input"
  log "Output: $output"
  conda activate pigz
  if pigz -t -p "$nthrd" "$input" 2>/dev/null
  then
    pigz -cd -p "$nthrd" -- "$input" > "$output"
  else
    cp -- "$input" "$output"
  fi
  conda deactivate
  [[ -s "$output" ]] || die "Created empty annotation output: $output"
  log "Successfully retrieved annotation: $output"
}

get_annotation () {
  local label="$1" source="$2" method="$3" locator="$4" annotation_gff="$5"
  local output="annotations/$label.gff3"
  local download="$work_dir/$label.annotation"
  local json record url member

  log "$label: preparing annotation ($source/$method)"

  # An explicit local path or URL overrides the assembly source. This permits,
  # for example, a local GFF paired with a remotely retrieved assembly.
  if [[ "$annotation_gff" == /* ]]
  then
    [[ -r "$annotation_gff" ]] \
      || die "$label: cannot read annotation: $annotation_gff"
    norm_annotation "$annotation_gff" "$output"
    return
  elif [[ "$annotation_gff" =~ ^https?:// ]]
  then
    module load usc curl
    curl -fL --retry 3 --output "$download" "$annotation_gff"
    module purge
    norm_annotation "$download" "$output"
    return
  fi

  case "$source:$method" in
    local:file)
      [[ -r "$annotation_gff" ]] \
        || die "$label: cannot read annotation: $annotation_gff"
      norm_annotation "$annotation_gff" "$output"
      ;;
    remote:url)
      url="${locator%/*}/$annotation_gff"
      module load usc curl
      curl -fL --retry 3 --output "$download" "$url"
      module purge
      norm_annotation "$download" "$output"
      ;;
    remote:url_tar)
      module load usc curl
      curl -fL --retry 3 --output "$download.tar" "$locator"
      module purge
      tar -tf "$download.tar" \
        | awk -v name="$annotation_gff" 'BEGIN { FS="/" } $NF == name' \
        > "$work_dir/annotation_members"
      mapfile -t members < "$work_dir/annotation_members"
      (( ${#members[@]} == 1 )) \
        || die "$label: expected one GFF named $annotation_gff; found ${#members[@]}"
      member="${members[0]}"
      tar -xOf "$download.tar" "$member" > "$download"
      norm_annotation "$download" "$output"
      ;;
    remote:figshare)
      module load usc curl jq
      json="$(curl -fsSL --retry 3 "https://api.figshare.com/v2/articles/$locator/files")"
      record="$(jq -ce '[.[] | {name: .name, url: .download_url}]' <<<"$json" \
        | select_gff "$annotation_gff")" \
        || die "$label: could not select one Figshare GFF"
      url="$(jq -r .url <<<"$record")"
      curl -fL --retry 3 --output "$download" "$url"
      module purge
      norm_annotation "$download" "$output"
      ;;
    remote:dataverse)
      module load usc curl jq
      json="$(curl -fsSL --retry 3 \
        "$DATAVERSE_API_BASE/api/datasets/:persistentId/?persistentId=doi:$locator")"
      record="$(jq -ce --arg base "$DATAVERSE_API_BASE" '
        [.data.latestVersion.files[].dataFile
          | {name: .filename, url: ($base + "/api/access/datafile/" + (.id | tostring))}]
      ' <<<"$json" | select_gff "$annotation_gff")" \
        || die "$label: could not select one Dataverse GFF"
      url="$(jq -r .url <<<"$record")"
      curl -fL --retry 3 --output "$download" "$url"
      module purge
      norm_annotation "$download" "$output"
      ;;
    *)
      die "$label: annotations are unsupported for source/method: $source/$method"
      ;;
  esac
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

get_assembly_reads () {
  local clean_label="$1" source="$2" method="$3" platform="$4" locator="$5"
  local single_output="reads/${clean_label}_${platform}.fastq.gz"
  local pe1_output="reads/${clean_label}_${platform}_1.fastq.gz"
  local pe2_output="reads/${clean_label}_${platform}_2.fastq.gz"
  local accession fastq dump_tmp single pe1 pe2 locator_r1 locator_r2 extra
  local -a inputs inputs_r1 inputs_r2 accessions outputs

  [[ "$source" == local || "$source" == remote ]] \
    || die "$clean_label: Assembly_Reads_Source must be local or remote"
  log "$clean_label: preparing $platform assembly reads ($source/$method)"

  case "$source:$method" in
    local:glob)
      if [[ "$platform" == illumina ]]
      then
        IFS=',' read -r locator_r1 locator_r2 extra <<<"$locator"
        [[ -n "$locator_r1" && -n "$locator_r2" && -z "${extra:-}" ]] \
          || die "$clean_label: local Illumina assembly reads require two locator globs (R1,R2)"
        compgen -G "$locator_r1" | sort > "$work_dir/${clean_label}_assembly_r1" || true
        compgen -G "$locator_r2" | sort > "$work_dir/${clean_label}_assembly_r2" || true
        mapfile -t inputs_r1 < "$work_dir/${clean_label}_assembly_r1"
        mapfile -t inputs_r2 < "$work_dir/${clean_label}_assembly_r2"
        (( ${#inputs_r1[@]} > 0 && ${#inputs_r1[@]} == ${#inputs_r2[@]} )) \
          || die "$clean_label: Illumina R1/R2 globs must match the same nonzero number of files"
        : > "$pe1_output"
        : > "$pe2_output"
        for fastq in "${inputs_r1[@]}"; do append_fastq "$fastq" "$pe1_output"; done
        for fastq in "${inputs_r2[@]}"; do append_fastq "$fastq" "$pe2_output"; done
        outputs=("$pe1_output" "$pe2_output")
      else
        compgen -G "$locator" | sort > "$work_dir/${clean_label}_assembly_inputs" || true
        mapfile -t inputs < "$work_dir/${clean_label}_assembly_inputs"
        (( ${#inputs[@]} > 0 )) \
          || die "$clean_label: glob matched no files: $locator"
        : > "$single_output"
        for fastq in "${inputs[@]}"; do append_fastq "$fastq" "$single_output"; done
        outputs=("$single_output")
      fi
      ;;
    remote:sra)
      IFS=',' read -r -a accessions <<<"$locator"
      if [[ "$platform" == illumina ]]
      then
        : > "$pe1_output"
        : > "$pe2_output"
        outputs=("$pe1_output" "$pe2_output")
      else
        : > "$single_output"
        outputs=("$single_output")
      fi

      for accession in "${accessions[@]}"
      do
        accession="${accession//[[:space:]]/}"
        [[ "$accession" =~ ^[SED]RR[0-9]+$ ]] \
          || die "$clean_label: invalid assembly-read SRA accession: $accession"
        conda activate sratools
        prefetch --progress --max-size "${max_sra_mem}M" \
          --output-directory "$TMPDIR" "$accession"
        dump_tmp="$TMPDIR/dump_tmp_$accession"
        mkdir -p "$dump_tmp"
        if [[ "$platform" == illumina ]]
        then
          fasterq-dump "$TMPDIR/$accession" \
            --progress --threads "$nthrd" --temp "$dump_tmp" --mem 15GB \
            --outdir "$TMPDIR" --split-files
        else
          fasterq-dump "$TMPDIR/$accession" \
            --progress --threads "$nthrd" --temp "$dump_tmp" --mem 15GB \
            --outdir "$TMPDIR"
        fi
        conda deactivate

        if [[ "$platform" == illumina ]]
        then
          pe1="$TMPDIR/${accession}_1.fastq"
          pe2="$TMPDIR/${accession}_2.fastq"
          [[ -s "$pe1" && -s "$pe2" ]] \
            || die "$clean_label: expected paired Illumina FASTQs for $accession"
          append_fastq "$pe1" "$pe1_output"
          append_fastq "$pe2" "$pe2_output"
          rm -rf -- "${TMPDIR:?}/$accession" "$dump_tmp" "$pe1" "$pe2"
        else
          single="$TMPDIR/${accession}.fastq"
          [[ -s "$single" ]] \
            || die "$clean_label: expected single-end $platform FASTQ for $accession"
          append_fastq "$single" "$single_output"
          rm -rf -- "${TMPDIR:?}/$accession" "$dump_tmp" "$single"
        fi
      done
      ;;
    *)
      die "$clean_label: unsupported assembly-read source/method: $source/$method"
      ;;
  esac

  for output in "${outputs[@]}"
  do
    [[ -s "$output" ]] || die "Created empty output: $output"
    # Validate compression
    conda activate pigz
    pigz -t -p "$nthrd" "$output" \
      || die "$clean_label: invalid gzip output: $output"
    conda deactivate
    log "Successfully retrieved assembly reads: $output"
  done
}

# Run a dataset-specific script that writes paths to existing analysis FASTQs.
# Interface: SCRIPT LABEL SAMPLE_ID SCRIPTS_DIR OUTPUT_DIR
make_analysis_lists () {
  local label="$1" sample_id="$2" locator="$3"
  local script="$locator"
  local female_list="read_lists/${label}_female.list"
  local male_list="read_lists/${label}_male.list"

  if [[ ! -r "$script" && -r "$scripts_dir/$script" ]]
  then
    script="$scripts_dir/$script"
  fi
  [[ -r "$script" ]] || die "$label: cannot read analysis-read script: $locator"

  if [[ -s "$female_list" && -s "$male_list" ]]
  then
    log "$label: analysis read lists already exist. Skipping..."
    return
  fi

  module --latest load r
  Rscript --version || die "$label: Rscript is unavailable"
  Rscript "$script" "$label" "$sample_id" "$scripts_dir" "read_lists"
  module purge

  [[ -s "$female_list" ]] || die "$label: script did not create $female_list"
  [[ -s "$male_list" ]] || die "$label: script did not create $male_list"
  log "$label: created analysis read lists without copying local FASTQs"
}

# Download paired-end Illumina SRA runs without merging biological samples.
# Each accession retains its own R1/R2 files and is added to the sex-specific
# analysis read list.
get_analysis_sra () {
  local label="$1" sx="$2" locator="$3"
  local list_file="read_lists/${label}_${sx}.list"
  local accession dump_tmp pe1 pe2 out1 out2
  local -a accessions

  IFS=',' read -r -a accessions <<<"$locator"
  : > "$list_file"

  for accession in "${accessions[@]}"
  do
    accession="${accession//[[:space:]]/}"
    [[ "$accession" =~ ^[SED]RR[0-9]+$ ]] \
      || die "$label: invalid $sx SRA run accession: $accession"
    
    # Save accession metadata to CSV in reads/
    conda activate sratools
    NQUIRE_HELPER=wget efetch -db sra -id "$accession" -format runinfo \
      > "reads/${accession}_runinfo.csv" \
      || die "Could not get metatdata for $accession"
    conda deactivate

    out1="reads/${accession}_1.fastq.gz"
    out2="reads/${accession}_2.fastq.gz"
    if [[ -s "$out1" && -s "$out2" ]]
    then
      log "$label: $accession already downloaded. Skipping..."
    elif [[ -e "$out1" || -e "$out2" ]]
    then
      die "$label: incomplete existing paired output for $accession"
    else
      conda activate sratools
      prefetch --progress --max-size "${max_sra_mem}M" \
        --output-directory "$TMPDIR" "$accession"
      dump_tmp="reads/dump_tmp_$accession"
      mkdir -p "$dump_tmp"
      # Larger analysis FASTQs ouput directly to reads/ directory
      fasterq-dump "$TMPDIR/$accession" \
        --progress --threads "$nthrd" --temp "$dump_tmp" --mem 15GB \
        --outdir reads --split-files
      conda deactivate
      # Check for fasterq-dump output FASTQs
      pe1="reads/${accession}_1.fastq"
      pe2="reads/${accession}_2.fastq"
      [[ -s "$pe1" && -s "$pe2" ]] \
        || die "$label: expected paired FASTQs for $sx run $accession"
      rm -rf -- "${TMPDIR:?}/$accession" "$dump_tmp" # clean up TMPDIR
      # Compress FASTQs
      conda activate pigz
      pigz -c -p "$nthrd" "$pe1" > "$out1" \
        || die "$label: compression failed for $pe1"
      pigz -c -p "$nthrd" "$pe2" > "$out2" \
        || die "$label: compression failed for $pe2"
      # Validate compression
      pigz -t -p "$nthrd" "$out1" "$out2" \
        || die "$label: invalid gzip output for $accession"
      conda deactivate
      rm -f -- "$pe1" "$pe2" # remove uncompressed FASTQs
    fi

    printf '%s\n%s\n' "$(realpath -e "$out1")" "$(realpath -e "$out2")" \
      >> "$list_file"
  done

  [[ -s "$list_file" ]] || die "$label: created empty list: $list_file"
  log "$label: wrote $sx analysis-read paths to $list_file"
}

# Check TSV column names
exp_header=$'# Label\tPretty_Label\tGenome_Source\tGenome_Method\tGenome_Locator\tGenome_Selector\tAssembly_Reads_Source\tAssembly_Reads_Method\tAssembly_Reads_Platform\tAssembly_Reads_Locator\tAnalysis_Reads_Platform\tAnalysis_Reads_Method\tAnalysis_Reads_Locator\tAnalysis_Reads_Female_Locator\tAnalysis_Reads_Male_Locator\tAnnotation_GFF\tSample_ID\tRegion\tSex\tPloidy'
IFS= read -r header < "$in_tsv" # read just first line of TSV
header="${header%$'\r'}" # strip carriage return
[[ "$header" == "$exp_header" ]] || die "Expected TSV header: $exp_header"
# Verify TSV format
awk -F '\t' '
  NR > 1 && $0 !~ /^#/ && NF > 0 && NF != 20 {
    printf "ERROR: line %d has %d TSV columns; expected 20\n", NR, NF > "/dev/stderr"
    bad = 1
  }
  END { exit bad }
' "$in_tsv" || exit 1
# Parse TSV using file separator IFS for "read"
awk -F '\t' -v OFS='\034' '
  NR > 1 && $0 !~ /^#/ && NF > 0 {
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20
  }
' "$in_tsv" > "$parsed"

# Iterate through parsed TSV
declare -A seen_labels=()
while IFS=$'\034' read -r label pretty_label genome_source genome_method \
  genome_locator genome_selector assembly_reads_source assembly_reads_method \
  assembly_reads_platform assembly_reads_locator analysis_reads_platform \
  analysis_reads_method analysis_reads_locator analysis_reads_female_locator \
  analysis_reads_male_locator annotation_gff sample_id region sex ploidy extra
do
  label="${label%$'\r'}"
  assembly_reads_locator="${assembly_reads_locator%$'\r'}"
  analysis_reads_locator="${analysis_reads_locator%$'\r'}"
  analysis_reads_female_locator="${analysis_reads_female_locator%$'\r'}"
  analysis_reads_male_locator="${analysis_reads_male_locator%$'\r'}"
  annotation_gff="${annotation_gff%$'\r'}"
  ploidy="${ploidy%$'\r'}"
  [[ -n "$label" ]] || continue
  [[ "$label" != \#* ]] || continue
  [[ "$label" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || die "Invalid label: $label"
  [[ -z "${seen_labels[$label]+x}" ]] || die "Duplicate label: $label"
  seen_labels[$label]=1
  [[ -z "${extra:-}" ]] || die "$label: too many TSV columns"
  echo
  log "Starting log for $label [$pretty_label]"
  # Skip assembly download if FASTA exists in assemblies/
  if [[ -f "assemblies/$label.fa" ]]
  then
    log "Already retrieved assembly for $label. Skipping..."
  else
    # Fetch genome assemblies
    if [[ -n "$genome_source" || -n "$genome_method" || -n "$genome_locator" ]]
    then
      log "Getting genome for $label..."
      [[ -n "$genome_source" && -n "$genome_method" && -n "$genome_locator" ]] \
        || die "$label: incomplete genome fields"
      get_genome "$label" "$genome_source" "$genome_method" "$genome_locator" "$genome_selector"
    fi
  fi

  # Fetch annotation when one is specified and normalized output is absent
  if [[ -n "$annotation_gff" ]]
  then
    if [[ -f "annotations/$label.gff3" ]]
    then
      log "Already retrieved annotation for $label. Skipping..."
    else
      log "Getting annotation for $label..."
      get_annotation "$label" "$genome_source" "$genome_method" \
        "$genome_locator" "$annotation_gff"
    fi
  fi

  # Assembly-read metadata must be complete or entirely blank.
  assembly_read_fields=("$assembly_reads_source" "$assembly_reads_method" \
    "$assembly_reads_platform" "$assembly_reads_locator")
  n_assembly_read_fields=0
  for field in "${assembly_read_fields[@]}"
  do
    [[ -n "$field" ]] && (( n_assembly_read_fields += 1 ))
  done
  (( n_assembly_read_fields == 0 || n_assembly_read_fields == 4 )) \
    || die "$label: all four Assembly_Reads fields must be filled or all blank"

  if (( n_assembly_read_fields == 0 ))
  then
    log "$label: no assembly reads specified. Skipping..."
  else
    # Strip version from label (same reads for all versions)
    clean_label="${label%%_v[0-9]*}"
    if [[ "$assembly_reads_platform" == illumina ]]
    then
      assembly_outputs=(
        "reads/${clean_label}_illumina_1.fastq.gz"
        "reads/${clean_label}_illumina_2.fastq.gz"
      )
    else
      assembly_outputs=("reads/${clean_label}_${assembly_reads_platform}.fastq.gz")
    fi
    if [[ -s "${assembly_outputs[0]}" && \
          ( ${#assembly_outputs[@]} -eq 1 || -s "${assembly_outputs[1]}" ) ]]
    then
      log "Already retrieved $assembly_reads_platform assembly reads for $clean_label. Skipping..."
    else
      log "Getting $assembly_reads_platform assembly reads for $clean_label..."
      get_assembly_reads "$clean_label" "$assembly_reads_source" \
        "$assembly_reads_method" "$assembly_reads_platform" \
        "$assembly_reads_locator"
    fi
  fi

  # Analysis reads are currently paired-end Illumina. Local script methods
  # generate lists without copying FASTQs; SRA retains sex assignments.
  if [[ -z "$analysis_reads_platform" && -z "$analysis_reads_method" && \
        -z "$analysis_reads_locator" && -z "$analysis_reads_female_locator" && \
        -z "$analysis_reads_male_locator" ]]
  then
    log "$label: no analysis reads specified. Skipping..."
  else
    [[ "$analysis_reads_platform" == illumina ]] \
      || die "$label: Analysis_Reads_Platform must be illumina"
    case "$analysis_reads_method" in
      script)
        [[ -n "$analysis_reads_locator" && \
           -z "$analysis_reads_female_locator" && \
           -z "$analysis_reads_male_locator" ]] \
          || die "$label: analysis method=script requires only Analysis_Reads_Locator"
        make_analysis_lists "$label" "$sample_id" "$analysis_reads_locator"
        ;;
      sra)
        [[ -z "$analysis_reads_locator" && \
           -n "$analysis_reads_female_locator" && \
           -n "$analysis_reads_male_locator" ]] \
          || die "$label: analysis method=sra requires female and male locators"
        get_analysis_sra "$label" female "$analysis_reads_female_locator"
        get_analysis_sra "$label" male "$analysis_reads_male_locator"
        ;;
      *)
        die "$label: unsupported analysis-read method: $analysis_reads_method"
        ;;
    esac
  fi
done < "$parsed"

# Clean up working directory
rm -rf -- "$work_dir"

# Log output
echo
echo "Finished. Output:"
tree -n --charset=UTF-8 assemblies reads read_lists annotations

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo
