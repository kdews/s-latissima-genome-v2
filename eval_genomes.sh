#!/bin/bash
#SBATCH -J eval_genomes
#SBATCH -o %x_logs/%x_%j.log
#SBATCH --time=05:00:00
#SBATCH --mem=3gb

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]; then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
ref_label=corteva_v2 # reference assembly label
ref_asm="assemblies/$ref_label.fa" # reference assembly
ref_chr="assemblies/CHR_$ref_label.fa" # chromosome-only reference assembly

# Help message
show_help () {
    echo "\
Runs BUSCO and QUAST on assemblies found from TSV manifest of genome metadata.

Usage:
  Shell:
    bash $job_name.sh <in_tsv> [path/to/busco.sbatch] [path/to/quast.sbatch]
    bash $job_name.sh <in_tsv> [path/to/busco_compare.sbatch]
  Slurm:
    sbatch $job_name.sh <in_tsv> [path/to/busco.sbatch] [path/to/quast.sbatch]
    sbatch $job_name.sh <in_tsv> [path/to/busco_compare.sbatch]

Requires:
 - BUSCO (https://busco.ezlab.org)
 - QUAST (https://quast.sourceforge.net)
 - R (https://www.r-project.org)"
}
# Catch help flag or missing argument(s)
if [[ $1 == "-h" || $1 == "--help" ]]; then
  show_help
  exit 0
elif (( $# < 2 )); then
  echo "Error: 2 positional arguments required."
  echo
  show_help
  exit 1
fi

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

# Catch errors in pipes
set -o pipefail

# Input
# Assembly metadata
in_tsv="${1:-genomes_metadata.tsv}"
if [[ -f "$in_tsv" ]]; then
  echo "Using TSV: $in_tsv"
  # Assembly labels
  mapfile -t labels < <(grep -v "#" "$in_tsv" | cut -f1)
  mapfile -t pretty_labels < <(grep -v "#" "$in_tsv" | cut -f2)
  echo "Found ${#labels[@]} assembly labels."
else
  echo "Error: input TSV ($in_tsv) not found."
  exit 1
fi
# Scripts to run (all positional arguments after first one)
# script_list=("${@:2}")
declare -A script_list
for script_path in "${@:2}"; do
  script_name="$(basename "$script_path" .sbatch)"
  script_list["$script_name"]="$script_path"
done
echo "Running ${#script_list[@]} scripts: ${script_list[*]}"

# Run BUSCO results comparison visualization script
if [[ -f "${script_list[busco_compare]}" ]]; then
  script="${script_list[busco_compare]}"
  echo "BUSCO COMPARE"
  echo "Submitting BUSCO compare job for assembly versions..."
  for f in "${!labels[@]}"; do
    printf '  %s [%s]\n' "${labels[f]}" "${pretty_labels[f]}"
  done
  echo
  submission=(sbatch "$script" "$in_tsv")
  echo "${submission[*]}"
  "${submission[@]}"
  echo
  
  # Print date and time
  echo
  date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
  date +"$date_fmt"
  exit 0
fi
echo

# Run BUSCO jobs
if [[ -f "${script_list[busco]}" ]]; then
  echo "BUSCO"
  script="${script_list[busco]}"
  echo "Submitting BUSCO jobs..."
  for i in "${!labels[@]}"; do
    label="${labels[i]}"
    pretty_label="${pretty_labels[i]}"
    log="busco_logs/busco_$label.log" # per-assembly log
    assembly="assemblies/$label.fa"
    [[ -f "$assembly" ]] || { echo "Error: assembly ($assembly) not found."; exit 1; }
    echo "Assembly version: $label [$pretty_label]"
    submission=(
      sbatch
      -o "$log"
      "$script"
      "$label"
      "$assembly"
    )
    echo "${submission[*]}"
    "${submission[@]}"
    echo
  done
fi
echo

# Run QUAST jobs
# Declare associative arrays to store grouped labels
declare -A assemblies_by_clean_label
declare -A seen_clean_label
clean_labels=()
if [[ -f "${script_list[quast]}" ]]; then
  echo "QUAST"
  script="${script_list[quast]}"
  if [[ -f "$ref_chr" ]]; then
    echo "Running QUAST with chromosome-only reference for $ref_label: $ref_chr"
  else
    echo "Creating QUAST chromosome-only reference for $ref_label: $ref_chr"
    if [[ -f "$ref_asm" ]]; then
      echo "Using $ref_asm..."
      # Load conda
      cond="$HOME/.conda_for_sbatch.sh"
      if [[ -f "$cond" ]]; then
        source "$cond"
      else
        echo "Error sourcing $cond"
        exit 1
      fi
      # Load seqkit env
      conda activate seqkit
      seqkit version || exit 1
      set -x
      seqkit head -n 31 "$ref_asm" > "$ref_chr"
      set +x
      conda deactivate
      # Load samtools env
      conda activate samtools
      samtools --version | head -n1 || exit 1
      set -x
      samtools faidx "$ref_chr"
      set +x
      conda deactivate
    else
      echo "Error: Reference FASTA not found ($ref_asm)."
    fi
  fi
  echo
  echo "Submitting per-version QUAST jobs..."
  for n in "${!labels[@]}"; do
    label="${labels[n]}"
    pretty_label="${pretty_labels[n]}"
    assembly="assemblies/$label.fa"
    [[ -f "$assembly" ]] \
      || { echo "Error: Assembly ($assembly) not found."; exit 1; }
    log="quast_logs/quast_$label.log" # per-assembly log
    # Preserve first-seen ordering of clean labels
    clean_label="${label%%_v[0-9]*}"
    if [[ -z ${seen_clean_label[$clean_label]+x} ]]; then
      clean_labels+=("$clean_label")
      seen_clean_label["$clean_label"]=1
    fi
    if echo "$assembly" | grep -Pq "lindell|corteva|roscoff|ilvo"; then
      assemblies_by_clean_label[latissima]+="$assembly "
    fi
    # Labels and paths contain no whitespace, use space-delimited values
    assemblies_by_clean_label["$clean_label"]+="$assembly "
    # Check for assembly reads
    # Find FASTQ(s) named by label stripped of version number
    pe1="reads/${clean_label}_1.fastq.gz"
    pe2="reads/${clean_label}_2.fastq.gz"
    long_pb="reads/${clean_label}_pacbio.fastq.gz"
    long_ont="reads/${clean_label}_nanopore.fastq.gz"
    # partition=main # default Slurm partition
    if [[ -f "$long_pb" || -f "$long_ont" || -f "$pe1" && -f "$pe2" ]]; then
      reads_opt=reads
      # read_bytes=0
      # for read_file in "$pe1" "$pe2" "$long_pb" "$long_ont"; do
      #   [[ -f "$read_file" ]] \
      #     && (( read_bytes += $(stat -c '%s' "$read_file") ))
      # done
      # read_gib=$(( (read_bytes + 1024**3 - 1) / 1024**3 ))
      # if (( read_gib == 0 )); then
      #   cpus=12
      #   mem="15g"
      #   walltime="05:00:00"
      # elif (( read_gib <= 50 )); then
      #   cpus=24
      #   mem="64g"
      #   walltime="1-00:00:00"
      # elif (( read_gib <= 125 )); then
      #   partition=largemem
      #   cpus=24
      #   mem="96g"
      #   walltime="2-00:00:00"
      # elif (( read_gib <= 250 )); then
      #   partition=largemem
      #   cpus=32
      #   mem="128g"
      #   walltime="3-00:00:00"
      # else
      #   partition=largemem
      #   cpus=32
      #   mem="192g"
      #   walltime="4-00:00:00"
      # fi
    else
      reads_opt=no_reads
    fi
    echo "Assembly version: $label [$pretty_label]"
    submission=(
      sbatch
      -o "$log"
      # -p "$partition"
      # --cpus-per-task="$cpus"
      # --mem="$mem"
      # --time="$walltime"
      "$script"
      "$label"
      "$ref_chr"
      "$reads_opt"
      "$assembly"
    )
    echo "${submission[*]}"
    "${submission[@]}"
    echo
  done
  echo "Submitting per-assembly QUAST jobs..."
  clean_labels+=(latissima)
  for clean_label in "${clean_labels[@]}"; do
    read -r -a assemblies <<< "${assemblies_by_clean_label[$clean_label]}"
    # Skip groups containing only one assembly
    (( ${#assemblies[@]} > 1 )) || continue
    if [[ "$clean_label" == latissima ]]; then
      echo "Running QUAST on all Saccharina latissima assemblies"
      # Use chromosomal reference for all S. latissima assemblies
      ref="$ref_chr"
    else
      # Use latest assembly version as reference
      last_idx="$(( ${#assemblies[@]} - 1 ))"
      ref="${assemblies[last_idx]}"
      # assemblies=("${assemblies[@]:0:$last_idx}")
    fi
    echo "Assembly: $clean_label (reference for QUAST: $ref)"
    log="quast_logs/quast_$clean_label.log"
    submission=(
      sbatch
      -o "$log"
      -p largemem # use largemem partition
      "$script"
      "$clean_label"
      "$ref"
      no_reads
      "${assemblies[@]}"
    )
    echo "${submission[*]}"
    "${submission[@]}"
    echo
  done
fi
echo

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
