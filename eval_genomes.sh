#!/bin/bash
#SBATCH -J eval_genomes
#SBATCH -o %x_logs/%x_%j.log
#SBATCH --time=05:00:00
#SBATCH --mem=3gb

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi

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
if [[ $1 == "-h" || $1 == "--help" ]]
then
  show_help
  exit 0
elif (( $# < 2 ))
then
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

# Input
# Assembly metadata
in_tsv="${1:-genomes_metadata.tsv}"
if [[ -f "$in_tsv" ]]
then
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
script_list=("${@:2}")
echo "Running ${#script_list[@]} scripts: ${script_list[*]}"

# Arrays to store grouped labels
declare -A assemblies_by_clean_label
declare -A seen_clean_label
clean_labels=()
for script in "${script_list[@]}"
do
  echo
  echo
  script_no_ext="$(basename "$script")"
  script_no_ext="${script_no_ext%%.*}"
  if [[ -f "$script" ]]
  then
    if [[ "$script_no_ext" == busco_compare ]]
    then
      echo "Submitting BUSCO compare job for assembly versions..."
      for f in "${!labels[@]}"
      do
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
    else
      if [[ "$script_no_ext" == busco ]]
      then
        echo "Submitting BUSCO jobs..."
        for i in "${!labels[@]}"
        do
          label="${labels[i]}"
          pretty_label="${pretty_labels[i]}"
          log="${script_no_ext}_logs/${script_no_ext}_$label.log" # per-assembly log
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
      elif [[ "$script_no_ext" == quast ]]
      then
        echo "Submitting per-version QUAST jobs..."
        for n in "${!labels[@]}"
        do
          label="${labels[n]}"
          pretty_label="${pretty_labels[n]}"
          assembly="assemblies/$label.fa"
          [[ -f "$assembly" ]] \
            || { echo "Error: assembly ($assembly) not found."; exit 1; }
          log="${script_no_ext}_logs/${script_no_ext}_$label.log" # per-assembly log
          partition=main # default Slurm partition
          # Preserve first-seen ordering of clean labels
          clean_label="${label%%_v[0-9]*}"
          if [[ -z ${seen_clean_label[$clean_label]+x} ]]
          then
            clean_labels+=("$clean_label")
            seen_clean_label["$clean_label"]=1
          fi
          # Labels and paths contain no whitespace, use space-delimited values
          assemblies_by_clean_label["$clean_label"]+="$assembly "
          # Check for assembly reads
          # Find FASTQ(s) named by label stripped of version number
          pe1="reads/${clean_label}_1.fastq.gz"
          pe2="reads/${clean_label}_2.fastq.gz"
          long_pb="reads/${clean_label}_pacbio.fastq.gz"
          long_ont="reads/${clean_label}_nanopore.fastq.gz"
          if [[ -f "$long_pb" || -f "$long_ont" || -f "$pe1" && -f "$pe2" ]]
          then
            reads_opt=reads
            read_bytes=0
            for read_file in "$pe1" "$pe2" "$long_pb" "$long_ont"
            do
              [[ -f "$read_file" ]] \
                && (( read_bytes += $(stat -c '%s' "$read_file") ))
            done
            read_gib=$(( (read_bytes + 1024**3 - 1) / 1024**3 ))
            if (( read_gib == 0 ))
            then
              cpus=12
              mem="15G"
              walltime="05:00:00"
            elif (( read_gib <= 50 ))
            then
              cpus=24
              mem="64G"
              walltime="1-00:00:00"
            elif (( read_gib <= 125 ))
            then
              partition=largemem
              cpus=24
              mem="96G"
              walltime="2-00:00:00"
            elif (( read_gib <= 250 ))
            then
              partition=largemem
              cpus=32
              mem="128G"
              walltime="3-00:00:00"
            else
              partition=largemem
              cpus=32
              mem="192G"
              walltime="4-00:00:00"
            fi
          else
            reads_opt=no_reads
          fi
          echo "Assembly version: $label [$pretty_label]"
          submission=(
            sbatch
            -o "$log"
            -p "$partition"
            --cpus-per-task="$cpus"
            --mem="$mem"
            --time="$walltime"
            "$script"
            "$label"
            "$reads_opt"
            "$assembly"
          )
          echo "${submission[*]}"
          "${submission[@]}"
          echo
        done
      fi
    fi
  else
    echo "Error: script $script not found."
    exit 1
  fi
done
echo

# Submit one QUAST job per clean label
for script in "${script_list[@]}"
do
  script_no_ext="$(basename "$script")"
  script_no_ext="${script_no_ext%%.*}"
  if [[ "$script_no_ext" == quast ]]
  then
    echo "Submitting per-assembly QUAST jobs..."
    for clean_label in "${clean_labels[@]}"
    do
      read -r -a assemblies <<< "${assemblies_by_clean_label[$clean_label]}"
      # Skip groups containing only one assembly
      (( ${#assemblies[@]} > 1 )) || continue
      echo "Assembly: $clean_label"
      log="${script_no_ext}_logs/${script_no_ext}_$clean_label.log"
      submission=(
        sbatch
        -o "$log"
        "$script"
        "$clean_label"
        no_reads
        "${assemblies[@]}"
      )
      echo "${submission[*]}"
      "${submission[@]}"
      echo
    done
  fi
done

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
