#!/bin/bash

# Script variables
# Extract job name from script
if [[ -z "${SLURM_JOB_NAME}" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
rscript="parse_ena_tsv.R" # R script to run
scripts_dir="s-latissima-genome-v2/" # scripts directory
outdir="seq_divergence" # output directory
# Input ENA TSV from https://www.ebi.ac.uk/ena/browser/view/PRJEB72149
ena_tsv="dna_reads_PRJEB72149.tsv"
# Prepend scripts directory to Rscript and TSV (if exists)
if [[ -d "$scripts_dir" ]]
then
  rscript="${scripts_dir}$rscript"
  ena_tsv="${scripts_dir}$ena_tsv"
fi
out_tsv="$outdir/parsed_$(basename "$ena_tsv")" # parsed output TSV

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Parses TSV from ENA and downloads FASTQs via FTP protocol.
Usage: bash $job_name.sbatch

Output:
  └── $outdir
      ├── $(basename $out_tsv)      parsed ENA TSV
      └── <ENA_ID>.fastq.gz         locally downloaded FASTQs

Requires:
 - R v4.6.0 (https://www.r-project.org/)"
  exit 0
fi

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"

# Create output directory
mkdir -p "$outdir"

# Load modules for Rscript
module purge
module load r/4.6.0
Rscript --version || exit 1

# Run Rscript to parse ENA TSV
set -x
Rscript "$rscript" "$ena_tsv" "$outdir"
set +x

# Unload modules
module purge

mapfile -t ftps < <(cut -f1 "$out_tsv")
mapfile -t md5s < <(cut -f2 "$out_tsv")
mapfile -t fastqs < <(cut -f3 "$out_tsv")

for i in "${!ftps[@]}"
do
  md5="${md5[i]}"
  ftp="${ftps[i]}"
  fastq="${fastqs[i]}"
  set -x
  curl -o "$fastq" "$ftp"
  set +x
  echo "$md5  file.txt" | md5sum -c
done

# Print date and time
echo
date +"$date_fmt"

