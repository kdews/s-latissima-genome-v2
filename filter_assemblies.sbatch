#!/bin/bash
#SBATCH -J filter_assemblies
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=15gb
#SBATCH --time=05:00:00

# Script variables
# Extract job name from script
if [[ -z "${SLURM_JOB_NAME}" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Filters each assembly for contigs longer than provided minimum length (minlen).
Updates seqFile with paths to filtered assembly FASTAs.
Usage: sbatch $job_name.sbatch <seqFile> <minlen>

Required input:
  seqFile   Cactus-formatted seqFile containing a Newick tree and species list
            i.e.:   NEWICK tree
                    name1 path1
                    name2 path2
                    ...
                    nameN pathN
  minlen   minimum length to filter contigs out below

Output:
  <path1>.filt_<minlen>.fasta        size-filtered assembly FASTAs
  <path2>.filt_<minlen>.fasta
  ...
  <pathN>.filt_<minlen>.fasta

Requires:
 - seqkit (https://bioinf.shenwei.me/seqkit/)"
  exit 0
fi

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"

# Load conda
cond="$HOME/.conda_for_sbatch.sh"
if [[ -f "$cond" ]]
then
  source "$cond"
else
  echo "Error sourcing $cond"
  exit 1
fi

# Load seqkit env
module purge
conda activate seqkit
seqkit version || exit 1

# Input
# Cactus seqFile
# seqFile=s_latissima_prog_align.txt
seqFile="$1"
# Minimum length to filter assemblies, e.g., >50kb
# minlen=50000
minlen="$2"
if [[ -f "$seqFile" ]]
then
  echo "Using seqFile: $seqFile"
  # Assembly labels
  mapfile -t lbl_list < <(grep -v "(" "$seqFile" | awk -F '\t' '{print $1}')
  # Paths to assembly FASTAs
  mapfile -t assembly_list < <(grep -v "(" "$seqFile" | awk -F '\t' '{print $2}')
  echo "Found ${#assembly_list[@]} assemblies."
else
  echo "Error: seqFile ($seqFile) not found."
  exit 1
fi
if [[ -n "$minlen" ]]
then
  echo "Filtering for contigs >$minlen bp"
else
  echo "Error: minimum length not set."
  exit 1
fi

# Output
# Create new filtered seqFile
filt_seqFile="${seqFile%%.*}.filt_$minlen.txt"
cp "$seqFile" "$filt_seqFile"

# Filter each assembly for contigs >minimum length
filt_asms=()
for i in "${!assembly_list[@]}"
do
  lbl="${lbl_list[i]}"
  asm="${assembly_list[i]}"
  asm_no_ext="${asm%%.fasta}"
  asm_no_ext="${asm_no_ext%%.fa}"
  filt_asm="$asm_no_ext.filt_$minlen.fasta"
  echo "Updating filtered seqFile ($filt_seqFile) with filtered $lbl assembly..."
  sed -i "s,$asm,$filt_asm,g" "$filt_seqFile"
  if [[ -f "$asm" ]]
  then
    echo "Filtering $lbl assembly: $asm"
    # seqkit seq -m "$minlen" "$asm" > "$filt_asm"
  else
    echo "Error: $lbl assembly not found at path: $asm"
    exit 1
  fi
done

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
