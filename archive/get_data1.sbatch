#!/bin/bash
#SBATCH -J get_data
#SBATCH -o %x_logs/%x_%j.log
#SBATCH --time=05:00:00
#SBATCH --mem=20g

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
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || (( $# < 2 ))
then
  echo "\
Copies from local or downloads from remote all genome FASTAs and reads.
Uses curl or NCBI CLI to download genome FASTAs and SRA Toolkit to download
WGS reads.
Usage: sbatch $job_name.sbatch <>
Requires:
 - curl v8+
 - jq v1.7+
 - NCBI Datasets command-line tools (CLI) v18.31.0 (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/)
 - UnZip v6.00+"
  exit 0
fi

# Input
# GENOMES
# Local path to directory containing assemblies
local_asm_dir='/project2/noujdine_61/kdeweese/latissima/all_s_latissima_assemblies'
# Genomes of Saccharina latissima Lindell (DeWeese et al. 2025)
lindell0_fasta="$local_asm_dir/lindell_v0/SugarkelpHifiasm.raconPolished.fasta.gz"
lindell1_fasta="$local_asm_dir/lindell_v1/SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.fasta"
# Genomes of Saccharina latissima Corteva
corteva1_fasta="$local_asm_dir/corteva_v1/Draft_ONT_TelomereContigs_Saccharina.fasta.gz"
corteva2_fasta="$local_asm_dir/corteva_v2/Slattisima_v0.2_Genome_Polished_111925.fasta"
corteva2p_fasta="$local_asm_dir/corteva_v2/Genome_Polished_013026_Illumina_Polished.fasta.gz"
# Genomes of Saccharina latissima Roscoff (Denoeud et al. 2024)
roscoff1_url='https://phaeoexplorer.sb-roscoff.fr/download/Saccharina-latissima.tar.gz'
roscoff2_url='https://phaeoexplorer.sb-roscoff.fr/download/Saccharina-latissima-v2.tar.gz'
roscoff1_fasta='Saccharina-latissima_FEMALE.fa'
roscoff2_fasta='S-latissima_SLPER63f_v2.0.fa'
# Genome of Ectocarpus sp. v5 (Liu et al. 2024) https://doi.org/10.17617/3.NQDSLW
ec_doi='10.17617/3.NQDSLW' # dataset DOI
ec_api_base='https://edmond.mpg.de' # Dataverse-based Edmond API base URL
ec_api="$ec_api_base/api/datasets/:persistentId/?persistentId=doi:$ec_doi" # Edmond API
# Genome of Saccharina japonica (Li et al. 2025)
# https://figshare.com/articles/dataset/_i_Saccharina_japonica_i_female_sporophyte/27569472
sj_id=27569472 # figshare article ID
sj_api_base='https://api.figshare.com/v2/articles/' # figshare API base URL
sj_api="$sj_api_base/$sj_id/files" # figshare API
# Genome of Saccharina latissima ILVO (Bråtelund et al. 2024)
ilvo_acc='GCA_034768055.1' # NCBI genome accession
# READS
# Local paths to reads used to assemble genomes
# Reads of Saccharina latissima Lindell (DeWeese et al. 2025)
lindell_reads_dir='/project2/noujdine_61/kdeweese/latissima/longreads/ccs'
# Reads of Saccharina latissima Corteva
corteva_reads_dir="$local_asm_dir/corteva_v2/"
# Reads of Saccharina latissima Roscoff (Denoeud et al. 2024)
roscoff_sra=ERR13438528 # NCBI SRA for PacBio female reads
# Reads of Ectocarpus sp. v5 (Liu et al. 2024)
ec_sra1=SRR28880056 # NCBI SRA for ONT male reads
ec_sra2=SRR28831804 # NCBI SRA for ONT female reads
# Reads of Saccharina japonica (Li et al. 2025)
sj_sra=SRR31210069 # NCBI SRA for PacBio reads
# Reads of Saccharina latissima ILVO (Bråtelund et al. 2024)
ilvo_sra=SRR26857257 # NCBI SRA for Illumina reads

# Create output directories
mkdir -p assemblies
mkdir -p reads

### GENOME FASTAS
## Locally hosted genomes
# Saccharina latissima Lindell (DeWeese et al. 2025)
cp "$lindell0_fasta" "$lindell1_fasta" assemblies/
lindell0_fasta="$(basename "$lindell0_fasta")"
lindell1_fasta="$(basename "$lindell1_fasta")"

# Genomes of Saccharina latissima Corteva
cp "$lindell0_fasta" "$lindell1_fasta" assemblies/
corteva1_fasta="$(basename "$corteva1_fasta")"
corteva2_fasta="$(basename "$corteva2_fasta")"
corteva2p_fasta="$(basename "$corteva2p_fasta")"

## Remotely hosted genomes
# Load curl and jq modules
module load jq curl
jq --version || exit 1
curl --version || exit 1

# Saccharina laitissima Roscoff (Denoeud et al. 2024)
# Get tarball names
roscoff1_tar="$(basename "$roscoff1_url")"
roscoff2_tar="$(basename "$roscoff2_url")"
# Download each genome version tarball
set -x
curl -fL "$roscoff1_url" -o "assemblies/$roscoff1_tar"
curl -fL "$roscoff2_url" -o "assemblies/$roscoff2_tar"
set +x
# Extract genome FASTAs from tarballs
set -x
tar -xvf "assemblies/$roscoff1_tar" -C assemblies --wildcards '*'"$roscoff1_fasta"
tar -xvf "assemblies/$roscoff2_tar" -C assemblies --wildcards '*'"$roscoff2_fasta"
set +x
# Clean up
echo "Cleaning up tarballs..."
rm "assemblies/$roscoff1_tar" "assemblies/$roscoff2_tar"

# Saccharina japonica (Li et al. 2025)
# Get FASTA basename and corresponding URL (figshare API)
read -r sj_fasta sj_url < <(
    curl -fsSL "$sj_api" |
    jq -r '
        .[]
        | select(
            .name
            | test("\\.(fa|fna|fasta|fas)(\\.(gz|bgz|bz2|xz|zst))?$"; "i")
        )
        | "\(.name) \(.download_url)"
    '
)
# Download to assemblies/
echo "Downloading $sj_fasta from: $sj_url"
set -x
curl -fL "$sj_url" -o "assemblies/$sj_fasta"
set +x

# Ectocarpus sp. v5 (Liu et al. 2024)
# Get FASTA basename and corresponding URL (Dataverse API)
read -r ec_fasta ec_url < <(
  curl -fsSL "$ec_api" |
  jq -r --arg base "$ec_api_base" '
    .data.latestVersion.files[]
    | .dataFile as $f
    | select(
        $f.filename
        | test("\\.(fa|fna|fasta|fas)(\\.(gz|bgz|bz2|xz|zst))?$"; "i")
      )
    | [$f.filename, ($base + "/api/access/datafile/" + ($f.id | tostring))]
    | @tsv
  '
)
# Download to assemblies/
echo "Downloading $ec_fasta from: $ec_url"
set -x
curl -fL "$ec_url" -o "assemblies/$ec_fasta"
set +x

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

# Saccharina latissima ILVO (Bråtelund et al. 2024)
# Load NCBI CLI env
conda activate ncbi
datasets --version || exit 1
# Download genome accession with NCBI CLI
set -x
datasets download genome accession "$ilvo_acc" --include genome
set +x
conda deactivate
# Load unzip module
module load unzip
unzip -v || exit 1
# Unzip NCBI dataset
set -x
unzip -o ncbi_dataset.zip -x README.md
set +x
echo
# Verify integrity after unzip
md5sum -c md5sum.txt || { echo "Error: NCBI download unsuccessful."; exit 1; }
echo
# Get accession FASTA filename
ilvo_fasta="$(ls "ncbi_dataset/data/$ilvo_acc/"*'.fna')"
ilvo_fasta="$(basename "$ilvo_fasta")"
# Clean up directory
echo "Cleaning NCBI download $ilvo_acc..."
set -x
mv "ncbi_dataset/data/$ilvo_acc/"*'.fna' assemblies/
set -x
rm -rf ncbi_dataset ncbi_dataset.zip md5sum.txt
echo

# Uncompress genome FASTAs
set -x
gunzip "assemblies/$lindell0_fasta" \
  "assemblies/$lindell1_fasta" \
  "assemblies/$corteva1_fasta" \
  "assemblies/$corteva2_fasta" \
  "assemblies/$corteva2p_fasta" \
  "assemblies/$roscoff1_fasta" \
  "assemblies/$roscoff2_fasta" \
  "assemblies/$ilvo_fasta" \
  "assemblies/$ec_fasta"
set +x

# Update paths to downloaded genome FASTAs
lindell0_fasta="assemblies/$(basename "$lindell0_fasta" .gz)"
lindell1_fasta="assemblies/$(basename "$lindell1_fasta" .gz)"
corteva1_fasta="assemblies/$(basename "$corteva1_fasta" .gz)"
corteva2_fasta="assemblies/$(basename "$corteva2_fasta" .gz)"
corteva2p_fasta="assemblies/$(basename "$corteva2p_fasta" .gz)"
roscoff1_fasta="assemblies/$(basename "$roscoff1_fasta" .gz)"
roscoff2_fasta="assemblies/$(basename "$roscoff2_fasta" .gz)"
ilvo_fasta="assemblies/$(basename "$ilvo_fasta" .gz)"
ec_fasta="assemblies/$(basename "$ec_fasta" .gz)"

# Load samtools env
module purge
conda activate samtools
samtools --version || exit 1

# Index genome FASTAs
samtools faidx "$lindell0_fasta"
samtools faidx "$lindell1_fasta"
samtools faidx "$corteva1_fasta"
samtools faidx "$corteva2_fasta"
samtools faidx "$corteva2_fasta"
samtools faidx "$corteva2p_fasta"
samtools faidx "$roscoff1_fasta"
samtools faidx "$roscoff2_fasta"
samtools faidx "$ilvo_fasta"
samtools faidx "$ec_fasta"
conda deactivate

### READS
## Remotely hosted reads
# Load sra-tools env
conda activate sratools
prefetch --version || exit 1
fasterq-dump --version || exit 1

# Function to download FASTQ(s) from SRA(s) and compress FASTQ(s)
get_sra () {
  if [[ "$#" -eq 0 ]]
  then
    echo "Error: No SRA provided."
  else
    # First prefetch SRA, then extract with fasterq-dump
    set -x
    prefetch --max-size 60000000 "$@"
    set +x
    for sra in "$@"
    do
      set -x
      fasterq-dump "./reads/$sra" -F fastq -v -p -t "$TMPDIR" -m 15000MB -O reads
      set +x
      # Clean up
      echo "Cleaning up and compressing SRA download $sra..."
      rm -r "reads/$sra"
      gzip "reads/$sra.fastq"
    done
  fi
}

# Download and compress remote reads
get_sra "$roscoff_sra" "$ec_sra1" "$ec_sra2" "$sj_sra" "$ilvo_sra"
# Combine Ectocarpus reads
set -x
cat "reads/$ec_sra1.fastq.gz" "reads/$ec_sra2.fastq.gz" \
  > reads/ectocarpus_reads.fastq.gz
set +x

## Locally hosted reads
# Combine all Lindell PacBio CCS reads into one FASTQ
set -x
cat "$lindell_reads_dir"/*fastq.gz > reads/lindell_reads.fastq.gz
set +x
# Combine all Corteva ONT reads into one compressed FASTQ
set -x
cat "$corteva_reads_dir"/Slattisima_v0.2_ONT_raw_reads_Part[12].fastq \
  | gzip -c \
  > reads/corteva_reads.fastq.gz
set +x

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo
