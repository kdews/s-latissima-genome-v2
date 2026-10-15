#!/bin/bash
#SBATCH -J depth_sdr
#SBATCH -o %x_logs/%x_%a_%j.log
#SBATCH -p gpu
#SBATCH --time=05:00:00
#SBATCH --mem=65g
#SBATCH --cpus-per-task=8
#SBATCH -a 1-20

# Script variables
# Extract job name from script
if [[ -n ${SLURM_JOB_NAME:-} ]]; then
  job_name="$SLURM_JOB_NAME"
elif [[ -n ${BASH_SOURCE[0]:-} && ${BASH_SOURCE[0]} != "bash" ]]; then
  job_name="$(basename "${BASH_SOURCE[0]}")"
  job_name="${job_name%%.*}"
else
  job_name="depth_sdr"
fi
main_outdir="${2:-sdr_id_results}" # main output directory
depth_outdir="$main_outdir/${job_name}_results" # job type output directory
depth_reads_file="os_ss.txt" # paths to same and opposite sex reads
depth_tab_file="os_ss_table.tsv" # table of sample IDs

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
Identifies SDR with sex-biased depth analysis from Barrera-Redondo et al. (2025).
Reads must be include a roughly equal number of reads from opposite sex (os)
and same sex (ss) with respect to the assembly.

Usage: sbatch $job_name.sbatch <genome.fasta(.gz)> [main_outdir]

Input:
  genome.fasta(.gz)     genome assembly FASTA
  main_outdir           # optional: main output directory (default: sdr_id_results)

Output:
  $depth_outdir/<genome>/
  └── bams/
      ├── <indiv_1>.sorted.bam
      └── <indiv_n>.sorted.bam
  └── mosdepth/
      ├── <indiv_1>.
      └── <indiv_n>.
  ├── 
  ├── 
  ├── 
  └── 

Requires:
 - bwa-mem2 v2.3 (https://github.com/bwa-mem2/bwa-mem2)
 - samtools v1.24 (https://www.htslib.org/)
 - mosdepth v0.3.14 (https://github.com/brentp/mosdepth)
 - R v6+ (https://www.r-project.org)"
}
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]; then
  show_help
  exit 0
elif [[ "$#" -lt 1 ]]; then
  echo "Error: at least 1 positional argument is required."
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
if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
  echo "Slurm array job ID: ${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID"
fi

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

# Input
genome="${1:-"assemblies/corteva_v2_ha.fa"}" # genome FASTA of one sex (e.g., female, U)
asm_id="${2:-"SL-CT1-FG-3"}" # genome assembly ID

# Output
genome_basename="$(basename "$genome")" # genome filename basename
genome_base="${genome_basename%%.fa*}" # genome filename base
depth_out="$depth_outdir/$genome_base" # job output subdirectory

# Parse individual ID and read paths
indiv_id="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$depth_out/$depth_tab_file" | cut -f1)"
mapfile -t reads < <(grep "${indiv_id}_" "$depth_out/$depth_reads_file")

# Validate that at least one complete R1/R2 pair was found
if (( ${#reads[@]} == 0 )); then
  echo "Error: No reads found for individual $indiv_id."
  exit 1
elif (( ${#reads[@]} % 2 != 0 )); then
  echo "Error: Found an odd number of read files for $indiv_id: ${#reads[@]}"
  printf '%s\n' "${reads[@]}"
  exit 1
fi

# Align samples to genome assembly
conda activate bwa-mem2
echo -n "bwa-mem2 version: "
bwa-mem2 version 2>&1 | tail -n +3 || exit 1
samtools --version | head -n2 || exit 1
echo
# Process each consecutive R1/R2 pair
bam_args=()
for ((i = 0; i < ${#reads[@]}; i += 2)); do
  R1="${reads[$i]}" # forward reads
  R2="${reads[$((i + 1))]}" # reverse reads
  R1_basename="$(basename "$R1")" # read filename basename
  sample_id="${R1_basename%%_R1*}" # read filename base (= sample ID)
  sample_bam="$TMPDIR/$sample_id.sorted.bam" # per sample library BAM
  # If multiple libraries, sort sample library BAM by query
  sample_bam_q="$TMPDIR/$sample_id.query_sorted.bam"
  
  # Define read groups for SAM headers from sample ID and FASTQ file read IDs
  # Input sample ID format = SampleID_Library_RunID_Barcode_Plate_Lane_Read
  # Create array 'rgs1' from sample ID by splitting on '_'
  IFS='_' read -r -a rgs1 <<< "$sample_id"
  # Read first line of FASTQ to get read ID
  read_id="$(zcat "$R1" | head -n1 | awk '{print $1}' | sed 's/@//g')"
  # Create array 'rgs2' from read ID by splitting on ':'
  IFS=':' read -r -a rgs2 <<< "$read_id"
  # Index rgs1 (sample ID array) and rgs2 (read ID array) to generate read group values
  # Read group ID = {FLOWCELL_BARCODE}.{LANE}
  rg_id="${rgs2[2]}.${rgs2[3]}"
  echo "Read group ID = ID:$rg_id"
  # Read group SaMple name (SM) = shortened SampleID within input sample ID (no UniqueID)
  rg_sm="${rgs1[0]}"
  echo "Read group sample name = SM:$rg_sm"
  # Read group Platform Unit (PU) = {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_NAME}
  rg_pu="${rgs2[2]}.${rgs2[3]}.${rgs1[0]}"
  echo "Read group platform unit = PU:$rg_pu"
  # Read group LiBrary identifier = SampleID-UniqueID
  rg_lb="${rgs1[0]}-${rgs1[1]}"
  echo "Read group library identifier = LB:$rg_lb"
  # Read group PLatform (PL) = ILLUMINA
  rg_pl="ILLUMINA"
  echo "Read group platform = PL:$rg_pl"
  
  # Allocate from remaining memory (/tmp stores in RAM)
  tmp_mem="$(du -m "$TMPDIR" | awk '{print $1}')" # unit=M
  mem="$(( SLURM_MEM_PER_NODE - tmp_mem ))" # unit=M
  # Memory per thread for samtools sort (75% of total to avoid OOM error)
  mem_per_thread="$(( 75 * mem / nthrd / 100 ))" # unit=M
  
  # Use bwa-mem2 to align PE reads for each sample library
  set -x
  bwa-mem2 mem -t "$nthrd" \
    -R "@RG\tID:${rg_id}\tSM:${rg_sm}\tLB:${rg_lb}\tPL:${rg_pl}\tPU:${rg_pu}" \
    "$depth_out/$genome_base" "$R1" "$R2" |
    samtools sort -@ "$nthrd" -m "${mem_per_thread}M" -O bam -o "$sample_bam"
  { set +x; } 2>/dev/null
  if [[ "${#reads[@]}" -gt 2 ]]; then
    set -x
    samtools sort -n -@ "$nthrd" -m "${mem_per_thread}M" -O bam \
      -o "$sample_bam_q" "$sample_bam"
    { set +x; } 2>/dev/null
    sample_bam="$sample_bam_q" # udpate sample BAM to query sorted
  fi

  # Save updated file path and add required "-I" flag before each filename
  check_output "$sample_bam"
  bam_args+=("-I" "$sample_bam")
done
conda deactivate

# GATK4 Spark parameters
# Allocate from remaining memory (/tmp stores in RAM)
tmp_mem="$(du -m "$TMPDIR" | awk '{print $1}')" # force unit=M
mem="$(( SLURM_MEM_PER_NODE - ( tmp_mem * 2 ) ))" # double for output file
# Save 4gb for driver/off-heap and other threads
driver_mem=4000 # unit=M
exec_mem="$(( mem - driver_mem ))" # unit=M
offheap_mem=2000 # off-heap shuffle buffer (unit=M)
# Partitons used for shuffling data during joins/aggregations
shuffle_parts="$(( nthrd * 3 ))" # avoid default 200

# Load GATK4
module purge
module load gcc/13.3.0 gatk/4.5.0.0 python/3.13.2 openjdk/21.0.0_35
gatk --version || exit 1

# Run GATK4's MarkDuplicates
indiv_bam="${indiv_id}.sorted.bam" # one BAM for biological individual
set -x
gatk MarkDuplicatesSpark \
  -R "$genome" \
  --tmp-dir "$TMPDIR" \
  "${bam_args[@]}" \
  -O "$TMPDIR/$indiv_bam" \
  -VS LENIENT \
  --spark-master "local[$nthrd]" \
  --conf spark.local.dir="$TMPDIR" \
  --conf spark.task.cpus=1 \
  --conf spark.executor.cores="$nthrd" \
  --conf spark.driver.cores=1 \
  --conf spark.executor.memory="${exec_mem}M" \
  --conf spark.driver.memory="${driver_mem}M" \
  --conf spark.memory.offHeap.enabled=true \
  --conf spark.memory.offHeap.size="${offheap_mem}M" \
  --conf spark.sql.shuffle.partitions="$shuffle_parts" \
  --conf spark.shuffle.file.buffer=32k \
  --conf spark.executor.heartbeatInterval="60s" \
  --conf spark.network.timeout="300s"
{ set +x; } 2>/dev/null
module purge
check_output "$TMPDIR/$indiv_bam"

# Copy results from temporary storage
set -x
mkdir -p "$depth_out/bams" # create BAM output directory
rsync --info=progress2 -ah "$TMPDIR/$indiv_bam" "$TMPDIR/$indiv_bam.bai" \
  "$depth_out/bams/"
{ set +x; } 2>/dev/null
check_output "$depth_out/bams/$indiv_bam"
check_output "$depth_out/bams/$indiv_bam.bai"

# Calculate read depth across genome in 10-kb windows
conda activate mosdepth
mkdir -p "$depth_out/mosdepth"
set -x
mosdepth \
  --threads "$nthrd" \
  --by 10000 \
  --fast-mode \
  --mapq 30 \
  --no-per-base \
  "$depth_out/mosdepth/$indiv_id" \
  "$depth_out/bams/$indiv_bam"
{ set +x; } 2>/dev/null
conda deactivate

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
