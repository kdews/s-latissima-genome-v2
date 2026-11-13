#!/bin/bash
#SBATCH -J sdr_id
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p htcondor
#SBATCH --time=01:00:00
#SBATCH --mem=40g
#SBATCH --cpus-per-task=16

# Script variables
# Extract job name from script
if [[ -n ${SLURM_JOB_NAME:-} ]]; then
  job_name="$SLURM_JOB_NAME"
elif [[ -n ${BASH_SOURCE[0]:-} && ${BASH_SOURCE[0]} != "bash" ]]; then
  job_name="$(basename "${BASH_SOURCE[0]}")"
  job_name="${job_name%%.*}"
else
  job_name="sdr_id"
fi
main_outdir="${job_name}_results" # main SDR output directory
prep_job="sdr_prepare.py" # Python script to parse metadata for analysis
# YGS SDR job variables
ygs_job_name="ygs_sdr"
ygs_job="$ygs_job_name.sbatch" # Slurm job submission script
ygs_outdir="${ygs_job_name}_results" # output subdirectory
# URL to YGS.pl from Carvalho and Clark (2013)
ygs_url="https://genome.cshlp.org/content/genome/suppl/2013/09/24/gr.156034.113.DC1/YGS.txt"
ygs_pl="YGS.pl" # local name of YGS script
ygs_plot="ygs_plot.sbatch" # job script to plot YGS results
depth_plot="depth_plot.sbatch" # job script to plot depths results
# Read depth SDR job variables
depth_job_name="depth_sdr"
depth_job="$depth_job_name.sbatch" # Slurm job submission script
depth_outdir="${depth_job_name}_results" # output subdirectory
scripts_dir="s-latissima-genome-v2" # scripts directory
plan_tsv="sdr_plan.tsv" # SDR analysis TSV parsed from input TSV

# Help message
show_help () {
  echo "\
Runs both k-mer-based SDR identification (Carvalho and Clark 2013) and read
depth-based SDR identification from Barrera-Redondo et al. (2025). Requires
reads from opposite sex (os) with respect to the assembly, with a set of same
sex (ss) reads used to validate assembly k-mers. Assembly metadata from
genomes_metadata.tsv is used to source genomes and reads from
read_lists/<label>_{female,male}.list. Ploidy is used to configure the YGS 
analysis as follows:
  - U: target_U YGS
  - V: target_V YGS
  - UV: both
Assemblies without designated female/male read lists are skipped.
--dry-run validates inputs, writes manifests, and prints jobs without submitting.

Usage: sbatch $job_name.sbatch [genomes_metadata.tsv] [--dry-run]


Input:
  genomes_metadata.tsv  20-column assembly metadata TSV generated from the
                        per-assembly YAML files. Defaults to
                        s-latissima-genome-v2/genomes_metadata.tsv when that
                        path exists, otherwise genomes_metadata.tsv.

Output:
  $main_outdir/
  ├── $plan_tsv
  ├── $ygs_outdir/
  │   └── <label>/target_<U|V>/
  │       ├── <windowsize>_windows.contig_result                  # log for YGS contig mode
  │       ├── <windowsize>_windows.ctg_rep.gz                     # YGS bit-array of genomic repeat <k>-mers in > 1 contig
  │       ├── <windowsize>_windows.fasta.gz                       # genome FASTA split into <windowsize> sliding windows
  │       ├── <windowsize>_windows.gen_rep.gz                     # YGS bit-array of genomic repeat <k>-mers
  │       ├── <windowsize>_windows_os_k<k>_ss_k<k>.final_result   # speaks for itself
  │       ├── os.txt                                              # list of read file paths (os)
  │       ├── os_k<k>.fasta.gz                                    # <k>-mers extracted into FASTA (os)
  │       ├── os_k<k>.jf                                          # <k>-mer Jellyfish binary hash (os)
  │       ├── os_k<k>.trace.gz                                    # YGS <k>-mer bit-array (ss)
  │       ├── os_k<k>.trace_result                                # log for YGS trace mode (os)
  │       ├── ss.txt                                              # list of same sex read file paths (ss)
  │       ├── ss_k<k>.fasta.gz                                    # <k>-mers extracted into FASTA (ss)
  │       ├── ss_k<k>.jf                                          # <k>-mer Jellyfish binary hash (ss)
  │       ├── ss_k<k>.trace.gz                                    # YGS <k>-mer bit-array (ss)
  │       └── ss_k<k>.trace_result                                # log for YGS trace mode (ss)
  ├── $depth_outdir/
  │   └── <label>/
  │       ├── samples.tsv                                         # 'Sample Sex'
  │       ├── read_pairs.tsv                                      # 'Sample  Sex  Library  Pair_ID  R1  R2'
  │       ├── assembly_metadata.tsv                               # 'Label	Pretty_Label	Genome_Source	Genome_Method	Genome_Locator	Genome_Selector	Assembly_Reads_Source	Assembly_Reads_Method	Assembly_Reads_Platform	Assembly_Reads_Locator	Analysis_Reads_Platform	Analysis_Reads_Method	Analysis_Reads_Locator	Analysis_Reads_Female_Locator	Analysis_Reads_Male_Locator	Annotation_GFF	Sample_ID	Region	Sex	Ploidy'
  │       ├── <label>.0123                                        # bwa-mem2 index files
  │       ├── <label>.amb                                         # bwa-mem2 index files
  │       ├── <label>.ann                                         # bwa-mem2 index files
  │       ├── <label>.bwt.2bit.64                                 # bwa-mem2 index files
  │       ├── bams/                                               # individual-level BAMs
  │       │   ├── <indiv_n>.sorted.bam
  │       │   ├── <indiv_n>.sorted.bam.bai
  │       │   └── ...
  │       ├── mosdepth/                                           # mosdepth result files
  │       │   ├── <indiv_n>.mosdepth.global.dist.txt
  │       │   ├── <indiv_n>.mosdepth.region.dist.txt
  │       │   ├── <indiv_n>.mosdepth.summary.txt
  │       │   ├── <indiv_n>.regions.bed.gz
  │       │   ├── <indiv_n>.regions.bed.gz.csi
  │       │   └── ...
  │       └── depth_filt.png                                      # plots of sample read depths
  └── $scripts_dir/
      └── figures/
          ├── depth_all_<label>.png
          ├── depth_candidate_<label>.png
          ├── depth_diff_all_<label>.png
          ├── depth_log2_all_<label>.png
          ├── depth_log2_candidate_<label>.png
          ├── ygs_all_<label>.png
          └── ygs_candidate_<label>.png

Requires:
 - Python 3 (https://www.python.org)
 - Jellyfish v2.3.1+ (https://github.com/gmarcais/jellyfish)
 - pigz v2.8+ (https://zlib.net/pigz/)
 - Perl v5.23.1 (https://www.perl.org/)
 - bwa-mem2 v2.2.1 (https://bio-bwa.sourceforge.net/)
 - samtools v1.24 (https://www.htslib.org/)
 - mosdepth v0.3.14 (https://github.com/brentp/mosdepth)"
}
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]; then
  show_help
  exit 0
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

# Input
dry_run=false # set up optional dry run
in_tsv=""
for arg in "$@"; do
  if [[ "$arg" == --dry-run ]]; then
    dry_run=true
  elif [[ -z "$in_tsv" ]]; then
    in_tsv="$arg"
  else
    echo "Error: unexpected argument: $arg" >&2
    exit 1
  fi
done
if [[ -z "$in_tsv" ]]; then
  if [[ -f "$scripts_dir/genomes_metadata.tsv" ]]; then
    in_tsv="$scripts_dir/genomes_metadata.tsv"
  else
    in_tsv="genomes_metadata.tsv"
  fi
fi

# Download YGS.pl to scripts_dir (if needed), then resolve its absolute path
if [[ ! -s "$scripts_dir/$ygs_pl" ]]; then
  wget -O "$scripts_dir/$ygs_pl" "$ygs_url" || exit 1
fi
ygs_pl="$(realpath -e "$scripts_dir/$ygs_pl")" || exit 1


# Verify input and script dependencies
[[ -r "$in_tsv" ]] || { echo "Cannot read metadata: $in_tsv" >&2; exit 1; }
for script in "$prep_job" "$ygs_job" "$depth_job" "$ygs_plot" "$depth_plot"; do
  if [[ ! -r "$scripts_dir/$script" ]]; then
    echo "Missing script dependency: $scripts_dir/$script" >&2
    exit 1
  fi
done

# Parse metadata for job submission
mkdir -p "$main_outdir"
set -x
if python3 "$scripts_dir/$prep_job" "$in_tsv" "$main_outdir" "$ygs_outdir" \
  "$depth_outdir" "$plan_tsv"; then
  { set +x; } 2>/dev/null
  if [[ ! -s "$main_outdir/$plan_tsv" ]]; then
    echo 'No configured analyses.'
    exit 0
  fi
else
  { set +x; } 2>/dev/null
  echo 'Error: SDR preparation failed.' >&2
  exit 1
fi

# Iterate through $plan_tsv
while IFS=$'\t' read -r label genome sex n_samples; do
  depth_out="$main_outdir/$depth_outdir/$label"
  if [[ "$dry_run" == true ]]; then
    echo "PLAN depth: $label --array=1-$n_samples; depth_plot afterok"
  else
    index_missing=false
    for suffix in bwt.2bit.64 0123 amb ann pac; do
      [[ -s "$depth_out/$label.$suffix" ]] || index_missing=true
    done
    if [[ "$index_missing" == true ]]; then
      conda activate bwa-mem2
      set -x
      bwa-mem2 index -p "$depth_out/$label" "$genome" || exit 1
      { set +x; } 2>/dev/null
      conda deactivate
    fi
    conda activate samtools
    set -x
    samtools faidx "$genome" || exit 1
    { set +x; } 2>/dev/null
    conda deactivate
    # Load GATK4
    module load gcc/13.3.0 gatk/4.5.0.0 python/3.13.2 openjdk/21.0.0_35
    gatk --version || exit 1
    # Build GATK4 sequence dictionary from genome FASTA
    if [[ ! -s "${genome%.fa}.dict" ]]; then
      set -x
      gatk CreateSequenceDictionary -R "$genome" -O "${genome%.fa}.dict" \
        || exit 1
      { set +x; } 2>/dev/null
    fi
    module purge
    # set -x
    # d_id="$(sbatch --parsable --array="1-$n_samples" "$scripts_dir/$depth_job" \
    #   "$genome")" \
    #   || exit 1
    # { set +x; } 2>/dev/null
    # set -x
    # sbatch --dependency="afterok:$d_id" "$scripts_dir/$depth_plot" "$depth_out" \
    #   "$scripts_dir"
    # { set +x; } 2>/dev/null
  fi
  targets=("$sex")
  [[ "$sex" != UV ]] || targets=(U V)
  for target in "${targets[@]}"; do
    ygs_out="$main_outdir/$ygs_outdir/$label/target_$target"
    if [[ "$dry_run" == true ]]; then
      echo "PLAN YGS: $label target_$target; ygs_plot afterok"
    # else
    #   set -x
    #   y_id="$(sbatch --parsable "$scripts_dir/$ygs_job" "$genome" "$target" \
    #     "$scripts_dir")" \
    #     || exit 1
    #   { set +x; } 2>/dev/null
    #   set -x
    #   sbatch --dependency="afterok:$y_id" "$scripts_dir/$ygs_plot" "$ygs_out" \
    #     "$scripts_dir"
    #   { set +x; } 2>/dev/null
    fi
  done
done < "$main_outdir/$plan_tsv"

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
