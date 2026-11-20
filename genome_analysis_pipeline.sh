#!/bin/bash
#SBATCH -J genome_analysis_pipeline
#SBATCH -o %x_logs/%x_%j.log
#SBATCH --time=05:00:00
#SBATCH --mem=3gb

# Script variables
# Extract job name from script
if [[ -n ${SLURM_JOB_NAME:-} ]]; then
  job_name="$SLURM_JOB_NAME"
elif [[ -n ${BASH_SOURCE[0]:-} && ${BASH_SOURCE[0]} != "bash" ]]; then
  job_name="$(basename "${BASH_SOURCE[0]}")"
  job_name="${job_name%%.*}"
else
  job_name="genome_analysis_pipeline"
fi

# Help message
show_help () {
  echo "\
Usage:
  bash $job_name.sh
  sbatch $job_name.sh"
}
# Catch help flag
[[ $1 == "-h" || $1 == "--help" ]] && { show_help; exit 0; }

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

# Prints (to sdterr) and submits Slurm batch job
submit_job() {
  printf '%q ' sbatch "$@" >&2
  printf '\n' >&2
  sbatch "$@"
}

# Variables
in_tsv='genomes_metadata.tsv'
scripts_dir='s-latissima-genome-v2'
# assembly_file='s_latissima_genome_list.txt'
pan_seqFile='s_latissima_align1.txt'
prog_seqFile='s_latissima_prog_align.txt'
ref_id='corteva_v2'
quer_id='roscoff_v2'
pangenome_dir='cactus_pangenome_results/01-28-26_6039802'

# Input data
# Generate genome metadata TSV from YAMLs in <scripts_dir>
submit_job "$scripts_dir/yaml_to_metadata.sbatch"
# Download assembly FASTAs and sequencing reads and name by assembly labels
submit_job "$scripts_dir/get_data.sbatch" "$scripts_dir/$in_tsv"

# Evaluate basic assembly statistics
# Optional: download BUSCO lineage first
submit_job "$scripts_dir/lineage_download.sbatch"
# Run QUAST and BUSCO
submit_job "$scripts_dir/eval_genomes.sh" "$scripts_dir/$in_tsv" \
  "$scripts_dir/busco.sbatch" "$scripts_dir/quast.sbatch"
# Generate BUSCO summary files and visualize
submit_job "$scripts_dir/eval_genomes.sh" "$scripts_dir/$in_tsv" \
  "$scripts_dir/busco_compare.sbatch"
# # Visualization (NEED TO ADD)
# "$scripts_dir"/scaffold_eval.R

# Cactus alignments
# Install Cactus
submit_job "$scripts_dir/cactus_install_conda.sh"

# 1. Cactus Pangenome
# Run Minigraph-Cactus alignment
submit_job "$scripts_dir/cactus_pangenome.sbatch" "$scripts_dir/$pan_seqFile" "$ref_id"
# Extract HAL from pangenome PAFs
submit_job "$scripts_dir/hal_extract.sbatch" "$pangenome_dir" "$ref_id" "$quer_id"
# Generate pangenome plots with R
submit_job "$scripts_dir/pangenome_plots.sbatch"
# Generate per-chromosome ODGI graphs
submit_job "$scripts_dir/odgi_viz.sbatch" "$scripts_dir/$pan_seqFile" "$pangenome_dir"
# Generate ntSynt-viz pangenome graph
submit_job "$scripts_dir/ntSynt_viz.sbatch" "$scripts_dir/$pan_seqFile" "$pangenome_dir"

# 2. Progressive Cactus
# Make Cactus guide tree
submit_job "$scripts_dir/make_guide_tree.sbatch" "$scripts_dir/$prog_seqFile"
# Filter assemblies for Cactus (minimum contig length = 1000000)
submit_job "$scripts_dir/filter_assemblies.sbatch" "$scripts_dir/$prog_seqFile" 1000000
# Run Progressive Cactus alignment
submit_job "$scripts_dir/prog_cactus.sbatch" "$scripts_dir/$prog_seqFile"
# Visualize resulting HAL (SCRIPT IN PROGRESS)
submit_job "$scripts_dir/hal_visualize.sbatch"

# 3. Identify genome SDR
submit_job "$scripts_dir/sdr_id.sbatch" "$scripts_dir/$in_tsv"

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
