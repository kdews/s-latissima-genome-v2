#!/bin/bash
#SBATCH -J genome_analysis_pipeline
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
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || (( $# < 2 ))
then
  echo "\
bash $job_name.sh
sbatch $job_name.sh
"
  exit 0
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

# Variables
in_tsv='genomes_metadata.tsv'
scripts_dir='s-latissima-genome-v2'
# assembly_file='s_latissima_genome_list.txt'
pan_seqFile='s_latissima_align1.txt'
prog_seqFile='s_latissima_prog_align.txt'
ref_id='corteva_v2'
quer_id='roscoff_v2'
pangenome_dir='cactus_pangenome_01-28-26_6039802'

# Input data
# Download assembly FASTAs and sequencing reads and name by assembly labels
sbatch "$scripts_dir"/get_data.sbatch "$in_tsv"

# Evaluate basic assembly statistics
# Optional: download BUSCO lineage first
sbatch "$scripts_dir"/lineage_download.sbatch
# Run QUAST and BUSCO
sbatch "$scripts_dir"/eval_genomes.sh "$in_tsv" \
  "$scripts_dir"/busco.sbatch "$scripts_dir"/quast.sbatch
# Generate BUSCO summary files and visualize
sbatch "$scripts_dir"/eval_genomes.sh "$in_tsv" \
  "$scripts_dir"/busco_compare.sbatch
# # Visualization (NEED TO ADD)
# "$scripts_dir"/scaffold_eval.R

# Cactus alignments
# Install Cactus
sbatch "$scripts_dir"/cactus_install_conda.sh

# 1. Cactus Pangenome
# Run Minigraph-Cactus alignment
sbatch "$scripts_dir"/cactus_pangenome.sbatch "$pan_seqFile" "$ref_id"
# Extract HAL from pangenome PAFs
sbatch "$scripts_dir"/hal_extract.sbatch "$pangenome_dir" "$ref_id" "$quer_id"
# Generate pangenome plots with R
sbatch "$scripts_dir"/pangenome_plots.sbatch
# Generate per-chromosome ODGI graphs
sbatch "$scripts_dir"/odgi_viz.sbatch "$pan_seqFile" "$pangenome_dir"
# Generate ntSynt-viz pangenome graph
sbatch "$scripts_dir"/ntSynt_viz.sbatch "$pan_seqFile" "$pangenome_dir"

# 2. Progressive Cactus
# Make Cactus guide tree
sbatch "$scripts_dir"/make_guide_tree.sbatch "$prog_seqFile"
# Filter assemblies for Cactus (minimum contig length = 1000000)
sbatch "$scripts_dir"/filter_assemblies.sbatch "$prog_seqFile" 1000000
# Run Progressive Cactus alignment
sbatch "$scripts_dir"/prog_cactus.sbatch "$prog_seqFile"
# Visualize resulting HAL (SCRIPT IN PROGRESS)
sbatch "$scripts_dir"/hal_visualize.sbatch

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
