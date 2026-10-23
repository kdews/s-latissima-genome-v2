## Initialization
# Clear environment
rm(list = ls())
# Load required packages
suppressPackageStartupMessages(library(tidyverse))

# DELETE?
# Set working directory
wd <- "/scratch2/kdeweese/corteva_genome"
setwd(wd)

# Scripts directory
scripts_dir <- "s-latissima-genome-v2"
# Directory containing renamed S. latissima WGS reads
popgen_dir <- "/project2/noujdine_61/kdeweese/latissima/popgen_with_feeling"
ids_file <- file.path(popgen_dir, "indexer/individuals_file.txt")
wgs_dir <- file.path(popgen_dir, "wgs")

# Check inputs exist
if (!file.exists(ids_file)) stop("File ", ids_file, "not found.")
if (!dir.exists(wgs_dir)) stop("Directory ", wgs_dir, " not found.")

# Metadata annotation
pop_codes_file <- "Saccharina_pop_codes_WHOI.csv" # S. latissima population info
metaDf_file <- "metaDf.R" # R script that contains 'metaDf' 
# Source from scripts directory, if exists
if (dir.exists(scripts_dir)) {
  pop_codes_file <- file.path(scripts_dir, pop_codes_file)
  metaDf_file <- file.path(scripts_dir, metaDf_file)
}
source(metaDf_file)
# Import population codes
pop_codes <- read_csv(pop_codes_file, show_col_types = F)

# Import
ids <- read_table(ids_file, col_names = c("ID"), show_col_types = F)

metaDf(ids, "ID", pop_codes)

