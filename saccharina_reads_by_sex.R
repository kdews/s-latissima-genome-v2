## Initialization
# Clear environment
rm(list = ls())
# Load required packages
suppressPackageStartupMessages(library(tidyverse))

## Functions
# Match assembly metadata in sample set & select sex
matchAsm <- function(df, asm_id, sx) {
  asm_df <- df %>%
    filter(ID == asm_id)
  # Keep only samples that match assembly sample location & subspecies
  df %>%
    filter(Location == asm_df$Location &
             subspecies == asm_df$subspecies &
             # Select female or male
             sex == sx)
}
# Write each set of reads
writeReads <- function(df, asm_id, sx) {
  
}

## Input
# Only take command line input if not running interactively
if (interactive()) {
  # Set working directory
  wd <- "/scratch2/kdeweese/corteva_genome"
  setwd(wd)
  sl_asm_id <- "SL-CT1-FG-3"
  sa_asm_id <- "SA-CB-7-MG-3"
  # Scripts directory
  scripts_dir <- "s-latissima-genome-v2"
} else {
  line_args <- commandArgs(trailingOnly = T)
  sl_asm_id <- line_args[1]
  sa_asm_id <- line_args[2]
  scripts_dir <- line_args[3]
}
pop_codes_file <- "Saccharina_pop_codes_WHOI.csv" # S. latissima population info
metaDf_file <- "metaDf.R" # R script that contains 'metaDf' 
# Source from scripts directory, if exists
if (dir.exists(scripts_dir)) {
  pop_codes_file <- file.path(scripts_dir, pop_codes_file)
  metaDf_file <- file.path(scripts_dir, metaDf_file)
}
# Directory containing renamed S. latissima WGS reads
popgen_dir <- "/project2/noujdine_61/kdeweese/latissima/popgen_with_feeling"
ids_file <- file.path(popgen_dir, "indexer/individuals_file.txt")
wgs_dir <- file.path(popgen_dir, "wgs")

# Check inputs exist
if (!file.exists(ids_file)) stop("File ", ids_file, "not found.")
if (!dir.exists(wgs_dir)) stop("Directory ", wgs_dir, " not found.")

## Annotate metadata
source(metaDf_file)
# Import population codes
pop_codes <- read_csv(pop_codes_file, show_col_types = F)
# Import individual IDs
ids_df <- read_table(ids_file, col_names = c("ID"), show_col_types = F)
# Add metadata
ids_df <- metaDf(ids_df, "ID", pop_codes)
# Add read paths
ids_df <- ids_df %>%
  rowwise() %>%
  mutate(Reads = paste(list.files(path = paste0(popgen_dir, "/wgs"),
                                  pattern = paste0(ID, "_"), full.names = T),
                       collapse = ";")) %>%
  ungroup() %>%
  separate_longer_delim(Reads, delim = ";")
# Match assembly IDs to sample data and separate female and male
sl_f_df <- matchAsm(ids_df, sl_asm_id, "female")
sl_m_df <- matchAsm(ids_df, sl_asm_id, "male")
sa_f_df <- matchAsm(ids_df, sa_asm_id, "female")
sa_m_df <- matchAsm(ids_df, sa_asm_id, "male")
# Separate out S. latissima from S. angustissima and female from male
sl_f_reads <- sl_f_df %>% select(Reads)
sl_m_reads <- sl_m_df %>% select(Reads)
sa_f_reads <- sa_f_df %>% select(Reads)
sa_m_reads <- sa_m_df %>% select(Reads)


## Write output 
write_tsv(sl_f_reads, "reads/")

