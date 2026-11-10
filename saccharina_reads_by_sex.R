## Initialization
# Clear environment
rm(list = ls())
# Load required packages
suppressPackageStartupMessages(library(tidyverse))

## Functions
# Match assembly metadata in sample set and select one sex
matchAsm <- function(df, asm_id, sx) {
  asm_df <- df %>%
    filter(ID == asm_id)

  if (nrow(asm_df) != 1) {
    stop("Expected exactly one match for assembly sample ", asm_id,
         "; found ", nrow(asm_df), ".", call. = F)
  }

  # Keep only samples that match the assembly location and subspecies
  df %>%
    filter(
      Location == asm_df$Location[[1]],
      subspecies == asm_df$subspecies[[1]],
      sex == sx
    )
}

# Add paths for all read files belonging to each selected individual
addReads <- function(df, wgs_files) {
  reads_df <- df %>%
    mutate(
      Reads = map(ID, ~ wgs_files[
        startsWith(basename(wgs_files), paste0(.x, "_"))
      ])
    )

  missing_ids <- reads_df %>%
    filter(lengths(Reads) == 0) %>%
    pull(ID)
  if (length(missing_ids) > 0) {
    warning("No WGS files found for: ", paste(missing_ids, collapse = ", "))
  }

  reads_df %>%
    filter(lengths(Reads) > 0) %>%
    unnest_longer(Reads)
}

# Write one read path per line without a header
writeReads <- function(df, label, sx, out_dir) {
  out_file <- file.path(out_dir, paste0(label, "_", sx, ".list"))
  reads <- df %>%
    pull(Reads) %>%
    unique() %>%
    sort()

  if (length(reads) == 0) {
    stop("No ", sx, " reads found for ", label, ".", call. = F)
  }
  missing_reads <- reads[!file.exists(reads)]
  if (length(missing_reads) > 0) {
    stop("Read files not found:\n", paste(missing_reads, collapse = "\n"),
         call. = F)
  }

  writeLines(reads, out_file)
  message("Wrote ", length(reads), " ", sx, " read paths to ", out_file)
}

## Input
# Only take command-line input if not running interactively
if (interactive()) {
  # Set working directory
  wd <- "/scratch2/kdeweese/corteva_genome"
  setwd(wd)
  label <- "lindell_v0"
  asm_id <- "SL-CT1-FG-3"
  scripts_dir <- "s-latissima-genome-v2"
  out_dir <- "read_lists"
} else {
  line_args <- commandArgs(trailingOnly = T)
  if (length(line_args) != 4) {
    stop("Usage: Rscript saccharina_reads_by_sex.R ",
         "<label> <assembly_sample_id> <scripts_dir> <output_dir>",
         call. = F)
  }
  label <- line_args[1]
  asm_id <- line_args[2]
  scripts_dir <- line_args[3]
  out_dir <- line_args[4]
}

pop_codes_file <- "Saccharina_pop_codes_WHOI.csv"
metaDf_file <- "metaDf.R"
# Source project files from the scripts directory when it exists
if (dir.exists(scripts_dir)) {
  pop_codes_file <- file.path(scripts_dir, pop_codes_file)
  metaDf_file <- file.path(scripts_dir, metaDf_file)
}

# Directory containing renamed Saccharina WGS reads
popgen_dir <- "/project2/noujdine_61/kdeweese/latissima/popgen_with_feeling"
ids_file <- file.path(popgen_dir, "indexer/individuals_file.txt")
wgs_dir <- file.path(popgen_dir, "wgs")

# Check inputs exist
if (!file.exists(ids_file)) stop("File ", ids_file, " not found.")
if (!dir.exists(wgs_dir)) stop("Directory ", wgs_dir, " not found.")
if (!file.exists(pop_codes_file)) stop("File ", pop_codes_file, " not found.")
if (!file.exists(metaDf_file)) stop("File ", metaDf_file, " not found.")
dir.create(out_dir, recursive = T, showWarnings = F)

## Annotate metadata
source(metaDf_file)
pop_codes <- read_csv(pop_codes_file, show_col_types = F)
ids_df <- read_table(ids_file, col_names = c("ID"), show_col_types = F)
ids_df <- metaDf(ids_df, "ID", pop_codes)

# Select samples matching the assembly location and subspecies
female_df <- matchAsm(ids_df, asm_id, "female")
male_df <- matchAsm(ids_df, asm_id, "male")

# Add existing local FASTQ paths only after selecting candidate samples
wgs_files <- list.files(wgs_dir, full.names = T)
female_reads_df <- addReads(female_df, wgs_files)
male_reads_df <- addReads(male_df, wgs_files)

## Write output
writeReads(female_reads_df, label, "female", out_dir)
writeReads(male_reads_df, label, "male", out_dir)
