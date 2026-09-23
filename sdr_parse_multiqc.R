## Initialization
# Clear environment
rm(list = ls())
# Load required packages
suppressPackageStartupMessages(library(tidyverse))

## Functions
# Annotate dataframe with metadata, where column "label" contains individual IDs
metaDf <- function(df, label) {
  sex_convert <- c("FG" = "female", "MG" = "male")
  subsp_convert <- c("SL" = "latissima", "SA" = "angustissima")
  # Extract metadata from sample IDs
  df_annot <- df %>%
    mutate(
      rename_id = gsub("LIS-F1-3", "SL-SNE-1-FG-3", .data[[label]]),
      rename_id = gsub("CT1", "CT1-0", rename_id),
      rename_id = gsub("Female", "FG-0", rename_id)
    ) %>%
    separate_wider_delim(
      rename_id, delim = "-",
      names = c("subspecies", "pop", "n_samp", "sex", "n_gam")
    ) %>%
    mutate(
      subspecies = subsp_convert[subspecies],
      sex = sex_convert[sex]
    ) %>%
    # Merge with table of population code metadata
    left_join(pop_codes, by = join_by(pop == Abbreviation))
  return(df_annot)
}

## Input
# Only take command line input if not running interactively
if (interactive()) {
  # Set working directory
  wd <- "/scratch2/kdeweese/corteva_genome"
  setwd(wd)
  # Target coverage of genome with reads
  target_cov <- 60
  # Genome assembly size (in bp)
  asm_size <- 712691060
  # Genome assembly ID
  asm_id <- "SL-CT1-FG-3"
  # Output directory
  outdir <- "sdr_id_results"
} else {
  line_args <- commandArgs(trailingOnly = T)
  target_cov <- line_args[1]
  asm_size <- line_args[2]
  asm_id <- line_args[3]
  outdir <- line_args[4]
}
target_cov <- as.numeric(target_cov)
asm_size <- as.numeric(asm_size)
# Directory containing popgen results
popgen_dir <- "/project2/noujdine_61/kdeweese/latissima/popgen_with_feeling"
# Find most recent MultiQC result directory name
mqc_dir <- list.files(path = popgen_dir, pattern = "multiqc_data",
                      full.names = T) %>%
  file.info() %>%
  arrange(desc(mtime)) %>%
  slice(1) %>%
  rownames()
# MultiQC general stats TSV
mqc_stats_file <- list.files(path = mqc_dir,
                             pattern = "multiqc_general_stats.txt",
                             recursive = T, full.names = T)
# Picard alignment summary, parsed into TSV by MultiQC
aln_file <- list.files(path = mqc_dir,
                       pattern = "picard_alignment_summary_Aligned_Bases.txt",
                       recursive = T, full.names = T)
# Directory containing VCF analysis
vcf_dir <- "/project2/noujdine_61/kdeweese/latissima/compare_vcfs/compare-s-lat-vcfs"
# S. latissima population codes CSV
pop_codes_file <- list.files(path = vcf_dir, recursive = T,
                             pattern = "Saccharina_pop_codes_WHOI.csv",
                             full.names = T)

## Output
os_ids_file <- "os_ids_table.tsv"
os_reads_file <- "os_reads.txt"
ss_ids_file <- "ss_ids_table.tsv"
ss_reads_file <- "ss_reads.txt"
if (dir.exists(outdir)) {
  os_ids_file <- file.path(outdir, os_ids_file)
  os_reads_file <- file.path(outdir, os_reads_file)
  ss_ids_file <- file.path(outdir, ss_ids_file)
  ss_reads_file <- file.path(outdir, ss_reads_file)
}

## Data filtering
# Filter MultiQC general stats table for relevant columns
mqc_df <- read_tsv(mqc_stats_file, show_col_types = F)
mqc_df <- mqc_df %>%
  select(Sample,
         Mean_Read_Depth = `samtools_coverage-meandepth`,
         Coverage_Percent = `samtools_coverage-coverage`,
         Mapped_Percent = `samtools_stats-reads_mapped_percent`) %>%
  # Remove non-sample level rows
  filter(!is.na(Mapped_Percent)) %>%
  # Sort by mean read depth
  arrange(desc(Mean_Read_Depth))
# Add aligned bases column
aln_df <- read_tsv(aln_file, show_col_types = F) %>%
  rename(Aligned_Bases = `Aligned Bases`)
ids_df <- mqc_df %>%
  left_join(aln_df, by = join_by(Sample))
# Import population codes
pop_codes <- read_csv(pop_codes_file, show_col_types = F)
# Annotate full data frame with metadata
ids_df <- metaDf(ids_df, "Sample")

# Keep only samples from assembly sample location
asm_df <- ids_df %>%
  filter(Sample == asm_id) %>%
  select(sex, Location)
filt_ids_df <- ids_df %>%
  filter(Location == asm_df$Location)

# Filter for samples of opposite sex as assembly
os_ids_df <- filt_ids_df %>%
  filter(sex != asm_df$sex)
# Calculate estimated genome coverage with reads from aligned bases
os_bases <- sum(os_ids_df$Aligned_Bases, na.rm = T)
os_cov <- os_bases/asm_size
while (os_cov > target_cov && nrow(os_ids_df) > 0) {
  os_ids_df <- os_ids_df %>%
    slice(-n())
  os_bases <- sum(os_ids_df$Aligned_Bases, na.rm = T)
  os_cov <- os_bases/asm_size
}
# Log filtering
cat("Log for opposite sex reads...\n")
cat("Estimated genome coverage: ", round(os_cov, digits = 1),
    "x (target=", target_cov, "x).\n", sep = "")
# Create data frame of sample IDs and read paths
os_reads <- os_ids_df %>%
  select(Sample) %>%
  mutate(Reads = paste(list.files(path = paste0(popgen_dir, "/trimmed_reads"),
                                  pattern = paste0(Sample, "_"), full.names = T),
                       collapse = ";")) %>%
  separate_longer_delim(Reads, delim = ";") %>%
  select(Reads)
cat("Total (expected) aligned bases:", round(os_bases*1e-9, digits = 1), "Gb\n")

# Filter for samples of same sex as assembly
ss_ids_df <- filt_ids_df %>%
  filter(sex == asm_df$sex)
# Calculate estimated genome coverage with reads from aligned bases
ss_bases <- sum(ss_ids_df$Aligned_Bases, na.rm = T)
ss_cov <- ss_bases/asm_size
while (ss_cov > target_cov && nrow(ss_ids_df) > 0) {
  ss_ids_df <- ss_ids_df %>%
    slice(-n())
  ss_bases <- sum(ss_ids_df$Aligned_Bases, na.rm = T)
  ss_cov <- ss_bases/asm_size
}
# Log filtering
cat("Log for same sex reads...\n")
cat("Estimated genome coverage: ", round(ss_cov, digits = 1),
    "x (target=", target_cov, "x).\n", sep = "")
# Create data frame of sample IDs and read paths
ss_reads <- ss_ids_df %>%
  select(Sample) %>%
  mutate(Reads = paste(list.files(path = paste0(popgen_dir, "/trimmed_reads"),
                                  pattern = paste0(Sample, "_"), full.names = T),
                       collapse = ";")) %>%
  separate_longer_delim(Reads, delim = ";") %>%
  select(Reads)
cat("Total (expected) aligned bases:", round(ss_bases*1e-9, digits = 1), "Gb\n")

## Write output files
# ID metadata tables
write_tsv(os_ids_df, os_ids_file) # opposite sex
write_tsv(ss_ids_df, ss_ids_file) # same sex
# Paths to reads
write_tsv(os_reads, os_reads_file, col_names = F) # opposite sex
write_tsv(ss_reads, ss_reads_file, col_names = F) # same sex
