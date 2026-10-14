## Initialization
# Clear environment
rm(list = ls())
# Load required packages
suppressPackageStartupMessages(library(tidyverse))
library(geosphere)
library(ggpubr)


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
  target_cov <- 25
  # Genome assembly size (in bp)
  # asm_size <- 653974634 # corteva_v2_ha
  asm_size <- 712691060 # corteva_v2
  # Genome assembly ID
  asm_id <- "SL-CT1-FG-3"
  # Output directories
  ygs_out <- "sdr_id_results/ygs_sdr_results/corteva_v2_ha"
  depth_out <- "sdr_id_results/depth_sdr_results/corteva_v2_ha"
} else {
  line_args <- commandArgs(trailingOnly = T)
  target_cov <- line_args[1]
  asm_size <- line_args[2]
  asm_id <- line_args[3]
  ygs_out <- line_args[4]
  depth_out <- line_args[5]
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
ygs_os_tab_file <- "os_table.tsv"
ygs_ss_tab_file <- "ss_table.tsv"
ygs_os_reads_file <- "os.txt"
ygs_ss_reads_file <- "ss.txt"
depth_tab_file <- "os_ss_table.tsv"
depth_reads_file <- "os_ss.txt"
if (dir.exists(ygs_out)) {
  ygs_os_tab_file <- file.path(ygs_out, ygs_os_tab_file)
  ygs_os_reads_file <- file.path(ygs_out, ygs_os_reads_file)
  ygs_ss_tab_file <- file.path(ygs_out, ygs_ss_tab_file)
  ygs_ss_reads_file <- file.path(ygs_out, ygs_ss_reads_file)
}
if (dir.exists(depth_out)) {
  depth_tab_file <- file.path(depth_out, depth_tab_file)
  depth_reads_file <- file.path(depth_out, depth_reads_file)
}


## Data filtering
# Filter MultiQC general stats table for relevant columns
mqc_raw <- read_tsv(mqc_stats_file, show_col_types = F)
mqc_df <- mqc_raw %>%
  select(Sample,
         Mean_Read_Depth = `samtools_coverage-meandepth`,
         Coverage_Percent = `samtools_coverage-coverage`,
         Mapped_Percent = `samtools_stats-reads_mapped_percent`) %>%
  # Remove non-sample level rows
  filter(!is.na(Mapped_Percent))
# Add aligned bases column
aln_df <- read_tsv(aln_file, show_col_types = F) %>%
  rename(Aligned_Bases = `Aligned Bases`)
ids_df <- mqc_df %>%
  left_join(aln_df, by = join_by(Sample))
# Import population codes
pop_codes <- read_csv(pop_codes_file, show_col_types = F)
# Annotate full data frame with metadata
ids_df <- metaDf(ids_df, "Sample")

## Assembly metadata matching
# Keep only samples from assembly sample location
asm_df <- ids_df %>%
  filter(Sample == asm_id)
filt_ids_df <- ids_df %>%
  filter(Location == asm_df$Location) %>%
  # Sort by mean read depth
  arrange(desc(Mean_Read_Depth))


## YGS SDR ID
cat("YGS read filtering\n")
# YGS opposite sex reads
os_ids_df <- filt_ids_df %>%
  # Filter for samples of opposite sex as assembly
  filter(sex != asm_df$sex) %>%
  # Sort by mean read depth
  arrange(desc(Mean_Read_Depth))
# Calculate estimated genome coverage with reads from aligned bases
os_bases <- sum(os_ids_df$Aligned_Bases, na.rm = T)
os_cov <- os_bases/asm_size
# Filter for YGS
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
  rowwise() %>%
  mutate(Reads = paste(list.files(path = paste0(popgen_dir, "/trimmed_reads"),
                                  pattern = paste0(Sample, "_"), full.names = T),
                       collapse = ";")) %>%
  ungroup() %>%
  separate_longer_delim(Reads, delim = ";") %>%
  select(Reads)
cat("Total (expected) aligned bases:", round(os_bases*1e-9, digits = 1), "Gb\n")
cat("Number of samples:", length(unique(os_ids_df$Sample)), "\n")
cat("Number of sequencing runs:", dim(os_reads)[1]/2, "\n")

# YGS same sex reads
ss_ids_df <- filt_ids_df %>%
  # Filter for samples of same sex as assembly
  filter(sex == asm_df$sex) %>%
  # Sort by mean read depth
  arrange(desc(Mean_Read_Depth))
# Calculate estimated genome coverage with reads from aligned bases
ss_bases <- sum(ss_ids_df$Aligned_Bases, na.rm = T)
ss_cov <- ss_bases/asm_size
# Filter for YGS
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
  rowwise() %>%
  mutate(Reads = paste(list.files(path = paste0(popgen_dir, "/trimmed_reads"),
                                  pattern = paste0(Sample, "_"), full.names = T),
                       collapse = ";")) %>%
  ungroup() %>%
  separate_longer_delim(Reads, delim = ";") %>%
  select(Reads)
cat("Total (expected) aligned bases:", round(ss_bases*1e-9, digits = 1), "Gb\n")
cat("Number of samples:", length(unique(ss_ids_df$Sample)), "\n")
cat("Number of sequencing runs:", dim(ss_reads)[1]/2, "\n")

# Write YGS read files
write_tsv(os_ids_df, ygs_os_tab_file) # opposite sex metadata
write_tsv(ss_ids_df, ygs_ss_tab_file) # same sex metadata
write_tsv(os_reads, ygs_os_reads_file, col_names = F) # paths to opposite sex reads
write_tsv(ss_reads, ygs_ss_reads_file, col_names = F) # paths to same sex reads


## Sex-based depth comparison SDR ID
# Filter for 10 samples per sex with min. coverage & closest to assembly sample
sort_ids_df <- ids_df %>%
  rowwise() %>%
  # Use Haversine formula to compute distances (in km)
  mutate(Km_from_Asm = distHaversine(c(asm_df$Longitude, asm_df$Latitude),
                                       c(Longitude, Latitude)) / 1e3,
         Km_str = as.character(round(Km_from_Asm))) %>%
  ungroup() %>%
  arrange(Km_from_Asm) %>%
  # mutate(Location = factor(Location, levels = unique(Location))) %>%
  mutate(Km_str = factor(Km_str, levels = unique(Km_str)))

# Plots
p_distance_box <- ggplot(sort_ids_df,
                         aes(x = Km_str, y = Mean_Read_Depth, fill = sex)) +
  geom_boxplot() +
  labs(x = "Distance from assembly sample (km)", y = "Mean read depth")
p_sex_vio <- ggplot(sort_ids_df,
                         aes(x = sex, y = Mean_Read_Depth, color = sex)) +
  geom_violin(aes(fill = sex), color = "black", alpha = 0.8, trim = F) +
  geom_boxplot(color = "black", width = 0.2) +
  geom_jitter(color = "black", width = 0.05) +
  labs(y = "Mean read depth")

cutoffs <- sort(unique(sort_ids_df$Km_from_Asm))
depth_by_cutoff <- map_dfr(
  cutoffs,
  function(cutoff) {
    sort_ids_df %>%
      ungroup() %>%
      filter(Km_from_Asm <= cutoff) %>%
      mutate(cutoff_str = as.character(round(cutoff)),
             cutoff_str = factor(cutoff_str, levels = unique(cutoff_str))) %>%
      summarise(cutoff_km = cutoff, cutoff_str = unique(cutoff_str), n = n(),
                `Median Read Depth` = median(Mean_Read_Depth, na.rm = T),
                # `Total Read Depth` = sum(Mean_Read_Depth, na.rm = T),
                # `Median Aligned Bases (Gb)` = median(Aligned_Bases)/1e9,
                `Total Aligned Bases (Gb)` = sum(Aligned_Bases)/1e9,
                .by = sex) %>%
      pivot_longer(-c(sex, contains("cutoff"), n),
                   names_to = "stat", values_to = "value")
  }
)
p_cutoff <- ggplot(depth_by_cutoff, aes(x = cutoff_km, y = value, color = sex)) +
  geom_line() +
  geom_point(aes(fill = sex, size = n), color = "black", shape = 21,
             alpha = 0.5) +
  labs(x = "Maximum distance from assembly sample (km)", y = "",
       size = "Number of samples") +
  facet_wrap(~ stat, scales = "free", nrow = 1) +
  scale_x_continuous(n.breaks = length(cutoffs))

p1 <- ggarrange(p_sex_vio, p_distance_box, nrow = 1, widths = c(1, 2))

sort_ids_df %>%
  # Minimum mean depth cutoff of 2
  filter(Mean_Read_Depth >= 2) %>%
  # CHOOSE PROGRAMATICALLY
  filter(Km_from_Asm <= 65) %>%
  # Sort by least distance from ref. location & highest depth
  arrange(Km_from_Asm, desc(Mean_Read_Depth)) %>%
  slice_head(n = 10, by = sex) %>%
  summarize(total_bases_gb = sum(Aligned_Bases)/1e9, .by = sex)

# Write read files for sex-based depth comparison
write_tsv(os_ss_df, depth_tab_file, col_names = F) # table of sample metadata
write_tsv(os_reads, depth_reads_file, col_names = F) # paths to reads

