## Initialization
# Clear environment
rm(list = ls())
# Load required packages
suppressPackageStartupMessages(library(tidyverse))

## Script variables
depth_cols <- c("Chr", "Start", "End", "Mean_Depth")

## Input
# Only take command line input if not running interactively
if (interactive()) {
  # Set working directory
  wd <- "/scratch2/kdeweese/corteva_genome"
  setwd(wd)
  # Depth results directory
  depth_out <- "sdr_id_results/depth_sdr_results/corteva_v2_ha"
  # Output directory
  outdir <- "s-latissima-genome-v2"
} else {
  line_args <- commandArgs(trailingOnly = T)
  depth_out <- line_args[1]
  outdir <- line_args[2]
}
mos_dir <- file.path(depth_out, "mosdepth") # mosdepth results dir
# Locate regions BED files
regions_bed_files <- list.files(mos_dir, pattern = ".regions.bed.gz$",
                                full.names = T)


## Output
genome_base <- basename(depth_out) # parse genome name from YGS results directory
depth_all_plot_file <- paste0("depth_all_", genome_base, ".png")
depth_cand_plot_file <- paste0("depth_candidate_", genome_base, ".png")
depth_diff_all_plot_file <- paste0("depth_diff_all_", genome_base, ".png")
depth_log2_all_plot_file <- paste0("depth_log2_all_", genome_base, ".png")
depth_log2_cand_plot_file <- paste0("depth_log2_candidate_", genome_base, ".png")
if (dir.exists(outdir)) {
  outdir <- file.path(outdir, "figures")
  dir.create(outdir, showWarnings = F, recursive = T)
  depth_all_plot_file <- file.path(outdir, depth_all_plot_file)
  depth_cand_plot_file <- file.path(outdir, depth_cand_plot_file)
  depth_diff_all_plot_file <- file.path(outdir, depth_diff_all_plot_file)
  depth_log2_all_plot_file <- file.path(outdir, depth_log2_all_plot_file)
  depth_log2_cand_plot_file <- file.path(outdir, depth_log2_cand_plot_file)
}


## Analysis
depth_df <- read_tsv(regions_bed_files, col_names = depth_cols,
                       id = "Sample") %>%
  mutate(Sample = gsub(".regions.bed.gz$", "", basename(Sample)),
         Sex = case_when(grepl("-FG-", Sample) ~ "female",
                         grepl("-MG-", Sample) ~ "male"),
         Chr = gsub("Chr", "Chr ", Chr),
         Chr = gsub("Chr 0", "Chr ", Chr),
         n_Chr = as.integer(str_split_i(Chr, " ", 2)),
         Midpoint = (Start + End) / 2,
         Window_Length = End - Start) %>%
  # Normalize per sample by genome-wide mean depth
  mutate(Sample_Genome_Mean = weighted.mean(Mean_Depth, w = Window_Length,
                                            na.rm = T),
         Normalized_Depth = Mean_Depth / Sample_Genome_Mean,
         .by = Sample)
# Normalize per sex by genome-wide mean depth
sex_df <- depth_df %>%
  summarize(Mean_Normalized_Depth = mean(Normalized_Depth, na.rm = T),
            SE_Normalized_Depth = sd(Normalized_Depth, na.rm = T) /
              sqrt(sum(!is.na(Normalized_Depth))),
            n = n_distinct(Sample),
            .by = c(Sex, n_Chr, Chr, Start, End, Midpoint))
# Calculate sex-based difference in depth per window
bias_df <- sex_df %>%
  summarize(Male_Mean_Depth = mean(Mean_Normalized_Depth[Sex == "male"], na.rm = T),
    Female_Mean_Depth = mean(Mean_Normalized_Depth[Sex == "female"], na.rm = T),
    Sex_Biased_Difference = Female_Mean_Depth - Male_Mean_Depth,
    Log2_Sex_Ratio = log2((Female_Mean_Depth + 0.05) / (Male_Mean_Depth + 0.05)),
    .by = c(n_Chr, Chr, Start, End, Midpoint)
  )
# Choose candidate sex chromosome(s)
depth_cands <- bias_df %>%
  summarize(n = n(),
            total_bias = sum(Log2_Sex_Ratio),
            mean_bias = total_bias/n,
            .by = n_Chr) %>%
  filter(mean_bias == max(mean_bias))
# Filter dataframes
filt_sex_df <- sex_df %>%
  filter(n_Chr %in% depth_cands$n_Chr)
filt_bias_df <- bias_df %>%
  filter(n_Chr %in% depth_cands$n_Chr)
# Title for sex chromosome candidate plot
if (length(depth_cands$n_Chr) == 1) {
  cand_ttl <- "Candidate sex chromosome"
} else {
  cand_ttl <- "Candidate sex chromosomes"
}


## Plots
# Plot normalized female and male mean read depth over windows
p_depth_all <- ggplot(sex_df, aes(x = Midpoint, y = Mean_Normalized_Depth,
                                     color = Sex)) +
  geom_line(linewidth = 0.3) +
  facet_wrap(vars(n_Chr), scale = "free_x",
             labeller = labeller(n_Chr = \(x) paste("Chr", x))) +
  scale_x_continuous(name = "Position (Mb)", labels = ~ .x * 1e-6) +
  scale_y_continuous(name = "Normalized mean depth per 100-kb window",
                     limits = c(0, 5))
# Plot difference between female and male mean read depth over windows
p_depth_diff_all <- ggplot(bias_df, aes(x = Midpoint,
                                          y = Sex_Biased_Difference)) +
  geom_line(linewidth = 0.3) +
  facet_wrap(vars(n_Chr), scale = "free_x",
             labeller = labeller(n_Chr = \(x) paste("Chr", x))) +
  scale_x_continuous(name = "Position (Mb)", labels = ~ .x * 1e-6) +
  scale_y_continuous(name = str_wrap(
    "Female vs. male raw difference in mean read depth per 100-kb window",
    width = 45),
    limits = c(-3, 3))
# Plot log2 ratio between female and male mean read depth over windows
p_log2ratio_all <- ggplot(bias_df, aes(x = Midpoint, y = Log2_Sex_Ratio,
                          color = Log2_Sex_Ratio)) +
  geom_point(size = 0.3) +
  facet_wrap(vars(n_Chr), scale = "free_x",
             labeller = labeller(n_Chr = \(x) paste("Chr", x))) +
  scale_x_continuous(name = "Position (Mb)", labels = ~ .x * 1e-6) +
  scale_y_continuous(name = str_wrap(
    "Female vs. male log2 fold change in mean read depth per 100-kb window",
    width = 45)) +
  scale_color_gradient2(name = str_wrap("female log2 ratio", width = 10),
                       high = "red", mid = "black", low = "blue")
# Plot for candidate chr normalized female and male mean read depth over windows
p_depth_cand <- ggplot(filt_sex_df, aes(x = Midpoint,
                                               y = Mean_Normalized_Depth,
                                               color = Sex)) +
  geom_line(linewidth = 0.3, show.legend = F) +
  facet_grid(cols = vars(n_Chr), rows = vars(Sex), scale = "free_x",
             labeller = labeller(n_Chr = \(x) paste("Chr", x))) +
  scale_x_continuous(name = "Position (Mb)", labels = ~ .x * 1e-6) +
  scale_y_continuous(name = "Normalized mean depth per 100-kb window",
                     limits = c(0, 5))
# Plot for candidate chr log2 ratio between female and male mean read depth over windows
p_log2ratio_cand <- ggplot(filt_bias_df, aes(x = Midpoint, y = Log2_Sex_Ratio,
                          color = Log2_Sex_Ratio)) +
  geom_point(size = 0.3) +
  facet_wrap(vars(n_Chr), scale = "free_x",
             labeller = labeller(n_Chr = \(x) paste("Chr", x))) +
  scale_x_continuous(name = "Position (Mb)", labels = ~ .x * 1e-6) +
  scale_y_continuous(name = str_wrap(
    "Female vs. male log2 fold change in mean read depth per 100-kb window",
    width = 45)) +
  scale_color_gradient2(name = str_wrap("female log2 ratio", width = 10),
                       high = "red", mid = "black", low = "blue")

# Save plots
ggsave(depth_all_plot_file, p_depth_all, dpi = 300, width = 11,
       height = 8)
ggsave(depth_cand_plot_file, p_depth_cand, dpi = 300, width = 7,
       height = 5)
ggsave(depth_diff_all_plot_file, p_depth_diff_all, dpi = 300, width = 11,
       height = 8)
ggsave(depth_log2_all_plot_file, p_log2ratio_all, dpi = 300, width = 11,
       height = 8)
ggsave(depth_log2_cand_plot_file, p_log2ratio_cand, dpi = 300, width = 7,
       height = 5)

cat("Plots saved to:", "\n")
cat(depth_all_plot_file, "\n")
cat(depth_cand_plot_file, "\n")
cat(depth_diff_all_plot_file, "\n")
cat(depth_log2_all_plot_file, "\n")
cat(depth_log2_cand_plot_file, "\n")

