## Initialization
# Clear environment
rm(list = ls())
# Load required packages
suppressPackageStartupMessages(library(tidyverse))

## Script variables
ygs_cols <- c(
  GI = "Window ID",
  NUM = "n",
  MAX_K = "Ungapped Window Length (bp)",
  K = "Total Distinct k-mers",
  UK = "Unmatched k-mers",
  SC_K = "Single-Copy k-mers",
  SC_UK = "Single-Copy Unmatched k-mers",
  P_SC_UK = "Percent Single-Copy Unmatched k-mers",
  VSC_K = "Validated Single-Copy k-mers",
  VSC_UK = "Validated Single-Copy Unmatched k-mers",
  P_VSC_UK = "Percent Validated Single-Copy Unmatched k-mers"
)

## Input
# Only take command line input if not running interactively
if (interactive()) {
  # Set working directory
  wd <- "/scratch2/kdeweese/corteva_genome"
  setwd(wd)
  # YGS results directory
  ygs_out <- "sdr_id_results/ygs_sdr_results/corteva_v2_ha"
  # Output directory
  outdir <- "s-latissima-genome-v2"
} else {
  line_args <- commandArgs(trailingOnly = T)
  ygs_out <- line_args[1]
  outdir <- line_args[2]
}
# Locate final results file
ygs_res_file <- list.files(ygs_out, pattern = ".final_result", full.names = T)


## Output
genome_base <- basename(ygs_out) # parse genome name from YGS results directory
ygs_all_plot_file <- paste0("ygs_all_", genome_base, ".png") # plot image file
ygs_cand_plot_file <- paste0("ygs_candidate_", genome_base, ".png") # plot image file
if (dir.exists(outdir)) {
  outdir <- file.path(outdir, "figures")
  dir.create(outdir, showWarnings = F, recursive = T)
  ygs_all_plot_file <- file.path(outdir, ygs_all_plot_file)
  ygs_cand_plot_file <- file.path(outdir, ygs_cand_plot_file)
}


## Analysis
lines <- trimws(readLines(ygs_res_file))
first_row <- grep("^>", lines)[1]
lines <- paste(c(lines[first_row - 1], grep("^>", lines, value = T)),
               collapse = "\n")
ygs_df <- read_table(lines, na = ".", show_col_types = F) %>%
  rename_with(~ ygs_cols[.x]) %>% # add detailed column names
  relocate(n) %>% # move n column to first
  # Parse window IDs
  separate_wider_delim(`Window ID`, delim = ":",
                       names = c("Chr", "Coordinates"), cols_remove = F) %>%
  separate_wider_delim(Coordinates, delim = "-", names = c("Start", "End")) %>%
  mutate(Start = as.numeric(Start),
         End = as.numeric(End),
         Chr = gsub("^>|_sliding$", "", Chr),
         Chr = factor(Chr, levels = unique(Chr))) %>%
  mutate(n_Chr = cur_group_id(), .by = Chr, .before = Chr) %>%
  mutate(`Ungapped Chr Length (bp)` = sum(`Ungapped Window Length (bp)`),
         .by = Chr, .after = `Ungapped Window Length (bp)`)
# Choose candidate sex chromosome(s)
ygs_cands <- ygs_df %>%
  summarize(n = n(),
            total_YGS = sum(`Percent Validated Single-Copy Unmatched k-mers`),
            mean_YGS = total_YGS/n,
            .by = n_Chr) %>%
  filter(mean_YGS == max(mean_YGS))
filt_ygs_df <- ygs_df %>%
  filter(n_Chr %in% ygs_cands$n_Chr)
# Title for sex chromosome candidate plot
if (length(ygs_cands$n_Chr) == 1) {
  cand_ttl <- "Candidate sex chromosome"
} else {
  cand_ttl <- "Candidate sex chromosomes"
}


## Plots
# Plot windows as line segments faceted by chromosome
p_all <- ggplot(ygs_df,
       aes(x = Start,
           y = `Percent Validated Single-Copy Unmatched k-mers`,
           color = `Percent Validated Single-Copy Unmatched k-mers`)) +
  geom_segment(aes(xend = End,
                   yend = `Percent Validated Single-Copy Unmatched k-mers`),
               linewidth = 1) +
  facet_wrap(vars(n_Chr), scale = "free_x",
             labeller = labeller(n_Chr = \(x) paste("Chr", x))) +
  scale_x_continuous(name = "Position (Mb)", labels = ~ .x * 1e-6) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_gradient(low = "black", high = "red", name = "%") +
  theme_bw()
# Plot for candidate sex chr only
p_cand <- ggplot(filt_ygs_df,
       aes(x = Start,
           y = `Percent Validated Single-Copy Unmatched k-mers`,
           color = `Percent Validated Single-Copy Unmatched k-mers`)) +
  geom_segment(aes(xend = End,
                   yend = `Percent Validated Single-Copy Unmatched k-mers`),
               linewidth = 1) +
  facet_wrap(vars(n_Chr), scale = "free_x",
             labeller = labeller(n_Chr = \(x) paste("Chr", x))) +
  scale_x_continuous(name = "Position (Mb)", labels = ~ .x * 1e-6) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_gradient(low = "black", high = "red", name = "%") +
  labs(title = cand_ttl) +
  theme_bw()

# Save plots
ggsave(ygs_all_plot_file, p_all, dpi = 300, width = 11, height = 8)
ggsave(ygs_cand_plot_file, p_cand, dpi = 300, width = 7, height = 5)
cat("Plots saved to:", "\n")
cat(ygs_all_plot_file, "\n")
cat(ygs_cand_plot_file, "\n")

