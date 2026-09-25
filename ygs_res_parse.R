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
ygs_res_file <- "sdr_id_results_1/corteva_v2_500kb_os_reads_k15_ss_reads_k15.final_result"

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
filt_ygs_df <- ygs_df %>%
  filter(n_Chr <= 31) %>%
  filter(any(`Percent Validated Single-Copy Unmatched k-mers` >= 50, na.rm = T),
         .by = n_Chr)

# Use 90th percentile to color graph
cutoff <- quantile(filt_ygs_df$`Percent Validated Single-Copy Unmatched k-mers`,
                   0.90, na.rm = T)

ggplot(filt_ygs_df,
       aes(x = Start,
           y = `Percent Validated Single-Copy Unmatched k-mers`,
           color = `Percent Validated Single-Copy Unmatched k-mers` >= cutoff)) +
  geom_segment(aes(xend = End,
                   yend = `Percent Validated Single-Copy Unmatched k-mers`),
               linewidth = 1,
               show.legend = F) +
  # geom_point(show.legend = F) +
  facet_wrap(~ n_Chr, scale = "free_x") +
  scale_x_continuous(name = "Position (bp)",
                     labels = ~ .x * 1e-3) +
  scale_color_manual(values = c("black", "red")) +
  # coord_flip() +
  theme_minimal()
  # theme(
  #   axis.ticks.y = element_blank(),
  #   axis.text.y = element_blank(),
  #   axis.title.y = element_blank()
  # )
  




