## Initialization
# Clear environment
rm(list = ls())
# Load required packages
library(tidyverse)

## Input
if (interactive()) {
  setwd("/project2/noujdine_61/kdeweese/latissima/corteva_genome")
}
corteva_v2_fai_file <- "assemblies/corteva_v2.fa.fai"
ha_polish_corr_file <- "US_v2_HA_polish_chrs.tsv"
out_tsv <- "US_v2_HA_polish_correspondence.tsv"

## Analysis
old_df <- read_table(corteva_v2_fai_file, col_names = F)
new_df <- read_tsv(ha_polish_corr_file)
old_df <- old_df %>%
  select(X1) %>%
  rename(OLD_CORR = X1) %>%
  rowid_to_column(var = "OLD_ID")
corr_df <- new_df %>%
  left_join(old_df, by = join_by(NEW_CORR == OLD_CORR))
write_tsv(corr_df, out_tsv)
