# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)

# Inputs
if (interactive()) {
  wd <- "/project2/noujdine_61/kdeweese/latissima/corteva_genome"
  setwd(wd)
}

pretty_labs_file <- "pretty_genome_labels.tsv"
pretty_labs <- c(
  "lindell_v0" = "Saccharina latissima (US v0)",
  "lindell_v1" = "Saccharina latissima (US v1)",
  "corteva_v1" = "Saccharina latissima (US v2 draft)",
  "corteva_v2" = "Saccharina latissima (US v2)",
  "corteva_v2_polished" = "Saccharina latissima (US v2 polished)",
  "roscoff_v1" = "Saccharina latissima (France v1)",
  "roscoff_v2" = "Saccharina latissima (France v2)",
  "ilvo" = "Saccharina latissima (Norway; Bråtelund et al. 2024)",
  "japonica_li2025" = "Saccharina japonica (Li et al. 2025)",
  "ectocarpus_liu2024" = "Ectocarpus sp. (Liu et al. 2024)"
)
pretty_df <- tibble(Label = names(pretty_labs), Pretty_Label = pretty_labs)
write_tsv(pretty_df, pretty_labs_file)


