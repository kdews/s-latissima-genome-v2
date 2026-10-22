## Initialization
# Clear environment
rm(list = ls())
# Load required packages
# BiocManager::install("rtracklayer")
library(rtracklayer)
library(tidyverse)

## Functions
# Fix any attribute keys that are "start"
fix_gff_attr <- function(gff_file) {
  fix_gff <- readLines(gff_file)
  fix_gff <- sub("(^|;)start=", "\\1complete_start=", fix_gff)
  return(fix_gff)
}
# Use temporary file to import modified lines as GFF
import_gff3 <- function(gff_file) {
  lines <- fix_gff_attr(gff_file)
  tmp <- tempfile(fileext = ".gff3")
  on.exit(unlink(tmp))
  writeLines(lines, tmp)
  # gff3 <- import.gff3(tmp)
  # return(gff3)
  readGFF(tmp, version = 3)
}

## Input
annot_dir <- "/scratch2/kdeweese/corteva_genome/annotations"
annot_files <- list.files(annot_dir, full.names = T)
# read_tsv(annot_files, comment = "#", col_names = F)

ecto_gff <- import_gff3(annot_files[1])
ecto_df <- ecto_gff %>%
  as.data.frame() %>%
  as_tibble()

sj_gff <- import_gff3(annot_files[2])
sj_df <- sj_gff %>%
  as.data.frame() %>%
  as_tibble()

us_v1_gff <- import_gff3(annot_files[3])
us_v1_df <- us_v1_gff %>%
  as.data.frame() %>%
  as_tibble()

eu_v1_gff <- import_gff3(annot_files[[4]])
eu_v1_df <- eu_v1_gff %>%
  as.data.frame() %>%
  as_tibble()
read.table(annot_files[4], sep = "\t") %>% sapply(., class)
read.table(annot_files[5], sep = "\t") %>% sapply(., class)

eu_v2_gff <- import_gff3(annot_files[5])
eu_v2_df <- eu_v2_gff %>%
  as.data.frame() %>%
  as_tibble()

read_tsv(annot_files[1], col_names = F, comment = "#")
read_tsv(annot_files[4], col_names = F, comment = "#")
read_tsv(annot_files[5], col_names = F, comment = "#")
