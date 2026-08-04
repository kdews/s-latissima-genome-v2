## Initialization
# Load required packages
library(tidyverse)

## Input
if (interactive()) {
  setwd("/project2/noujdine_61/kdeweese/latissima/corteva_genome")
  line_args <- c(
    "s-latissima-genome-v2/dna_reads_PRJEB72149.tsv",
    "seq_divergence"
  )
} else if (length(commandArgs(trailingOnly = T)) == 2) {
  line_args <- commandArgs(trailingOnly = T)
} else {
  stop("2 positional arguments expected.")
}
ena_tsv <- line_args[1]
outdir <- line_args[2]

# Output parsed TSV
out_tsv <- file.path(outdir, paste0("parsed_", basename(ena_tsv)))

# Parse ENA TSV
ena_df <- read_tsv(ena_tsv) %>%
  separate_longer_delim(everything(), delim = ";") %>%
  mutate(
    fastq_ftp = paste0("ftp://", fastq_ftp),
    fastq = paste(
      scientific_name,
      word(sample_title, -1),
      gsub("(.*)\\.fastq.*$", "\\1", basename(fastq_ftp)),
      sep = "_"
    ),
    fastq = gsub("\\.", "", fastq),
    fastq = gsub(" ", "_", fastq),
    fastq = paste0(fastq, gsub(".*(\\.fastq.*)$", "\\1", basename(fastq_ftp))),
    fastq = file.path(outdir, fastq)
  ) %>%
  select(fastq_ftp, fastq_md5, fastq)

# Write table of FTPs and expected md5s
write_tsv(ena_df, out_tsv, col_names = F)
