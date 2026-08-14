# Clear environment
rm(list = ls())
# Set timezone
Sys.setenv(TZ = "America/Los_Angeles")
# Load required packages
library(readr)
library(dplyr)
library(GENESPACE)

# Input
# Only take command line input if not running interactively
if (interactive()) {
  # Set working dir
  wd <- "/project2/noujdine_61/kdeweese/latissima/corteva_genome"
  setwd(wd)
  in_tsv <- "s-latissima-genome-v2/genomes_metadata.tsv"
  outdir <- "s-latissima-genome-v2"
} else {
  line_args <- commandArgs(trailingOnly = T)
  cat("Arguments:", "\n")
  cat(line_args, "\n", sep = "\n")
  in_tsv <- line_args[1]
  outdir <- line_args[2]
}

genomeRepo <- "/scratch1/kdeweese/corteva_genome/assemblies"
wd <- "/scratch1/kdeweese/corteva_genome"
path2mcscanx <- "~/miniforge3/envs/genespace/bin/MCScanX"

genome_tab <- read_tsv(in_tsv, col_names = T, comment = "")
colnames(genome_tab) <- gsub("# ", "", colnames(genome_tab))

# Derive assembly paths
genome_tab <- genome_tab %>%
  mutate(Assembly = file.path("assemblies", paste0(Label, ".fa")))

genomes2run <- genome_tab$Label_real

