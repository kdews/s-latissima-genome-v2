## Cactus Pangenome Alignment Visualisation

# Clear environment
rm(list = ls())
# Required packages
library(gggenomes)
library(tidyverse)
library(patchwork)
library(scales)
if (require(showtext)) {
  showtext_auto()
  if (interactive()) {
    my_dpi <- dev.size("px")[1]/dev.size("in")[1]
    showtext_opts(dpi = my_dpi)
  }
  else {
    showtext_opts(dpi = 300)
  }
}

# Functions
# Extract numbers from contig IDs for filtering
fixChrom <- function(contigs) {
  contigs <- contigs %>%
    str_remove_all(".*_") %>%
    str_remove_all("[^0-9]") %>%
    as.numeric() %>%
    as.character()
  return(contigs)
}
# Creates dataframe of sequence IDs and their respective numbers,
# ordered by size only if chromosome IDs are not specified
chrNums <- function(cactus_tab, label, type = c("ref", "query")) {
  type <- match.arg(type)
  # Get FASTA index filenames for each assembly
  faidx <- cactus_tab %>%
    mutate(Index = paste0(Assembly, ".fai")) %>%
    filter(Label == label) %>%
    pull(Index)
  if (!file.exists(faidx)) stop(paste0("Index (,", faidx,") not found."))
  num_df <- read_tsv(faidx, col_names = c("ID", "Length"), show_col_types = F)
  num_df <- select(num_df, ID, Length)
  # Arrange "SL_" named chromosomes first by fixChrom extracted chr. number
  if (any(grepl("^SL_", num_df$ID))) {
    num_df <- num_df %>%
      mutate(fix_chrom = case_when(grepl("^SL_", ID) ~ fixChrom(ID))) %>%
      arrange(as.numeric(fix_chrom), desc(Length))
  }
  num_df <- num_df %>% 
    rownames_to_column("ID_n") %>%
    select(ID, ID_n, Length)
  # Name output ID column based on input type
  if (type == "ref") {
    num_df <- num_df %>%
      rename_with(~ gsub("ID", "seq_id2", .x)) %>%
      rename_with(~ gsub("Length", "length2", .x))
  } else if (type == "query") {
    num_df <- num_df %>%
      rename_with(~ gsub("ID", "seq_id", .x)) %>%
      rename_with(~ gsub("Length", "length", .x))
  }
  return(num_df)
}
# Create dataframe for dotplot from PAF data
dotDf <- function(paf, ref_nums, query_nums) {
  dot_data <- paf %>%
    mutate(
      x = start2,
      xend = end2,
      y = ifelse(strand == "+", start, end),
      yend = ifelse(strand == "+", end, start)
    ) %>%
    left_join(ref_nums) %>%
    left_join(query_nums)
  # Convert sequence ID number to factor for facet ordering
  chr_ord <- dot_data %>%
    mutate(chr_ord = as.numeric(seq_id_n)) %>%
    arrange(chr_ord) %>%
    distinct(seq_id_n) %>%
    pull(seq_id_n)
  chr2_ord <- dot_data %>%
    mutate(chr2_ord = as.numeric(seq_id2_n)) %>%
    arrange(chr2_ord) %>%
    distinct(seq_id2_n) %>%
    pull(seq_id2_n)
  dot_data <- dot_data %>%
    mutate(
      seq_id_n = factor(seq_id_n, levels = chr_ord),
      seq_id2_n = factor(seq_id2_n, levels = chr2_ord)
    )
  return(dot_data)
}

# Inputs
if (interactive()) {
  # Set working directory
  wd <- "/project2/noujdine_61/kdeweese/latissima/corteva_genome"
  setwd(wd)
  paf_file <- "cactus_pangenome_01-28-26_6039802/synteny_PAFs/merged_synteny.paf"
  ref <- "corteva_v2" # the --reference assembly prefix (as it appears in PAF tname)
  query <- "roscoff_v2" # second assembly name (as in PAF qname)
  seqFile <- "s_latissima_align.txt" # Cactus seqFile
} else {
  line_args <- commandArgs(trailingOnly = T)
  paf_file <- line_args[1]
  ref <- line_args[2]
  query <- line_args[3]
  seqFile <- line_args[4]
}

# Outputs
outdir <- "pangenome_plots"
chrom_dir <- file.path(outdir, "per_chromosome")
dir.create(outdir, showWarnings = F)
dir.create(chrom_dir, showWarnings = F)

# Import Cactus seqFile as table
cactus_tab <- read_tsv(seqFile, col_names = c("Label", "Assembly"))
# Tables of numbered sequence IDs
ref_nums <- chrNums(cactus_tab, ref, "ref")
query_nums <- chrNums(cactus_tab, query, "query")

# Read & filter PAF
# Standard PAF column names in gggenomes:
#   seq_id, length, start, end,
#   strand, seq_id2, length2, start2,
#   end2, map_match, map_length, map_quality
message("Reading PAF: ", paf_file, "\n")
# Suppress gggenomes::read_paf() warnings for records with fewer optional tags
paf_raw <- suppressWarnings(read_paf(paf_file))
message(sprintf("Raw alignments : %d", nrow(paf_raw)), "\n")

# PLOT 1: Coverage summary
cov <- paf_raw %>%
  group_by(seq_id2) %>%
  summarise(
    aligned_mb = sum(end2 - start2) / 1e6,
    n_contigs = n_distinct(seq_id),
    n_alignments = n(),
    .groups = "drop"
  ) %>%
  left_join(ref_nums) %>%
  mutate(frac_covered = pmin(aligned_mb * 1e6 / length2, 1)) %>%
  arrange(desc(length2))

p_cov <- ggplot(cov, aes(x = reorder(seq_id2_n, aligned_mb), y = aligned_mb)) +
  geom_col(aes(fill = n_contigs)) +
  scale_fill_viridis_c(name = "n contigs\naligned\n(France v2)", option = "plasma") +
  coord_flip() +
  labs(x = "North America v2", y = "Bases (Mb) aligned to France v2") +
  theme_bw(base_size = 11)

p_frac <- ggplot(cov, aes(x = reorder(seq_id2_n, frac_covered), y = frac_covered)) +
  geom_col(fill = "#41ab5d") +
  scale_y_continuous(labels = percent) +
  coord_flip() +
  labs(x = "North America v2", y = "Fraction covered by France v2 alignments") +
  theme_bw(base_size = 11)

p_hist <- ggplot(paf_raw, aes(x = map_length / 1e3)) +
  geom_histogram(bins = 60, fill = "#4393c3", color = "white", linewidth = 0.2) +
  scale_x_log10(labels = label_number()) +
  labs(x = "Alignment length (Kb)", y = "Count") +
  theme_bw(base_size = 11)

# Combine all stats plots
p_cov_panel <- (p_cov | p_frac) / p_hist

# Save stats plots to one image
if (interactive()) {
  print(p_cov_panel)
} else {
  showtext_opts(dpi = 300)
  ggsave(file.path(outdir, "coverage_summary.png"),
         p_cov_panel, width = 14, height = 12)
  showtext_opts(dpi = 100)
  message("Saved: coverage_summary.png", "\n")
}

# PLOT 2: Whole-genome dotplot
# Filter out short alignments and chrOther (non-chromosomes in ref)
paf <- paf_raw %>%
  # filter(substr(seq_id, 1, 2) == "SL") %>% # keep only chr (query)
  filter(grepl("^SL_", seq_id)) %>% # keep only chr (query)
  filter(!grepl("chrOther", seq_id2, ignore.case = T)) %>% # keep only chr (ref)
  filter(map_length >= 100e3)
message(sprintf("After filtering for chromosomes: %d alignments", nrow(paf)), "\n")

# Report sequence IDs that aligned
message("Aligned reference chromsomes", "\n")
print(sort(unique(paf$seq_id2)))
message("Aligned query chromosomes", "\n")
print(sort(unique(paf$seq_id)))

# Create dataframe for dotplot
dot_data <- dotDf(paf, ref_nums, query_nums)
dot_data_raw <- dotDf(paf_raw, ref_nums, query_nums)
# Get summary metadata from dataframe
n_df <- dot_data %>%
  summarize(ref = length(unique(seq_id2)), query = length(unique(seq_id)))
n_df_raw <- dot_data_raw %>%
  summarize(ref = length(unique(seq_id2)), query = length(unique(seq_id)))

# Generate whole-genome dotplot
p_dot <- ggplot(dot_data,
       aes(x = x, xend = xend, y = y, yend = yend, color = strand)) +
  geom_segment(linewidth = 0.8) +
  scale_color_manual(
    values = c("+" = "#2166ac", "-" = "#d6604d"),
    labels = c("+" = "Forward", "-" = "Reverse"),
    name   = "Strand"
  ) +
  facet_grid(
    rows = vars(seq_id2_n), # y-axis
    cols = vars(seq_id_n),  # x-axis
    scales = "free", # each facet has scaled axis range
    space  = "free", # facet size proportional to chromosome length
    switch = "both" # move
  ) +
  coord_flip() +
  labs(
    title = expression(
      italic("Saccharina latissima") * "whole-genome alignment"
    ),
    subtitle = sprintf("North America v2 (n = %d) vs France v2 (n = %d)",
                       n_df$ref, n_df$query),
    x = "North America v2",
    y = "France v2"
  ) +
  theme_linedraw() +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    strip.text = element_text(size = 7), # readable facet labels
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p_dot_raw <- ggplot(
  dot_data_raw, aes(x = x, xend = xend, y = y, yend = yend, color = strand)
) +
  geom_segment(linewidth = 0.8) +
  scale_color_manual(
    values = c("+" = "#2166ac", "-" = "#d6604d"),
    labels = c("+" = "Forward", "-" = "Reverse"),
    name   = "Strand"
  ) +
  facet_grid(
    rows = vars(seq_id2_n), # y-axis
    cols = vars(seq_id_n),  # x-axis
    scales = "free", # each facet has scaled axis range
    space  = "free", # facet size proportional to chromosome length
    switch = "both" # move
  ) +
  coord_flip() +
  labs(
    title = expression(
      italic("Saccharina latissima") * "whole-genome alignment"
    ),
    subtitle = sprintf("North America v2 (n = %d) vs France v2 (n = %d)",
                       n_df_raw$ref, n_df_raw$query),
    x = "North America v2",
    y = "France v2"
  ) +  
  theme_linedraw() +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    strip.text = element_text(size = 7), # readable facet labels
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Save whole-genome dotplot
if(interactive()) {
  print(p_dot)
  print(p_dot_raw)
} else {
  showtext_opts(dpi = 300)
  ggsave(file.path(outdir, "chr_whole_genome_dotplot.png"),
         p_dot, width = 14, height = 10)
  ggsave(file.path(outdir, "whole_genome_dotplot.png"),
         p_dot_raw, width = 14, height = 10)
  showtext_opts(dpi = 100)
  message("Saved: chr_whole_genome_dotplot.png and whole_genome_dotplot.png", "\n")
}

# PLOT 3: Per-chromosome dotplots
# One panel per North American chromosome, query seqs on y-axis
for (chrom in unique(dot_data$seq_id2)) {
  paf_chr <- paf %>% filter(seq_id2 == chrom)
  if (nrow(paf_chr) == 0) next
  
  # Only keep query seqs that have alignments to this chromosome
  qry_chr <- paf_chr %>%
    summarise(length = max(length), total = sum(map_length), .by = seq_id) %>%
    arrange(desc(total)) %>%
    mutate(offset = cumsum(lag(length, default = 0)))
  
  chr_len <- ref_nums %>% filter(seq_id2 == chrom) %>% pull(length2)
  total_qry_chr <- sum(qry_chr$length)
  
  dot_chr <- paf_chr %>%
    left_join(qry_chr %>% select(seq_id, qry_offset = offset),
              by = "seq_id") %>%
    mutate(
      x    = start2,
      xend = end2,
      y    = ifelse(strand == "+", start + qry_offset, end   + qry_offset),
      yend = ifelse(strand == "+", end   + qry_offset, start + qry_offset)
    )
  
  # Contig boundary lines on y-axis
  qry_bounds <- qry_chr %>% filter(offset > 0)
  
  p_chr <- ggplot(dot_chr,
                  aes(x = x, xend = xend, y = y, yend = yend, color = strand)) +
    geom_hline(data = qry_bounds, aes(yintercept = offset),
               color = "grey88", linewidth = 0.2) +
    geom_segment(alpha = 0.8, linewidth = 0.5) +
    scale_color_manual(
      values = c("+" = "#2166ac", "-" = "#d6604d"), name = "Strand"
    ) +
    scale_x_continuous(
      labels = label_number(scale = 1e-6, suffix = " Mb"),
      limits = c(0, chr_len), expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      labels = label_number(scale = 1e-6, suffix = " Mb"),
      limits = c(0, total_qry_chr), expand = c(0.01, 0)
    ) +
    labs(
      title    = chrom,
      subtitle = sprintf("%d aligned  |  %d unique alignments",
                         nrow(qry_chr), nrow(paf_chr)),
      x = "North America v2 position (Mb)",
      y = "France v2 cumulative position (Mb, by alignment weight)"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
  
  # Save chromosome alignment plot
  safe_nm <- gsub("[^A-Za-z0-9_]", "_", chrom)
  if (!interactive()) {
    showtext_opts(dpi = 300)
    ggsave(file.path(chrom_dir, paste0("chr_", safe_nm, ".png")),
           p_chr, width = 7, height = 5)
    showtext_opts(dpi = 100)
  }
}
for (chrom in unique(dot_data_raw$seq_id2)) {
  paf_chr <- paf_raw %>% filter(seq_id2 == chrom)
  if (nrow(paf_chr) == 0) next
  
  # Only keep query seqs that have alignments to this chromosome
  qry_chr <- paf_chr %>%
    summarise(length = max(length), total = sum(map_length), .by = seq_id) %>%
    arrange(desc(total)) %>%
    mutate(offset = cumsum(lag(length, default = 0)))
  
  chr_len <- ref_nums %>% filter(seq_id2 == chrom) %>% pull(length2)
  total_qry_chr <- sum(qry_chr$length)
  
  dot_chr <- paf_chr %>%
    left_join(qry_chr %>% select(seq_id, qry_offset = offset),
              by = "seq_id") %>%
    mutate(
      x    = start2,
      xend = end2,
      y    = ifelse(strand == "+", start + qry_offset, end   + qry_offset),
      yend = ifelse(strand == "+", end   + qry_offset, start + qry_offset)
    )
  
  # Contig boundary lines on y-axis
  qry_bounds <- qry_chr %>% filter(offset > 0)
  
  p_chr <- ggplot(dot_chr,
                  aes(x = x, xend = xend, y = y, yend = yend, color = strand)) +
    geom_hline(data = qry_bounds, aes(yintercept = offset),
               color = "grey88", linewidth = 0.2) +
    geom_segment(alpha = 0.8, linewidth = 0.5) +
    scale_color_manual(
      values = c("+" = "#2166ac", "-" = "#d6604d"), name = "Strand"
    ) +
    scale_x_continuous(
      labels = label_number(scale = 1e-6, suffix = " Mb"),
      limits = c(0, chr_len), expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      labels = label_number(scale = 1e-6, suffix = " Mb"),
      limits = c(0, total_qry_chr), expand = c(0.01, 0)
    ) +
    labs(
      title    = chrom,
      subtitle = sprintf("%d aligned  |  %d unique alignments",
                         nrow(qry_chr), nrow(paf_chr)),
      x = "North America v2 position (Mb)",
      y = "France v2 cumulative position (Mb, by alignment weight)"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
  
  # Save chromosome alignment plot
  safe_nm <- gsub("[^A-Za-z0-9_]", "_", chrom)
  if (!interactive()) {
    showtext_opts(dpi = 300)
    ggsave(file.path(chrom_dir, paste0(safe_nm, ".png")),
           p_chr, width = 7, height = 5)
    showtext_opts(dpi = 100)
  }
}
message("Saved per-chromosome dotplots to: ", chrom_dir, "\n")

# PLOT 4: Ribbon plots
# Build links table
# gggenomes links: seq_id/start/end = bin 1 (corteva, top)
#                  seq_id2/start2/end2 = bin 2 (roscoff, bottom)
# Create links for ribbon plot
links <- paf %>%
  rename(
    corteva_id    = seq_id2,
    corteva_start = start2,
    corteva_end   = end2,
    roscoff_id    = seq_id,
    roscoff_start = start,
    roscoff_end   = end
  ) %>%
  transmute(
    seq_id  = corteva_id,
    start   = corteva_start,
    end     = corteva_end,
    seq_id2 = roscoff_id,
    start2  = roscoff_start,
    end2    = roscoff_end,
    strand,
    map_length
  )

# Seq files
seqs_ref <- paf %>%
  group_by(seq_id = seq_id2) %>%
  summarise(length = max(length2), .groups = "drop") %>%
  mutate(bin_id = "North\nAmerica") %>%
  arrange(desc(length))
seqs_qry <- paf %>%
  group_by(seq_id) %>%
  summarise(length = max(length), .groups = "drop") %>%
  mutate(bin_id = "France") %>%
  arrange(desc(length))

# Find best match for each chromosome
best_match <- links %>%
  group_by(seq_id2, seq_id) %>%
  summarise(total = sum(map_length), .groups = "drop") %>%
  slice_max(total, n = 1, by = seq_id2) %>%
  arrange(seq_id2) %>%
  mutate(chr = fixChrom(seq_id2))
best_ref <- best_match %>%
  pull(chr, seq_id)
best_qry <- best_match %>%
  pull(chr, seq_id2)
seqs_ref <- seqs_ref %>%
  mutate(seq_id = as.character(best_ref[seq_id]))
seqs_qry <- seqs_qry %>%
  mutate(seq_id = as.character(best_qry[seq_id]))
links <- links %>%
  mutate(
    seq_id = as.character(best_ref[seq_id]),
    seq_id2 = as.character(best_qry[seq_id2])
  )
# Apply ordering to seqs_all
seqs_all <- bind_rows(seqs_ref, seqs_qry) %>%
  arrange(bin_id, as.integer(seq_id))

p_ribbon <-
  gggenomes(seqs = seqs_all, links = links, spacing = 0.2) +
  geom_seq() +
  geom_bin_label(size = 3) +
  geom_seq_label(nudge_y = -0.05) +
  geom_link(aes(fill = strand, color = strand), offset = 0.04, curve = 15) +
  scale_fill_manual(values = c("+" = "#2166ac", "-" = "#d6604d"), name = "Strand") +
  scale_colour_manual(values = c("+" = "#2166ac", "-" = "#d6604d"), guide = "none") +
  labs(
    title = expression(italic("Saccharina latissima") * " chromosome synteny")
  )
p_built <- ggplot_build(p_ribbon)
x_range <- p_built$layout$panel_params[[1]]$x.range
y_range <- p_built$layout$panel_params[[1]]$y.range
p_ribbon <- p_ribbon +
  coord_cartesian(
    xlim = c(x_range[1]*0.6, x_range[2]),
    ylim = c(y_range[1]*3.5, y_range[2]*0.75),
    expand = F
  )
if (interactive()) {
  print(p_ribbon)
} else {
  showtext_opts(dpi = 300)
  ggsave(file.path(outdir, "whole_genome_ribbon.png"),
         p_ribbon, width = 7, height = 5)
  showtext_opts(dpi = 100)
  message(
    "Saved whole-genome ribbon plot to: ",
    outdir, "whole_genome_ribbon.png", "\n"
  )
}

# Summary tables
sum_df <- paf_raw %>%
  group_by(ref_chrom = seq_id2, qry_contig = seq_id) %>%
  summarise(n = n(), total_kb = round(sum(map_length)/1e3, 1), .groups = "drop") %>%
  arrange(ref_chrom, desc(total_kb))
chr_sum_df <- paf %>%
  group_by(ref_chrom = seq_id2, qry_contig = seq_id) %>%
  summarise(n = n(), total_kb = round(sum(map_length)/1e3, 1), .groups = "drop") %>%
  arrange(ref_chrom, desc(total_kb))
# Write TSV
if(interactive()) {
  print("Skipping write TSV step.")
} else {
  write.table(
    sum_df,
    file.path(outdir, "alignment_summary.tsv"),
    sep = "\t", quote = F, row.names = F
  )
  write.table(
    chr_sum_df,
    file.path(outdir, "chr_alignment_summary.tsv"),
    sep = "\t", quote = F, row.names = F
  )
  message(
    "Wrote summary tables: ",
    file.path(outdir, "alignment_summary.tsv"),
    file.path(outdir, "chr_alignment_summary.tsv"),
    "\n"
  )
  message("Outputs in: ", outdir, "\n")
}

