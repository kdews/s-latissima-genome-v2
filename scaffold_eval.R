# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(scales)

## Input
if (interactive()) {
  wd <- "/project2/noujdine_61/kdeweese/latissima/corteva_genome"
  setwd(wd)
  in_tsv <- "s-latissima-genome-v2/genomes_metadata.tsv" # input genome TSV
  scripts_dir <- "s-latissima-genome-v2"
} else {
  line_args <- commandArgs(trailingOnly = T)
  in_tsv <- line_args[1]
  scripts_dir <- line_args[2]
}
# Order for label factor
lbl_order <- c(
  "latissima.*US",
  "latissima.*France",
  "^(?!.*(US|France)).*latissima",
  "^(?!.*latissima)"
)
# Order for region factor
rg_order <- c("US", "France", "Norway", "China", "Peru")

## Output
# Plot filenames
vio_plot <- "scaffold_sizes_violin_log.png"
bar_plot <- "scaffold_sizes_bar.png"
comp_v2_plot <- "comp_v2_length.png"
comp_polish_plot <- "comp_polish_length.png"
# Prepend output directory to file names (if it exists)
outdir <- file.path(scripts_dir, "figures")
if (dir.exists(outdir)) {
  vio_plot <- file.path(outdir, vio_plot)
  bar_plot <- file.path(outdir, bar_plot)
  comp_v2_plot <- file.path(outdir, comp_v2_plot)
  comp_polish_plot <- file.path(outdir, comp_polish_plot)
}

## Functions
# Format species Latin name
formatSpc <- function(spc) {
  spc_f <- gsub("_", " ", spc)
  if (!grepl("Ecto", spc_f)) {
    genus <- word(spc_f, 1, 1)
    species <- word(spc_f, 2, 2)
    genus_abv <- paste0(substr(genus, 1, 1), ".")
    spc_f <- paste(genus_abv, species)
  }
  return(spc_f)
}
# Extract numbers from contig IDs for filtering
fixChrom <- function(contigs) {
  contigs <-
    as.character(as.numeric(str_remove_all(str_remove_all(contigs, ".*_"),
                                           "[^0-9]")))
  return(contigs)
}
# Find N_len-length gaps (N's) from FASTA file using Biostrings
findGaps <- function(fasta_file, N_len) {
  library(Biostrings)
  fasta <- readDNAStringSet(fasta_file)
  gap_ptn <- paste(rep("N", N_len), collapse = "")
  gap_matches <- vmatchPattern(gap_ptn, fasta)
  start_comp <- startIndex(gap_matches)
  fasta_gaps <- tibble(start_comp) %>%
    mutate(seqid_comp = names(fasta)) %>%
    unnest_longer(start_comp) %>%
    mutate(start_comp = as.numeric(start_comp),
           end_comp = start_comp + N_len,
           gap_length = N_len)
  # %>%
  #   # Convert genomic position columns from bp to Mb
  #   mutate_at(.vars = vars(grep("start|end|length", colnames(.), value = T)),
  #             .funs = ~ .x*1e-6)
  return(fasta_gaps)
}
# Summarize index statistics for each assembly
sumDf <- function(df) {
  sum_df <- df %>%
    # Convert "Length" column from bp to Mb
    mutate(`Length (Mb)` = Length*1e-6) %>%
    group_by(Label, Pretty_Label, Species, Region, Version) %>%
    summarize(N50 = Biostrings::N50(`Length (Mb)`),
              L50 = sum(`Length (Mb)` > N50),
              total = sum(`Length (Mb)`),
              n = n(),
              .groups = "drop")
  return(sum_df)
}
# Create annotated violin plot of contig lengths by assembly version
violinPlot <- function(idx, sum_idx, n50 = NULL) {
  # Positioning functions for annotations
  # Calculate maximum y in density distributions of each assembly version
  find_max_y <- function(label)
    idx %>% filter(paste0(Region, Version) == label) %>%
    pull(Length) %>%
    log10() %>%
    density() %>%
    .$x %>%
    max()
  max_y <- 10^(max(sapply(unique(
    pull(mutate(idx, label = paste0(Region, Version)), label)
  ), find_max_y)))
  # Plot annotations
  # Total size and number of scaffolds+contigs
  size_lab <- paste(paste(round(sum_idx$total, digits = 1), "Mb"),
                    paste("n =", prettyNum(sum_idx$n, big.mark = ",")),
                    sep = "\n")
  # N50/L50 of each assembly (in Mb)
  n50_lab <- paste(paste("N50 =", round(sum_idx$N50, digits = 1), "Mb"),
                   paste("L50 =", sum_idx$L50),
                   sep = "\n")
  # Scale to use for violin and sina plots
  vio_scale <- "width"
  # Version dictionary
  vers_dict <- sum_idx %>% pull(Version, Label)
  p <- ggplot(data = idx, mapping = aes(x = Label, y = Length)) +
    # Violin: equalize violin width between groups, lower alpha for points
    geom_violin(mapping = aes(color = Label, fill = Label), scale = vio_scale,
                trim = F, linewidth = 1, alpha = 0.4, show.legend = F) +
    # # Points: use ggforce to shape points within violins
    # ggforce::geom_sina(mapping = aes(fill = Label, alpha = Length),
    #                    scale = vio_scale) +
    # Boxplot: overlay statistical summary
    geom_boxplot(width = 0.1) +
    facet_wrap(~ Species + Region, scales = "free_x", space = "free_x",
               nrow = 1) +
    # Label assembly total size and n()
    geom_text(data = sum_idx, mapping = aes(x = Label, y = max_y * 5),
              label = size_lab) +
    # Format x-axis labels
    scale_x_discrete(labels = vers_dict, name = "Assembly version") +
    # Convert y-axis from bp to log10 bp scale
    scale_y_log10(labels = label_log()) +
    ylab("Length (bp)") +
    theme_bw() +
    theme(strip.text = element_text(size = 14, face = "italic"),
          panel.ontop = F, panel.border = element_blank(),
          panel.background = element_rect(colour = "black", fill = NA)) +
    coord_cartesian(clip = "off")
  # Label N50/L50 on top of violins and mark N50 with dashed line on graph
  p_n50 <- p +
    geom_crossbar(data = sum_idx,
                  aes(x = Label, y = N50 * 1e6, ymin = N50 * 1e6,
                      ymax = N50 * 1e6),
      linetype = "dashed", linewidth = 0.1, width = 0.8, alpha = 0.8) +
    geom_label(data = sum_idx, aes(x = Label, y = N50 * 1e6, label = n50_lab),
               position = position_nudge(x = -0.35, y = 0.3),
               size = 3, fill = "white")
  if (missing(n50)) {
    return(p)
  } else if (n50) {
    return(p_n50)
  } else {
    cat("Error: unrecognized argument to 'annot' variable.\n")
  }
}
# Create annotated graphs of contig length distribution by assembly version
distPlot <- function(idx, limit = 1e6) {
  if (is.numeric(limit)) {
    var_type <- "limit"
    legend_title <- c(paste("<", limit/1e6, "Mb"), paste(">", limit/1e6, "Mb"))
  } else if (is.character(limit) & limit == "n50") {
    var_type <- "n50"
    legend_title <- c("< n50", "> n50")
  } else {
    cat("Error: unrecognized argument to 'limit' variable.",
          "Expected either numeric or 'n50'.\n")
  }
  dist_idx <- idx %>%
    arrange(Label, desc(Length)) %>%
    group_by(Label) %>%
    mutate(n = n(), ID = 1:unique(n),
           limit = limit, n50 = Biostrings::N50(Length)) %>%
    ungroup()
  sum_dist_idx <- dist_idx %>%
    group_by(Label, Pretty_Label, Species, Region, Version) %>%
    summarize(x_pos = max(ID) * 0.9, n50 = unique(n50),
              limit = unique(limit), .groups = "drop")
  ggplot(data = dist_idx, mapping = aes(x = ID, y = Length)) +
    # Horizontal line at limit
    geom_hline(aes(yintercept = .data[[var_type]]), color = "red", lty = 2) +
    # Label limit
    geom_text(data = sum_dist_idx, color = "red",
              mapping = aes(x = x_pos, y = .data[[var_type]] * 4,
                            label = paste(round(.data[[var_type]]/1e6, 2),
                                          "Mb"))) +
    # Log-scale barplot of scaffold and contig lengths
    geom_col(aes(col = Length >= .data[[var_type]],
                 fill = Length >= .data[[var_type]])) +
    # Set colors above and below cutoff
    scale_color_manual(name = "Length", values = c("lightgrey", "darkblue"),
                       aesthetics = c("color", "fill"), labels = legend_title,
                       guide = guide_legend(reverse = T)) +
    # Remove space to the left of the bars
    scale_x_continuous(expand = c(0, 0)) +
    # Convert length to log10 scale
    scale_y_log10(labels = label_log()) +
    facet_wrap(~ Pretty_Label, scales = "free_x") +
    labs(x = "", y = "Length (bp)") +
    theme_bw() +
    theme(strip.text = element_text(face = "italic"), # italicize species
          strip.background = element_blank(),
          # legend.direction = "vertical",
          legend.position = "top", legend.title.position = "left")
}

# Read in data
genome_tab <- read.table(in_tsv, sep = "\t", header = T, fill = NA,
                         comment.char = "", check.names = F)
colnames(genome_tab) <- gsub("#| ", "", colnames(genome_tab))
genome_tab <- genome_tab %>%
  # Keep only labels
  select(Label, Pretty_Label) %>%
  # Derive assembly paths
  mutate(Assembly = file.path("assemblies", paste0(Label, ".fa"))) %>%
  mutate(
    # Mark polished and draft assemblies
    AssemblyType = case_when(
      grepl("polish", Pretty_Label) ~ ".p",
      grepl("draft", Pretty_Label) ~ ".d",
      .default = ""
    ),
    # Get region from assembly labels
    Region = gsub("^.*\\(([A-Za-z]+)[ ;].*\\)$", "\\1", Pretty_Label, perl = T),
    Region = factor(Region, levels = rg_order),
    .after = Pretty_Label
  ) %>%
  rowwise() %>%
  mutate(
    # Parse version from assembly labels
    Version = case_when(
      # Parse vN
      grepl("v[0-9]", Pretty_Label) ~ gsub("^.*(v[0-9]).*$", "\\1", Pretty_Label),
      # Parse year of cited publication
      grepl("et al\\.", Pretty_Label) ~
        gsub("^.*et al\\.\\ ([0-9]+)\\).*$", "\\1", Pretty_Label),
      .default = "none"
    ),
    Version = paste0(Version, AssemblyType),
    # Parse and factor species names
    Species = gsub(" \\(.*$", "", Pretty_Label),
    Species = formatSpc(Species),
    Species = factor(Species, levels = unique(Species)),
    .after = Pretty_Label
  ) %>%
  ungroup()
fns <- paste0(genome_tab$Assembly, ".fai")
for (i in 1:length(fns)) {
  idx_temp <- read_tsv(fns[i], col_names = F, col_select = c(X1, X2),
                       show_col_types = F)
  idx_temp <- idx_temp %>%
    rename(ID = X1, Length = X2) %>%
    mutate(Label = genome_tab %>% slice(i) %>% pull(Label))
  if (i == 1) {
    idx <- idx_temp
  } else {
    idx <- rbind(idx, idx_temp)
  }
}
# Sort "Pretty_Label" column with factor for plotting
lbl_lvls <- sapply(regex(lbl_order),
                   grep,
                   unique(genome_tab$Pretty_Label),
                   value = T,
                   perl = T) %>%
  unname() %>%
  unlist()
idx <- idx %>%
  # Merge by "Label" column
  left_join(genome_tab, by = join_by(Label)) %>%
  # Factor "Pretty_Label"
  mutate(Pretty_Label = factor(Pretty_Label, levels = lbl_lvls)) %>%
  arrange(Pretty_Label) %>%
  # Factor "Label" to match "Pretty_Label"
  mutate(Label = factor(Label, levels = unique(Label))) %>%
  group_by(Pretty_Label) %>%
  arrange(desc(Length), .by_group = T) %>%
  # Convert ID characters into numeric
  mutate(ID_num = as.numeric(str_remove_all(str_remove_all(ID, ".*_"),
                                            "[^0-9]")))
# Summarize data frame
sum_idx <- sumDf(idx)

# Plot assembly length distributions
p_l <- violinPlot(idx, sum_idx, n50 = T)
p_d <- distPlot(idx, 1e6)
p_d_n50 <- distPlot(idx, "n50")

# Compare all v2 assemblies by region 
comp_v2_idx <- idx %>%
  filter(Species == "S. latissima", Version == "v2",
         !grepl("bacteria", ID, ignore.case = T)) %>%
  group_by(Label) %>%
  arrange(desc(Length), .by_group = T) %>%
  mutate(n = as.numeric(factor(ID, levels = unique(ID))),
         Type = case_when(row_number() <= 31 ~ "chromosome",
                          .default = "contig"))
# Save tables of labeled v2 chromosome IDs for Minigraph-Cactus pangenome
v2_contig_type <- comp_v2_idx %>%
  ungroup() %>%
  select(Label, ID, Type) %>%
  filter(Type == "chromosome") %>%
  mutate(Filename = paste0(Label, "-", Type, "_ids.txt"))
for (contig_file in unique(v2_contig_type$Filename)) {
  tmp_df <- v2_contig_type %>%
    filter(Filename == contig_file) %>%
    select(ID)
  write_tsv(x = tmp_df, file = contig_file, col_names = F)
}
# Plot v2 assemblies showing chromosome & contig split (log10-scale)
p_comp_v2 <- ggplot(comp_v2_idx, aes(x = n, y = Length)) +
  geom_vline(xintercept = 31, linetype = "dashed", alpha = 0.4) +
  geom_point(aes(color = Region), size = 2, alpha = 0.8) +
  annotate(geom = "label", x = 31, y = max(comp_v2_idx$Length),
           label = "n = 31") +
  scale_x_log10() +
  scale_y_log10(name = "Length (bp)", label = label_log()) +
  labs(title = expression("Chromosome length in " * italic("S. latissima") *
                            " v2 assemblies")) +
  theme_bw() +
  theme(legend.position = "inside", legend.justification = c(0.9, 0.9),
        legend.background = element_rect(color = "darkgrey"))

# Compare polished and unpolished US assemblies
comp_polish_idx <- idx %>%
  filter(Species == "S. latissima", Region == "US",
         Version %in% c("v2", "v2.p")) %>%
  group_by(Label) %>%
  arrange(desc(Length), .by_group = T) %>%
  mutate(n = as.numeric(factor(ID, levels = unique(ID))),
         Type = case_when(row_number() <= 31 ~ "chromosomes",
                          .default = "contigs"))
p_comp_polish <- ggplot(comp_polish_idx,
                        aes(x = n, y = Length)) +
  geom_col(aes(fill = Version), position = "dodge") +
  scale_fill_discrete(palette = scales::pal_brewer(palette = "Paired")) +
  scale_x_continuous(expand = c(0.01, 0.01), breaks = c(1:31, seq(35, 120, 5))) +
  facet_wrap(~ Type, ncol = 1, scales = "free") +
  labs(x = NULL, y = "Length (bp)",
       title = expression("Impact of polishing on length (US " * 
                            italic("S. latissima") * ")")) +
  theme_bw()

# Save plots
# Violin plot
cat("Saving assembly violin plots to:", vio_plot, "\n")
vio_width <- 9 * dim(sum_idx)[1]/5
ggsave(filename = vio_plot, plot = p_l, bg = "white", width = vio_width, height = 7)
cat("Saving contig length barplot plots to:", bar_plot, "\n")
ggsave(filename = bar_plot, plot = p_d, bg = "white", width = 10, height = 6)
# US vs. France v2 comparison
cat("Saving US vs. France v2 comparison plot to:", comp_v2_plot, "\n")
ggsave(filename = comp_v2_plot, plot = p_comp_v2, width = 7, height = 5)
# US polishing comparison
cat("Saving US polishing comparison plot to:", comp_polish_plot, "\n")
ggsave(filename = comp_polish_plot, plot = p_comp_polish, width = 10,
       height = 10)
