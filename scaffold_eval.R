# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(scales)
if (require(showtext)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Inputs
if (interactive()) {
  wd <- "/project2/noujdine_61/kdeweese/latissima/corteva_genome"
  setwd(wd)
  genome_file <- "s_latissima_genome_list.txt"
} else {
  line_args <- commandArgs(trailingOnly = T)
  genome_file <- line_args[1]
}
# Order for label factor
lbl_order <- c("roscoff_v1",
               "roscoff_v2",
               "lindell_v0",
               "lindell_v1",
               "corteva_v2")
# Output plot filenames
vio_plot <- "scaffold_sizes_violin_log.png"
bar_plot <- "scaffold_sizes_bar.png"

# Functions
# Get version from assembly labels
getVersion <- function(label) {
  version <- gsub(".*(v[0-9]).*", "\\1", label)
  # Label polished assemblies
  has_polish <- grepl("polish", label, ignore.case = T)
  version[has_polish] <- paste0(version[has_polish], ".p")
  return(version)
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
    group_by(Label, Region, Version) %>%
    summarize(N50=Biostrings::N50(`Length (Mb)`),
              L50=sum(`Length (Mb)` > N50),
              total=sum(`Length (Mb)`),
              n=n())
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
  # size_fun <- function(x) max(density(x)$x)*1.1
  # Calculate N50 from original size distribution, then convert to log10
  n50_line_fun <- function(x) log10(Biostrings::N50(10^x))
  n50_fun <- function(x) n50_line_fun(x)*1.08
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
  p <- ggplot(data = idx, mapping = aes(x = Label, y = Length)) +
    # Format x-axis labels
    scale_x_discrete(labels = getVersion, name = "Assembly version") +
    # Convert y-axis from bp to log10 bp scale
    scale_y_log10(labels = label_log()) +
    labs(y = "Length (bp)") +
    facet_wrap(~ Region, scales = "free_x") +
    # Violin: equalize violin width between groups, lower alpha for points
    geom_violin(
      mapping = aes(
        col = Label,
        fill = Label
      ),
      scale = vio_scale,
      trim = F,
      linewidth = 1,
      alpha = 0.4,
      show.legend = F
    ) +
    # Points: shrink size and jitter width (not height) for visibility
    # geom_jitter(height = 0, width = 0.1, alpha = 0.05) + # size = 0.2, ) +
    # ggforce::geom_sina(mapping = aes(col = Version),
    #                    scale = vio_scale, alpha = 0.8, show.legend = F) +
                       # , alpha = Length
    geom_boxplot(width = 0.1) +
    # Label assembly total size and n()
    geom_text(
      data = sum_idx,
      mapping = aes(x = levels(Label), y = max_y * 5),
      label = size_lab
    ) +
    labs(title = expression(italic("Saccharina latissima") * " genomes")) +
    theme_bw() +
    theme(
      # axis.text.x = element_text(face = "italic"), # italicize species
      strip.text = element_text(size = 14),
      panel.ontop = F,
      # panel.border = element_rect(fill = NA)
    ) +
    coord_cartesian(clip = "off")
  # Label N50/L50 on top of violins and mark N50 with dashed line on graph
  p_n50 <- p +
    geom_crossbar(
      data = sum_idx,
      aes(
        x = Label,
        y = N50 * 1e6,
        ymin = N50 * 1e6,
        ymax = N50 * 1e6
      ),
      width = 0.8,
      linetype = "dashed",
      linewidth = 0.1,
      alpha = 0.8
    ) +
    geom_label(
      data = sum_idx,
      aes(
        x = Label,
        y = N50* 1e6,
        label = n50_lab
      ),
      position = position_nudge(x = -0.35, y = 0.3),
      size = 3,
      fill = "white"
      # fill = NA
    )
    # stat_summary(
    #   mapping = aes(x = as.numeric(Label) - 0.45, xend = as.numeric(Label) + 0.45),
    #   geom = "segment",
    #   linetype = "dashed",
    #   alpha = 0.8,
    #   fun = n50_line_fun
    # ) +
    # stat_summary(
    #   mapping = aes(x = as.numeric(Label) - 0.33),
    #   geom = "label",
    #   label = n50_lab,
    #   fill = "white",
    #   size = 2.5,
    #   fun = n50_fun
    # )
  if (missing(n50)) {
    return(p)
  } else if (n50) {
    return(p_n50)
  } else {
    print("Error: unrecognized argument to 'annot' variable.")
  }
}
# Create annotated graphs of contig length distribution by assembly version
distPlot <- function(idx, limit = 1e6) {
  idx <- idx %>%
    arrange(Label, desc(Length)) %>%
    group_by(Label) %>%
    mutate(n=n(),
           ID=1:unique(n),
           n50=Biostrings::N50(Length),
           high=case_when(Length >= n50 ~ "over_N50"),
           cutoff=case_when(Length >= limit ~ paste0(">", limit),
                            .default = paste0("<", limit))) %>%
    ungroup() %>%
    mutate(Label=factor(Label, levels = rev(levels(Label))))
  p <- ggplot(data = idx, mapping = aes(x = ID, y = Length)) +
    # geom_hline(aes(yintercept=n50), col = "red", lty = 2, ) +
    # geom_col(aes(col = high, fill = high), show.legend = F) +
    geom_hline(yintercept = limit, col = "red", lty = 2) +
    geom_col(aes(col = cutoff, fill = cutoff), show.legend = F) +
    scale_color_manual(name = "", aesthetics = c("color", "fill"),
                       values = c("lightgrey", "darkblue")) +
    # Convert length to log10 scale
    scale_y_log10(labels = label_log()) +
    # scale_y_log10(labels = label_log(), expand = c(0, 0)) +
    # scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Scaffolds + contigs", y = "Length (bp)") +
    facet_wrap(~Label, scales = "free_x") +
    theme_bw() +
    theme(
      # strip.text = element_text(face = "italic"), # italicize species
      strip.background = element_blank()
    )
  return(p)
}

# Read in data
genome_tab <- read_tsv(genome_file, col_names = F, comment = "#")
colnames(genome_tab) <- c("Label", "Assembly")
genome_tab <- genome_tab %>%
  mutate(
    Region = case_when(
      grepl("corteva|lindell", Label) ~ "North America",
      grepl("roscoff", Label) ~ "France"
    ),
    Version = getVersion(Label),
    .after = Label
  )
fns <- paste0(genome_tab$Assembly, ".fai")
for (i in 1:length(fns)) {
  idx_temp <- read_tsv(fns[i], col_names = F, col_select = c(X1, X2))
  idx_temp <- idx_temp %>%
    rename(ID = X1, Length = X2) %>%
    mutate(
      Region = genome_tab %>% slice(i) %>% pull(Region),
      Version = genome_tab %>% slice(i) %>% pull(Version),
      Label = genome_tab %>% slice(i) %>% pull(Label)
    )
  if (i == 1) {
    idx <- idx_temp
  } else {
    idx <- rbind(idx, idx_temp)
  }
}
lbl_lvls <- sapply(lbl_order, grep, unique(idx$Label), value = T) %>%
  unname() %>%
  unlist()
idx <- idx %>%
  # Sort "Label" column with factor for plotting
  mutate(Label = factor(Label, levels = lbl_lvls)) %>%
  group_by(Label) %>%
  arrange(desc(Length), .by_group = T) %>%
  # Convert ID characters into numeric
  mutate(ID_num = as.numeric(str_remove_all(str_remove_all(ID, ".*_"),
                                            "[^0-9]")))
# Summarize data frame
sum_idx <- sumDf(idx)

# Plots
# Plot length distributions using violin plots
p_l <- violinPlot(idx, sum_idx, n50 = T)
# p_d <- distPlot(idx, 1e6)
# Save plots
print(paste("Saving violin plot to:", vio_plot))
showtext_opts(dpi = 300)
ggsave(filename = vio_plot, plot = p_l, bg = "white", width = 9, height = 7)
showtext_opts(dpi = 100)
# print(paste("Saving barplot plot to:", bar_plot))
# ggsave(filename = bar_plot, plot = p_d, bg = "white", width = 10, height = 6)

# Compare all v2 assemblies by region 
comp_v2_idx <- idx %>%
  filter(
    grepl("v2", Version),
    !grepl("bacteria", ID, ignore.case = T)
  ) %>%
  mutate(Region = factor(Region, levels = c("North America", "France"))) %>%
  group_by(Label) %>%
  arrange(desc(Length), .by_group = T) %>%
  mutate(
    `Length (Mb)` = Length / 1e6,
    `Chromosomes + contigs` = as.numeric(factor(ID, levels = unique(ID))),
    Type = case_when(row_number() <= 31 ~ "chromosome", .default = "contig")
  )
# Plot v2 assemblies showing chromosome & contig split (log10-scale)
p_comp_v2 <- ggplot(comp_v2_idx, aes(x = `Chromosomes + contigs`, y = Length)) +
  geom_vline(xintercept = 31, linetype = "dashed", alpha = 0.4) +
  geom_point(aes(color = Region), size = 2) +
  annotate(geom = "label", x = 31, y = max(comp_v2_idx$Length), label = "n = 31") +
  scale_x_log10() +
  scale_y_log10(name = "Length (bp)", label = label_log()) +
  labs(title = expression(
    "Chromosome length in " * italic("S. latissima") * " (v2) assemblies"
  ))
# Save plot of region comparison
showtext_opts(dpi = 300)
ggsave(
  filename = "comp_v2_length.png",
  plot = p_comp_v2,
  width = 7,
  height = 5
)
showtext_opts(dpi = 100)

# Compare polished and unpolished North American assemblies
comp_polish_idx <- comp_v2_idx %>%
  filter(Region == "North America")
p_comp_polish <- 
  ggplot(comp_polish_idx, aes(x = `Chromosomes + contigs`, y = Length)) +
  geom_col(aes(fill = Version), position = "dodge") +
  scale_fill_discrete(palette = scales::pal_brewer(palette = "Paired")) +
  scale_x_continuous(expand = c(0.01, 0.01), breaks = c(1:31, seq(35, 120, 5))) +
  facet_wrap(~ Type, ncol = 1, scales = "free") +
  labs(
    x = NULL,
    y = "Length (bp)",
    title = expression(
    "Impact of polishing on length (North America " * italic("S. latissima") * " )"
    )
  ) +
  theme_bw()
# Save plot of polished comparison
showtext_opts(dpi = 300)
ggsave(
  filename = "comp_polish_length.png",
  plot = p_comp_polish,
  width = 10,
  height = 10
)
showtext_opts(dpi = 100)

# Save table of labeled chromosome and contig IDs for v2 assemblies
v2_contig_type <- comp_v2_idx %>%
  ungroup() %>%
  select(Label, ID, Type) %>%
  mutate(Filename = paste0(Label, "-", Type, "_ids.txt"))
for (contig_file in unique(v2_contig_type$Filename)) {
  tmp_df <- v2_contig_type %>%
    filter(Filename == contig_file) %>%
    select(ID)
  write_tsv(x = tmp_df, file = contig_file, col_names = F)
}
