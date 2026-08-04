## Initialization
# Load required packages
# library(ape)
library(tidytree)
library(treeio)
library(tidyverse)
library(ggtree)
library(ggpubr)
if (require(showtext)) {
  showtext_auto()
}

# Input
if(interactive()) {
  # Set working dir
  setwd("/scratch1/kdeweese/latissima_scratch/roscoff_comp")
  assembly_file <- "species_table.txt"
  tree_file <- "https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt"
  seq_file <- "s_lat_alignment.txt"
  outdir <- "./"
} else if(sourced && sourced == T) {
  tree_file <- "https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt"
  seq_file <- "s_lat_alignment.txt"
  if (dir.exists(outdir)) seq_file <- paste0(outdir, seq_file)
} else if(length(commandArgs(trailingOnly = T)) > 0) {
  line_args <- commandArgs(trailingOnly = T)
  assembly_file <- line_args[1]
  tree_file <- line_args[2]
  seq_file <- line_args[3]
  outdir <- line_args[4]
} else {
  stop("4 positional arguments expected.")
}

# Set filename(s)
extens <- c("png", "tiff")
plot_file <- "FS1_phylo_prune"
plot_file <- paste(plot_file, extens, sep = ".")
# Split filename for output naming
tree_name <- tools::file_path_sans_ext(basename(tree_file))
tree_ext <- tools::file_ext(tree_file)
out_tree <- paste0(tree_name, "_pruned.", tree_ext)
# Append output directory to output filenames (if it exists)
if (dir.exists(outdir)) {
  plot_file <- paste0(outdir, plot_file)
  out_tree <- paste0(outdir, out_tree)
}

# Format species to match tree names
unformatSpc <- function(spc) {
  spc <- gsub(" ", "_", spc)
  return(spc)
}
# Format species Latin name
formatSpc <- function(spc) {
  spc_f <- gsub("_", " ", spc)
  # Converts Ectocarpus siliculosus to Ectocarpus sp. Ec32
  spc_f <- gsub("Ectocarpus siliculosus", "Ectocarpus sp. Ec32", spc_f)
  # spc_f <- str_wrap(spc_f, width = 20)
  return(spc_f)
}

## Data
# Import species list
species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
colnames(species_tab) <- c("Species", "Assembly")
species_tab <- species_tab %>% select(Species, Assembly)
# Import tree
t <- read.tree(tree_file)
t <- as.treedata(t)
# Reformat tip labels and save original tree labels to dictionary
f_labs <- formatSpc(tip.label(t))
# Use reformatted labels to match to species
spc_match <- tibble(spc = species_tab$Species) %>%
  rowwise() %>%
  mutate(species = names(unlist(sapply(f_labs, grep, spc, value = T))))

# spc_match <- tibble(species = f_labs) %>%
#   rowwise() %>%
#   mutate(spc = case_when(length(grep(species, species_tab$Species)) > 0 ~
#                            paste(grep(species, species_tab$Species, value = T),
#                                  collapse = ";"),
#                          .default = NA))

# spc_match <- unlist(sapply(f_labs, grep, species_tab$Species, value = T))
# Error if not all species in tree
if (!(any(species_tab$Species %in% spc_match$spc))) {
# if (!(any(species_tab$Species %in% spc_match))) {
  stop("Error: not all species found in tree.")
}
labs_df <- tibble(labs=tip.label(t), f_labs=f_labs)
labs_df <- full_join(labs_df, spc_match, join_by(f_labs == species))
labs_df <- labs_df %>%
  mutate(label=unformatSpc(f_labs),
         highlight=case_when(!is.na(spc) ~ 1))
  # group_by(labs, f_labs, label, highlight) %>%
  # summarize(spc = paste(unique(spc), collapse = ";"), .groups = "drop")
# labs_df <- tibble(labs=tip.label(t), f_labs=f_labs) %>%
#   mutate(label=unformatSpc(f_labs),
#          spc=spc_match[f_labs],
#          highlight=case_when(!is.na(spc) ~ 1))
spc_dict <- labs_df %>%
  filter(!if_any(everything(), .fns = is.na)) %>%
  pull(label, spc)
# Create output species table
# (original species replaced with formatted tree labels)
species_tab2 <- species_tab %>%
  mutate(Species=spc_dict[Species])
# Relabel tree
tip_dict <- labs_df %>% pull(label, labs)
tip.label(t) <- unname(tip_dict[tip.label(t)])
for_merge <- labs_df %>%
  group_by(labs, f_labs, label, highlight) %>%
  summarize(spc = paste(unique(spc), collapse = ";"), .groups = "drop")
t <- full_join(t, for_merge, by = "label")
# t <- full_join(t, labs_df, by = "label")
# Prune tree
tp <- keep.tip(t, unique(unname(spc_dict)))
# write.tree(get.tree(tp), out_tree)
# print(paste("Pruned tree written to:", out_tree))
# # Generate formatted seqFile for Cactus, i.e.,
# # NEWICK tree
# # name1 path1
# # name2 path2
# # ...
# # nameN pathN
# write.tree(get.tree(tp), seq_file)
# write.table(species_tab2, seq_file,
#             append = T, sep = "\t",
#             col.names = F, row.names = F, quote = F)
# print(paste("Cactus-formated seqFile written to:", seq_file))
# 
# ## Plot phylogenetic trees
# t1 <- ggtree(t) +
#   hexpand(0.5) +
#   geom_tiplab(aes(label = f_labs, color = factor(highlight)),
#               fontface = "italic",
#               show.legend = F) +
#   scale_color_manual(values = "red", na.value = "black") +
#   theme_tree2()
# t2 <- ggtree(tp) +
#   hexpand(0.5) +
#   geom_tiplab(aes(label = f_labs), fontface = "italic") +
#   theme_tree2()
# p <- ggarrange(t1, t2, labels = "AUTO")
# sapply(plot_file, ggsave,
#        plot = p, height = 5, width = 10,
#        simplify = F)
# print(paste("Full and pruned phylogenetic trees plotted in:", plot_file))
