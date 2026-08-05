# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
library(ape)
library(phangorn)
library(Biostrings)
library(msa)
library(treeio)
library(ggtree)
library(ggmsa)
library(ggpubr)
# library(ggtreeExtra)
# library(ggtext)

set.seed(100) # set for reproducibility

# Input
# Only take command line input if not running interactively
if (interactive()) {
  # Set working dir
  wd <- "/project2/noujdine_61/kdeweese/latissima/corteva_genome"
  setwd(wd)
  seq_file <- "s-latissima-genome-v2/s_latissima_prog_align.txt"
  lineage <- "stramenopiles_odb12.2"
  outdir <- "s-latissima-genome-v2"
} else {
  line_args <- commandArgs(trailingOnly = T)
  cat("Arguments:", "\n")
  cat(line_args, "\n", sep = "\n")
  seq_file <- line_args[1]
  lineage <- line_args[2]
  outdir <- line_args[3]
}
pretty_labs_file <- "pretty_genome_labels.tsv" # pretty label table
# HMM directory in local BUSCO downloads
hmms_dir <- paste("busco_downloads/lineages", lineage, "hmms", sep = "/")
out_grp <- "ectocarpus_liu2024" # outgroup: Ectocarpus

# Output
outname <- tools::file_path_sans_ext(basename(seq_file)) # name from input file
og_dist_file <- "unique_og_dist.png"
nj_nwk_file <- "NJ_tree.txt"
ml_nwk_file <- "ML_tree.txt"
trees_plot_file <- "NJ_ML_trees.png"
msa_plot_file <- "msa_ML_tree.png"
og_dist_file <- paste(outname, og_dist_file, sep = "-")
nj_nwk_file <- paste(outname, nj_nwk_file, sep = "-")
ml_nwk_file <- paste(outname, ml_nwk_file, sep = "-")
trees_plot_file <- paste(outname, trees_plot_file, sep = "-")
msa_plot_file <- paste(outname, msa_plot_file, sep = "-")
if (dir.exists(outdir)) {
  nj_nwk_file <- file.path(outdir, nj_nwk_file)
  ml_nwk_file <- file.path(outdir, ml_nwk_file)
  og_dist_file <- file.path(outdir, "figures", og_dist_file)
  trees_plot_file <- file.path(outdir, "figures", trees_plot_file)
  msa_plot_file <- file.path(outdir,  "figures", msa_plot_file)
}

# Multithreading cores
my_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(my_cores)) my_cores <- as.integer(1)
cat("Cores set to:", my_cores, "\n")

# Read Cactus seqFile into data frame
seq_tab <- read_tsv(seq_file,
                    col_names = c("Label", "Assembly"),
                    comment = "(") # skip guide tree (if one exists)
# Annotate with pretty labels
pretty_labs <- read_tsv(pretty_labs_file)
seq_tab <- left_join(seq_tab, pretty_labs) %>%
  mutate(Pretty_Label = gsub(" ", "_", Pretty_Label)) # replace spaces
# Create dictionaries correlating labels and pretty labels
dict_to_pretty <- seq_tab %>%
  distinct(Label, Pretty_Label) %>%
  pull(Pretty_Label, Label)
dict_to_reg <- seq_tab %>%
  distinct(Label, Pretty_Label) %>%
  pull(Label, Pretty_Label)
# Read in list of all lineage ortholog IDs from HMM directory
og_ids <- sub(
  "\\.hmm$",
  "",
  list.files(path = hmms_dir, pattern = "\\.hmm$")
) %>%
  tibble(OG = .) %>%
  separate_wider_delim(OG, delim = "at", names = c("OG", "Lineage_ID")) %>%
  distinct(OG)
# Construct paths to BUSCO single-copy FASTAs
seq_tab <- seq_tab %>%
  mutate(
    BUSCO_SC = paste0(
      "busco_results/",
      lineage,
      "/",
      Label,
      "/run_",
      lineage,
      "/busco_sequences/single_copy_busco_sequences"
    )
  )
n_asm <- length(unique(seq_tab$Label)) # number of assemblies
# Check each BUSCO single-copy result directory exists
for (sc_dir in seq_tab$BUSCO_SC) {
  if (!file.exists(sc_dir)) {
    stop(paste0("BUSCO single-copy seq. directory (", sc_dir, ") not found."))
  }
}
sc_list <- sapply(
  seq_tab$BUSCO_SC,
  list.files,
  pattern = ".faa",
  full.names = T,
  simplify = F
)
sc_list <- sapply(sc_list,
                  as_tibble_col,
                  column_name = "OG_FAA",
                  simplify = F)
sc_df <- sc_list %>%
  bind_rows(.id = "BUSCO_SC") %>%
  # Join to assembly labels
  left_join(seq_tab) %>%
  select(-BUSCO_SC) %>%
  # Extract BUSCO orthogroup (OG) IDs
  mutate(OG = gsub("\\.faa", "", basename(OG_FAA)), .before = 1) %>%
  mutate(OG = paste0("OG", OG)) %>%
  separate_wider_delim(OG, delim = "at", names = c("OG", "Lineage_ID")) %>%
  # Filter for OG present in all assemblies
  filter(n() == length(seq_tab$Label), .by = OG) %>%
  # Read in predicted protein sequences of each OG from each assembly
  rowwise() %>%
  mutate(Protein_Seq = as.character(readAAStringSet(OG_FAA))) %>%
  ungroup() %>%
  # Save protein IDs separately
  mutate(Protein_ID = names(Protein_Seq),
         Protein_Seq = unname(Protein_Seq)) %>%
  # Strip terminal stop codons (*)
  mutate(Protein_Seq = gsub("\\*$", "", Protein_Seq)) %>%
  # Filter out OGs with internal stop codons
  filter(!any(grepl("\\*", Protein_Seq)), .by = OG) %>%
  # Parse information from protein IDs
  separate_wider_delim(Protein_ID,
                       delim = "|",
                       names = c("OG_Hit", "Hit", "Strand")) %>%
  mutate(Strand = gsub("^.*([+-]).*$", "\\1", Strand)) %>%
  separate_wider_delim(Hit, delim = ":", names = c("Chr", "Coordinates")) %>%
  separate_wider_delim(Coordinates,
                       delim = "-",
                       names = c("Start", "End"))

# Plot distribution of unique protein sequences per orthogroup
cat("Plotting distribution of unique BUSCOs across assemblies per OG...", "\n")
n_ogs <- length(unique(sc_df$OG)) # number of present orthogroups
n_odb_ogs <- length(unique(og_ids$OG)) # number of all orthogroups in ODB
p_dist <- sc_df %>%
  summarize(n_in_og = length(unique(Protein_Seq)), .by = OG) %>%
  ggplot() +
  geom_histogram(aes(x = n_in_og), bins = n_asm, col = "white") +
  annotate(geom = "text", x = n_asm/2, y = n_odb_ogs/3,
           label = paste("n =", n_ogs, "/", n_odb_ogs)) +
  labs(
    title = str_replace(lineage, "_odb", " ODB") %>%
      paste("Divergence of single-copy orthologs in", .),
    subtitle = gsub("_", " ", seq_tab$Pretty_Label) %>%
      gsub("Saccharina", "S.", .) %>% 
      paste(collapse = ", ") %>%
      str_wrap(., width = 80),
    x = "assemblies with unique protein sequence (per orthogroup)",
    y = "orthogroups"
  ) +
  theme(plot.subtitle = element_text(face = "italic"))
if (interactive()) {
  print(p_dist)
} else {
  ggsave(og_dist_file, p_dist)
  cat("Plot saved to:", og_dist_file, "\n")
}

# Align each OG
cat("Aligning", n_ogs, "orthogroups...", "\n")
sc_aln <- list()
sc_aa <- list()
for (i in 1:n_ogs) {
  cat(i, "/", n_ogs, "\n")
  og <- unique(sc_df$OG)[[i]]
  prots <- sc_df %>%
    filter(OG == og) %>%
    arrange(Label) %>%
    pull(Protein_Seq, Label) %>%
    AAStringSet()
  aln <- msa(prots,
             method = "Muscle",
             type = "protein",
             order = "input",
             verbose = F)
  sc_aln[[og]] <- aln
  sc_aa[[og]] <- msaConvert(aln, type = "ape::AAbin")
}
# Concatenate all alignments into supermatrix
super <- do.call(cbind, sc_aa)
cat("Alignment supermatrix length:", dim(super)[2], "bp", "\n")
super_phy <- as.phyDat(super) # convert to phyDat object

# Generate assembly phylogeny from single-copy BUSCO protein alignments
# Tree 1: Neighbor-joining (NJ)
cat("Creating neighbor-joining (NJ) tree...", "\n")
dm <- dist.ml(super_phy)
treeNJ <- NJ(dm)
# Root tree to outgroup
treeNJ_r <- root(treeNJ, outgroup = out_grp, resolve.root = T)
# Save rooted tree
cat("Writing NJ Newick tree...", "\n")
if (interactive()) {
  write.tree(treeNJ_r)
} else {
  write.tree(treeNJ_r, nj_nwk_file)
  cat("Newick tree saved to:", nj_nwk_file, "\n")
}
# Bootstrap (on unrooted)
cat("Boostrapping NJ tree...", "\n")
bsNJ <- bootstrap.phyDat(super_phy, function(x) NJ(dist.ml(x)))
treeNJ_r_bs <- plotBS(treeNJ_r, bsNJ, main = "NJ + standard bootstrap",
                      type = "none") # get bootstrap node labels
# Plot tree
cat("Plotting NJ tree phylogeny...", "\n")
pNJ <- data.frame(
  node = 1:Nnode(treeNJ_r_bs) + Ntip(treeNJ_r_bs),
  bootstrap = as.numeric(treeNJ_r_bs$node.label) * 100 # add bootstrap as df
) %>%
  full_join(treeNJ_r_bs, ., by = "node") %>%
  ggtree() +
  geom_tiplab(aes(label = dict_to_pretty[label] %>%
                    str_replace_all("_", " ") %>%
                    str_replace_all(" \\(", "\n\\(")
                  ),
              fontface = "italic") +
  geom_nodepoint(aes(fill = cut(bootstrap, c(0, 70, 90, 100))),
                 shape = 21, size = 4) +
  theme_tree(legend.position = c(0.8, 0.2)) +
  scale_fill_manual(
    values = c("white", "grey", "black"),
    guide = "legend",
    name = "Bootstrap (%)",
    breaks = c('(90,100]', '(70,90]', '(0,70]'),
    labels = expression(BP >= 90, 70 <= BP * " < 90", BP < 70)
  ) +
  labs(title = "Neighbor-joining (NJ) tree") +
  hexpand(0.2) +
  theme_tree2()

# Tree 2: Maximum-likelihood (ML)
cat("Creating maximum-likelihood (ML) tree...", "\n")
# Automatically choose best model
if (interactive()) {
  mt <- modelTest(
    super_phy,
    # tree = treeNJ
  )
} else {
  mt <- modelTest(
    super_phy,
    # tree = treeNJ,
    control = pml.control(trace = 1) # increase verbosity
  )
}
cat("Fitting best model...", "\n")
if (interactive()) {
  fit_mt <- pml_bb(
    mt,
    # start = treeNJ
  )
} else {
  fit_mt <- pml_bb(
    mt,
    # start = treeNJ,
    control = pml.control(trace = 1) # increase verbosity
  )
}
treeML <- fit_mt$tree # extract tree from model result
# Root tree to outgroup
treeML_r <- root(fit_mt$tree, outgroup = out_grp, resolve.root = T)
# Strip node labels
treeML_r_tmp <- treeML_r 
treeML_r_tmp$node.label <- NULL
# Save rooted tree
cat("Writing ML Newick tree...", "\n")
if (interactive()) {
  write.tree(treeML_r_tmp)
} else {
  write.tree(treeML_r_tmp, ml_nwk_file)
  cat("Newick tree saved to:", ml_nwk_file, "\n")
}
# Bootstrap (on unrooted)
cat("Boostrapping ML tree...", "\n")
if (interactive()) {
  bsML <- bootstrap.pml(fit_mt, optNni = T)
} else {
  bsML <- bootstrap.pml(
    fit_mt,
    optNni = T,
    control = pml.control(trace = 2) # increase verbosity
  )
}
treeML_r_bs <- plotBS(treeML_r, bsML, main = "ML + standard bootstrap",
                      type = "none") # get bootstrap node labels
# Plot tree
cat("Plotting ML tree phylogeny...", "\n")
pML <- data.frame(
  node = 1:Nnode(treeML_r_bs) + Ntip(treeML_r_bs),
  bootstrap = as.numeric(treeML_r_bs$node.label) * 100 # add bootstrap as df
) %>%
  full_join(treeML_r_bs, ., by = "node") %>%
  ggtree() +
  geom_tiplab(aes(label = dict_to_pretty[label] %>%
                    str_replace_all("_", " ") %>%
                    str_replace_all(" \\(", "\n\\(")
                  ),
              fontface = "italic") +
  geom_nodepoint(aes(fill = cut(bootstrap, c(0, 70, 90, 100))),
                 shape = 21, size = 4) +
  theme_tree(legend.position = c(0.8, 0.2)) +
  scale_fill_manual(
    values = c("white", "grey", "black"),
    guide = "legend",
    name = "Bootstrap (%)",
    breaks = c('(90,100]', '(70,90]', '(0,70]'),
    labels = expression(BP >= 90, 70 <= BP * " < 90", BP < 70)
  ) +
  labs(title = "Maximum-likelihood (ML) tree") +
  hexpand(0.2) +
  theme_tree2()

# Plot both trees
cat("Combining tree plots...", "\n")
p <- ggarrange(pNJ, pML, common.legend = T, legend = "right")
if (interactive()) {
  print(p)
} else {
  ggsave(trees_plot_file, p, width = 12, height = 7, bg = "white")
  cat("Plot saved to:", trees_plot_file, "\n")
}

# Update seqFile with ML guide tree
seq_lines <- c(
  write.tree(treeML_r_tmp),
  paste(seq_tab$Label, seq_tab$Assembly, sep = "\t")
)
cat("Constructing Cactus seqFile with ML tree...", "\n")
cat(paste(seq_lines, collapse = "\n"), "\n")
cat("Wrote Cactus seqFile to:", seq_file, "\n")
writeLines(seq_lines, seq_file)

# Plot tree + MSA supermatrix
cat("Plotting MSA with tree...", "\n")
if (!all(treeML_r_bs$tip.label %in% names(super_phy))) {
  stop("Names in ML tree don't match alignment names.")
}
# Convert phyDat supermatrix to AAStringSet
super_aa <- as.character(super_phy) %>%
  apply(1, paste, collapse = "") %>%
  AAStringSet()
ptree <- ggtree(treeML_r_bs) +
  geom_tiplab(aes(label = dict_to_pretty[label] %>%
                    str_replace_all("_", " ") %>%
                    str_replace_all(" \\(", "\n\\(")
                  ),
              fontface = "italic") +
  geom_nodelab()
if (interactive()) {
  pMSA <- msaplot(ptree, super_aa, offset = 0.1, width = 3, window=c(150, 200)) +
    labs(title = "ML tree from concatenated BUSCO protein alignment") +
    guides(fill = "none")  # don't display MSA fill legend
  print(pMSA)
} else {
  pMSA <- msaplot(ptree, super_aa, offset = 0.1, width = 3) +
    labs(title = "ML tree from concatenated BUSCO protein alignment") +
    guides(fill = "none")  # don't display MSA fill legend
  ggsave(msa_plot_file, pMSA, bg = "white", height = 4, width = 15)
  cat("Plot saved to:", msa_plot_file, "\n")
}
