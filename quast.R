## Initialization
# Load required packages
library(scales)
library(tidyverse)
# if (require(showtext)) {
#   showtext_auto()
#   if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
# }

# Input
# Only take command line input if not running interactively
if (interactive()) {
  setwd("/project2/noujdine_61/kdeweese/latissima/corteva_genome/")
  assembly_file <- "s_latissima_genome_list.txt"
  # Output directory
  outdir <- "./"
}
pretty_labs_file <- "pretty_genome_labels.tsv" # pretty label table
n_content_file <- "N_content_per_100_kb.png"
# Prepend output directory to file names (if it exists)
if (dir.exists(outdir)) n_content_file <- paste0(outdir, n_content_file)
# Order for label factor
lbl_order <- c(
  "latissima.*US",
  "latissima.*France",
  "^(?!.*(US|France)).*latissima",
  "^(?!.*latissima)"
)
# Order for region factor
rg_order <- c("US", "France", "Norway", "Other")

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

# Data wrangling
# Import data frame of assembly labels and FASTAs
genome_tab <- read.table(assembly_file, sep = "\t", header = T, fill = NA,
                         comment.char = "", check.names = F)
colnames(genome_tab) <- gsub("#| ", "", colnames(genome_tab))
# Annotate with pretty labels
pretty_labs <- read_tsv(pretty_labs_file)
genome_tab <- left_join(pretty_labs, genome_tab)
# Sort "Pretty_Label" column with factor for plotting
lbl_lvls <- sapply(regex(lbl_order),
                   grep,
                   unique(genome_tab$Pretty_Label),
                   value = T,
                   perl = T) %>%
  unname() %>%
  unlist()
genome_tab <- genome_tab %>%
  mutate(
    # Factor "Pretty_Label"
    Pretty_Label = factor(Pretty_Label, levels = rev(lbl_lvls)),
    # Mark polished and draft assemblies
    AssemblyType = case_when(
      grepl("polish", Pretty_Label) ~ ".p",
      grepl("draft", Pretty_Label) ~ ".d",
      .default = ""
    ),
    # Get region from assembly labels
    Region = case_when(
      grepl("US", Pretty_Label, ignore.case = F) ~ "US",
      grepl("France", Pretty_Label, ignore.case = F) ~ "France",
      grepl("Norway", Pretty_Label, ignore.case = F) ~ "Norway",
      .default = "Other"
    ),
    Region = factor(Region, levels = rg_order),
    # Add QUAST report names
    QUAST_report = paste0("quast_", Label, "/transposed_report.tsv"),
    .after = Label
  ) %>%
  rowwise() %>%
  mutate(
    # Get version from assembly labels
    Version = case_when(
      # Get vN
      grepl("v[0-9]", Pretty_Label) ~ gsub("^.*(v[0-9]).*$", "\\1", Pretty_Label),
      # Get year of cited publication
      grepl("et al\\.", Pretty_Label) ~
        gsub("^.*et al\\.\\ ([0-9]+)\\).*$", "\\1", Pretty_Label),
      .default = "none"
    ),
    Version = paste0(Version, AssemblyType),
    Species = gsub(" \\(.*$", "", Pretty_Label),
    Species = formatSpc(Species),
    .after = Label
  ) %>%
  ungroup()
# Verify each QUAST report exists
for (quast_file in genome_tab$QUAST_report) {
  if (!file.exists(quast_file)) {
    stop("QUAST output directory not found: ", quast_file)
  }
}
quast_df <- read_tsv(genome_tab$QUAST_report, na = "-") %>%
  rename(Label = Assembly) %>%
  pivot_longer(-Label, names_to = "Stat")
quast_filt <- quast_df %>%
  filter(!grepl("_broken$", Label)) %>%
  left_join(genome_tab, .) %>%
  filter(!grepl(">=", Stat))
quast_filt %>%
  pivot_wider(names_from = Stat, values_from = value) %>%
  select(Label, 9:last_col()) %>%
  View()
n_content <- quast_filt %>% filter(grepl("^# N", Stat))

# Plotting
# Get dictionary for plot labels
lbl_dict <- genome_tab %>%
  mutate(Plot_Label = paste(Species, Version)) %>%
  pull(Plot_Label, Pretty_Label)
# Plot N content
p_n <- ggplot(data = n_content,
       mapping = aes(x = value,
                     y = Pretty_Label,
                     fill = Pretty_Label, col = Pretty_Label)) +
  geom_col(show.legend = F) +
  facet_grid(
    rows = vars(Region),
    scales = "free",
    space = "free",
    switch = "both",
  ) +
  scale_y_discrete(labels = lbl_dict) +
  xlab(unique(n_content$Stat)) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "italic"),
    legend.text = element_text(face = "italic"),
    strip.placement = "outside"
  )
ggsave(p_n, filename = n_content_file, bg = "white", width = 9, height = 6)
