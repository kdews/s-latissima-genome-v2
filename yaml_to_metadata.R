## Initialization
# Clear environment
rm(list = ls())
# Load required packages
library(yaml)

# Functions
# Order for label factor
lbl_order <- c(
  "lindell",
  "corteva",
  "roscoff",
  "ilvo",
  "japonica",
  "ectocarpus"
)
# Return error per YAML file
metadata_error <- function(file, field, message) {
  stop(file, ": ", field, ": ", message, call. = F)
}
# Return a nested YAML value or report which required field is missing
get_field <- function(record, file, ...) {
  keys <- unlist(list(...))
  value <- record
  traversed <- character()
  
  for (key in keys) {
    traversed <- c(traversed, key)
    if (!is.list(value) || is.null(value[[key]])) {
      metadata_error(file, paste(traversed, collapse = "."),
                     "required field is missing")
    }
    value <- value[[key]]
  }
  return(value)
}
# Convert YAML scalar to TSV text
# YAML null becomes an empty TSV field when required = F
scalar_text <- function(value, file, field, required = T) {
  if (is.null(value)) {
    if (required) metadata_error(file, field, "cannot be null")
    return("")
  }
  if (is.list(value) || length(value) != 1) {
    metadata_error(file, field, "must be scalar value")
  }
  text <- trimws(as.character(value))
  if (required && !nzchar(text)) metadata_error(file, field, "cannot be blank")
  if (grepl("[\t\r\n]", text)) {
    metadata_error(file, field, "cannot contain tabs or newlines")
  }
  return(text)
}
# Convert YAML list fields (e.g., [acc1, acc2]) to comma-separated strings for TSV
locator_text <- function(value, file, field) {
  if (is.list(value) || length(value) > 1) {
    values <- unlist(value, use.names = FALSE)
    if (length(values) == 0) metadata_error(file, field, "list cannot be empty")
    values <- vapply(values, scalar_text, character(1),
                     file = file, field = field)
    if (anyDuplicated(values)) {
      metadata_error(file, field, "contains duplicate values")
    }
    if (any(grepl(",", values, fixed = T))) {
      metadata_error(file, field, "list entries cannot contain commas")
    }
    return(paste(values, collapse = ","))
  }
  return(scalar_text(value, file, field))
}
# Parse YAML metadata file
parse_yaml <- function(file) {
  cat("Reading:", file, "\n")
  record <- read_yaml(file)
  # Check YAML is dictionary list
  if (!is.list(record)) stop(file, ": YAML root must be a mapping")
  # Parse each data frame column from YAML
  label <- scalar_text(get_field(record, file, "label"), file, "label")
  pretty_label <- scalar_text(get_field(record, file, "pretty_label"),
                              file, "pretty_label")
  genome_location <- tolower(scalar_text(get_field(record, file, "assembly",
                                                   "source", "location"),
                                         file, "assembly.source.location"))
  genome_method <- tolower(scalar_text(get_field(record, file, "assembly",
                                                 "source", "method"),
                                       file, "assembly.source.method"))
  genome_locator <- locator_text(get_field(record, file, "assembly", "source",
                                           "locator"),
                                 file, "assembly.source.locator")
  # selector is a required YAML key, but its value may be null
  assembly_source <- get_field(record, file, "assembly", "source")
  if (!"selector" %in% names(assembly_source)) {
    metadata_error(file, "assembly.source.selector", "required field missing")
  }
  genome_selector <- scalar_text(assembly_source$selector, file,
                                 "assembly.source.selector", required = F)
  reads_location <- tolower(scalar_text(get_field(record, file, "reads",
                                                  "source", "location"),
                                        file, "reads.source.location"))
  reads_method <- tolower(scalar_text(get_field(record, file, "reads", "source",
                                                "method"),
                                      file, "reads.source.method"))
  reads_locator <- locator_text(get_field(record, file, "reads", "source",
                                          "locator"),
                                file, "reads.source.locator")
  reads_type <- tolower(scalar_text(get_field(record, file, "reads", "platform"),
                                    file, "reads.platform"))
  if (reads_method == "sra") {
    accessions <- trimws(strsplit(reads_locator, ",", fixed = T)[[1]])
    invalid <- accessions[!grepl("^[SED]RR[0-9]+$", accessions)]
    if (length(invalid) > 0) {
      metadata_error(file, "reads.source.locator",
                     paste0("invalid SRA run accession(s): ",
                            paste(invalid, collapse = ", ")))
    }
  }
  # Combine fields into data frame
  df <- data.frame(
    Label = label,
    Pretty_Label = pretty_label,
    Genome_Source = genome_location,
    Genome_Method = genome_method,
    Genome_Locator = genome_locator,
    Genome_Selector = genome_selector,
    Reads_Source = reads_location,
    Reads_Method = reads_method,
    Reads_Type = reads_type,
    Reads_Locator = reads_locator,
    check.names = F,
    stringsAsFactors = F
  )
  return(df)
}

# Input
# Only take command line input if not running interactively
if (interactive()) {
  # Set working dir
  wd <- "/scratch1/kdeweese/corteva_genome"
  setwd(wd)
  yaml_dir <- "s-latissima-genome-v2/assemblies_metadata"
  out_tsv <- "s-latissima-genome-v2/genomes_metadata.tsv"
} else {
  line_args <- commandArgs(trailingOnly = T)
  cat("Arguments:", "\n")
  cat(line_args, "\n", sep = "\n")
  yaml_dir <- line_args[1]
  out_tsv <- line_args[2]
}

# Check input directory exists
if (!dir.exists(yaml_dir)) stop("YAML directory does not exist: ", yaml_dir)
# List YAML files
yaml_files <- list.files(yaml_dir, pattern = "\\.ya?ml$",
                         full.names = T, ignore.case = T)

# Read all YAML files into combined data frame
metadata <- do.call(rbind, lapply(yaml_files, parse_yaml))
# Sort "Label" column with factor
lbl_lvls <- sapply(lbl_order, grep, sort(metadata$Label), value = T,
                   perl = T) |>
  unname() |>
  unlist()
metadata <- metadata[order(factor(metadata$Label, levels = lbl_lvls)), ]

# Error if duplicated labels are found
if (anyDuplicated(metadata$Label)) {
  duplicates <- unique(metadata$Label[duplicated(metadata$Label)])
  stop("Duplicate label(s): ", paste(duplicates, collapse = ", "), call. = F)
}

# Write output metadata TSV
colnames(metadata)[1] <- paste("#", colnames(metadata)[1])
write.table(metadata, file = out_tsv, row.names = F, quote = F, na = "",
            sep = "\t")

# Log details
cat("Parsed", nrow(metadata), "YAML files into",
    ncol(metadata), "columns and", nrow(metadata), "records.\n")
cat("Output written to:", out_tsv, "\n")
