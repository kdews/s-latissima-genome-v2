# Annotate sample data using sample IDs and a population-code table.
#
# df:        Data frame containing sample identifiers.
# label:     Name of the sample-ID column, supplied as a character string.
# pop_codes: Population metadata data frame containing `Abbreviation`.

metaDf <- function(df, label, pop_codes) {
  # Expand abbreviations
  sex_convert <- c("FG" = "female", "MG" = "male")
  subsp_convert <- c("SL" = "latissima", "SA" = "angustissima")
  # Correct aberrant IDs
  id_corrections = c(
    "LIS-F1-3" = "SL-SNE-1-FG-3",
    "CT1" = "CT1-0",
    "Female" = "FG-0"
  )
  sample_ids <- as.character(df[[label]])
  for (old in names(id_corrections)) {
    sample_ids <- gsub(old, id_corrections[[old]], sample_ids, fixed = T)
  }
  # Extract metadata from sample IDs
  df %>%
    dplyr::mutate(rename_id = sample_ids) %>%
    tidyr::separate_wider_delim(
      rename_id,
      delim = "-",
      names = c("subspecies", "pop", "n_samp", "sex", "n_gam")
    ) %>%
    dplyr::mutate(
      subspecies = unname(subsp_convert[subspecies]),
      sex = unname(sex_convert[sex])
    ) %>%
    dplyr::left_join(pop_codes, by = dplyr::join_by(pop == Abbreviation))
}
