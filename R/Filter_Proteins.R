#' Filtering protein function
#'
#' Function to specify how you wish to filter out proteins.
#'
#' @param data The dataset to filter
#' @param matrix The metadata matrix
#' @param min_frac The percent of proteins to filter out (0 - 1)
#'
#' @return A filtered dataset
#' @export
#'
#' @examples
#' \dontrun{
#'   filtered_df <- Filter_Proteins(df, matrix, min_frac = 0.5)
#' }
Filter_Proteins <- function(data, matrix, min_frac = 0.7) {
  # Ensure order matches and use only matrix$sample columns
  df <- data[, matrix$sample, drop = FALSE]

  # For each group, find samples
  group_list <- split(matrix$sample, matrix$group)

  # For each protein, calculate % present in each group
  is_expressed_enough <- apply(df, 1, function(protein_row) {
    sapply(group_list, function(samples_in_group) {
      subset <- protein_row[samples_in_group]
      # Define "present" as non-NA and > 0
      mean(!is.na(subset) & subset > 0) >= min_frac
    })
  })

  # Keep proteins that are present in enough samples in *all* groups
  keep <- apply(is_expressed_enough, 2, all)

  # Filter rows
  filtered_df <- df[keep, , drop = FALSE]

  # Retain rownames properly
  return(filtered_df)
}


