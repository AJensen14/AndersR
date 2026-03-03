#' Filtering protein function
#'
#' Function to specify how you wish to filter out proteins.
#'
#' @param data The dataset to filter
#' @param matrix The metadata matrix
#' @param min_samples The minimum number of samples the protein needs to be in.
#'
#' @return A filtered dataset
#' @export
#'
#' @examples
#' \dontrun{
#'   filtered_df <- Filter_Proteins(df, matrix, min_samples = 3)
#' }
#'
Filter_Proteins <- function(data, matrix, min_samples) {
  # Ensure order matches and use only matrix$sample columns
  df_check <- data[, matrix$sample, drop = FALSE] # makes sure you only get samples from your matrix
  rownames(df_check) <- rownames(data) # puts rownames back as they are lost

  # For each protein, calculate % present in each group
  long_data <- df_check %>%
    rownames_to_column(var = "protein") %>%
    melt() %>%
    dplyr::rename(sample = variable) %>%
    left_join(matrix, by = c("sample" = "sample")) %>%
    group_by(protein, group) %>%
    summarise(detected = count(!is.na(value))) %>%
    filter(detected >= min_samples) %>%
    ungroup() %>%
    group_by(protein) %>%
    filter(n() > 1) %>%
    ungroup()

  # make vector of good proteins
  proteins <- unique(long_data$protein)

  # filter df_check based on this
  filtered_df <- df_check %>%
    rownames_to_column(var = "protein") %>%
    filter(protein %in% proteins) %>%
    column_to_rownames(var = "protein")

  # Retain rownames properly
  return(filtered_df)
}
