#' Function to see the effects of filtering threshold
#'
#' @param data The non-filtered data
#' @param matrix The meta data
#' @param n_samples The number of samples per group
#'
#' @return A plot showing filtering effect
#' @export
#'
#' @examples
#' \dontrun{
#' Check_Filtering(df, matrix, 10)
#' }
Check_Filtering <- function(data, matrix, n_samples) {

  # Ensure order matches
  df <- data[, matrix$sample, drop = FALSE]

  # For each group, find samples
  group_list <- split(matrix$sample, matrix$group)

  # Prepare a vector of thresholds (from 3 to n_samples)
  thresholds <- 3:n_samples

  # Store counts
  protein_counts <- numeric(length(thresholds))

  # For each threshold, calculate how many proteins would be retained
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]

    # For each protein, check in each group how many replicates have expression
    is_expressed_enough <- apply(df, 1, function(protein_row) {
      sapply(group_list, function(samples_in_group) {
        subset <- protein_row[samples_in_group]
        sum(!is.na(subset) & subset > 0) >= threshold
      })
    })

    # Keep proteins present in enough samples in *all groups*
    keep <- apply(is_expressed_enough, 2, all)

    # Store how many proteins would remain at this threshold
    protein_counts[i] <- sum(keep)
  }

  # Create a data frame for plotting
  plot_df <- data.frame(
    Min_Replicates = thresholds,
    Num_Proteins = protein_counts
  )

  # Plot
  p <- ggplot(plot_df, aes(x = Min_Replicates, y = Num_Proteins)) +
    geom_line(size = 1, color = "steelblue") +
    geom_point(size = 2, color = "steelblue") +
    labs(
      title = "Effect of Filtering Threshold on Protein Count",
      x = "Minimum Replicates per Group",
      y = "Number of Proteins Retained"
    ) +
    theme_minimal()+
    scale_x_continuous(breaks = seq(3,n_samples, 1))+
    geom_text(aes(label = Num_Proteins), nudge_y = max((plot_df$Num_Proteins)/20))

  return(p)

}

