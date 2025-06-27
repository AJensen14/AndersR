#' Function to produce boxplots for proteins of interest
#'
#' @param data The normalised data
#' @param DA_result Results from DA analysis
#' @param matrix Metadata
#' @param proteins vector of proteins to produce boxplots
#' @param group1 First group in metadata (control)
#' @param group2 Second group in metadata (treated)
#' @param p_type The statistic to show ("raw" or "adj")
#'
#' @return Boxplot and the long data used to make the boxplot for more customisation
#' @export
#'
#' @examples
#' \dontrun{
#' Boxplot_Proteins(df, DA_result, matrix, c("AQP4", "ACTN3", "LRP1", "PRKA2"), "control", "treated", "adj")
#' }
Boxplot_Proteins <- function(data, DA_result, matrix, proteins, group1, group2, p_type = "raw") {
  # Subset expression data
  df <- as.data.frame(data) %>%
    rownames_to_column(var = "Protein") %>%
    filter(Protein %in% proteins)

  # Melt into long format
  df_long <- melt(df, id.vars = "Protein", variable.name = "Sample", value.name = "Expression")

  # Join group info
  df_long <- df_long %>%
    left_join(matrix, by = c("Sample" = "sample")) %>%
    mutate(group = factor(group, levels = c(group1, group2)))

  # Check for any unmatched samples
  if (any(is.na(df_long$group))) {
    warning("Some samples in data did not match the 'matrix' sample names.")
  }

  if (p_type == "raw"){

    pval_table <- DA_result %>%
      rownames_to_column(var = "Protein") %>%
      filter(Protein %in% proteins) %>%
      dplyr::select(Protein, P.Value) %>%
      mutate(label = paste0("p = ", signif(P.Value, 3)))

  } else if (p_type == "adj") {

    pval_table <- DA_result %>%
      rownames_to_column(var = "Protein") %>%
      filter(Protein %in% proteins) %>%
      dplyr::select(Protein, adj.P.Val) %>%
      mutate(label = paste0("p = ", signif(adj.P.Val, 3)))

  }

  # sort labels
  group_1_lab <- str_to_title(tolower(group1))
  group_2_lab <- str_to_title(tolower(group2))


  # Base plot
  p <- ggplot(df_long, aes(x = group, y = Expression, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 1) +
    facet_wrap(~ Protein, scales = "free_y") +
    labs(
      title = "Protein Expression by Group",
      x = "Group",
      y = "Expression"
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#1E90FF", "#DC143C"))+
    scale_x_discrete(labels = c(group_1_lab, group_2_lab))

  # Annotate p-values per facet
  p <- p + geom_text(
    data = pval_table,
    aes(
      x = 1.5,  # Middle between the groups
      y = Inf,  # Top of the panel
      label = label
    ),
    vjust = 1.2,
    hjust = 0.5,
    size = 3,
    inherit.aes = FALSE
  )

  print(p)

  return(df_long)
}
