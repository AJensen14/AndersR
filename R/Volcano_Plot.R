#' Produce a volcano plot function
#'
#' @param daa_result Result from limma differential abundance analysis
#' @param pval_threshold threshold for the p-value (0.05)
#' @param logFC_threshold threshold for log2FC (1)
#'
#' @return A volcano plot
#' @export
#'
#' @examples
#' \dontrun{
#'   Volcano_Plot(result, 0.05, 1)
#' }
Volcano_Plot <- function(daa_result, pval_threshold = 0.05, logFC_threshold = 1) {

  volc_data <- daa_result %>%
    rownames_to_column(var = "X") %>%
    mutate(
      negLogP = -log10(P.Value),
      significant = (P.Value < pval_threshold) & (abs(logFC) >= logFC_threshold),
      direction = case_when(
        logFC > logFC_threshold & significant ~ "Up",
        logFC < -logFC_threshold & significant ~ "Down",
        TRUE ~ "ns"
      ),
      Regulation_adj = case_when(
        logFC > logFC_threshold & adj.P.Val < pval_threshold ~ "Up",
        logFC < -logFC_threshold & adj.P.Val < pval_threshold ~ "Down",
        TRUE ~ "ns"
      )) %>%
    mutate(Color = case_when(
      direction == "Up"~"#DC143C",
      direction == "Down"~"#1E90FF",
      TRUE ~ "grey"
    )) %>%
    mutate(label = case_when(direction == "Up" | direction == "Down"~X)) %>%
    mutate(label_adj = case_when(Regulation_adj == "Up" | Regulation_adj == "Down"~X)) %>%
    mutate(adj_sig = case_when(adj.P.Val <= 0.05~T,
                               adj.P.Val > 0.05~F)) %>%
    mutate(Color_label = lighten(Color, 0.5))

  plot_static <- ggplot(volc_data, aes(x = logFC, y = negLogP)) +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
    geom_point(aes(fill = Color), shape = 21, color = "black", alpha = 0.9) +
    theme_minimal() +
    scale_fill_identity() +
    labs(
      title = "Volcano Plot",
      x = "Fold change (log2)",
      y = "-log10(P-value)",
      fill = "Direction"
    ) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.tag = element_text(size = 18, face = "bold")
    ) +
    geom_label_repel(aes(label = label,
                         fontface = ifelse(adj_sig, "bold", "plain"),
                         fill = Color_label),
                     alpha = ifelse(volc_data$adj_sig, 1, 0.9),
                     min.segment.length = 0,
                     max.overlaps = 10)

  return(plot_static)
}
