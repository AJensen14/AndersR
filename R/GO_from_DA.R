#' Function to input the results of differential expression into ClusterProfiler.
#'
#' Use SummariseDA function before this
#'
#' @param DA_summary The result of SummariseDA
#' @param ontology The type of enrichment analysis MF, CC or BP
#' @param use_adj_p Use proteins significant at p_adj (True) or raw_p (FALSE)
#' @param group_up The name of the group associated with up-regulation
#' @param group_down The name of the group associated with down-regulation
#' @param top_n The number of terms to show on the dotplot
#'
#' @return A dotplot from clusterprofiler
#' @export
#'
#' @examples
#' \dontrun{
#' GO_from_DA(DA_summary, "BP", TRUE, "Old", "Young", 10)
#' }
GO_from_DA <- function(DA_summary,
                       ontology = "BP",
                       use_adj_p = TRUE,
                       group_up = "Up",
                       group_down = "Down",
                       top_n = 10) {
  # Validate ontology
  ontology <- toupper(ontology)
  if (!ontology %in% c("BP", "MF", "CC")) {
    stop("Ontology must be one of: 'BP', 'MF', or 'CC'.")
  }

  # Select gene lists based on p-value or adj p-value
  up_genes <- if (use_adj_p) DA_summary$up_p_adj else DA_summary$up_p
  down_genes <- if (use_adj_p) DA_summary$down_p_adj else DA_summary$down_p

  # Run enrichment for upregulated genes
  ego_up <- enrichGO(
    gene          = up_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ontology,
    pAdjustMethod = "BH",
    readable      = TRUE
  )

  # Run enrichment for downregulated genes
  ego_down <- enrichGO(
    gene          = down_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ontology,
    pAdjustMethod = "BH",
    readable      = TRUE
  )

  # Create dotplots
  plot_up <- dotplot(ego_up, showCategory = top_n) +
    ggtitle(paste("More abundant in", group_up))

  plot_down <- dotplot(ego_down, showCategory = top_n) +
    ggtitle(paste("More abundant in", group_down))
  # Combine the two dotplots
  combined_plot <- cowplot::plot_grid(
    plot_up, plot_down,
    labels = c("A", "B"),
    ncol = 2,
    align = "h"
  )

  # Add a title to the whole plot
  title <- ggdraw() +
    draw_label(
      paste0("GO Enrichment (", ontology, ")"),
      fontface = 'bold',
      x = 0.5, hjust = 0.5
    )

  final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
  return(final_plot)
}


