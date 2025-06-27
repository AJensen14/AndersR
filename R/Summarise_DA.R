#' Summarise DA analysis function
#'
#'  Takes the results of differential abundance analysis and extracts significant proteins.
#' @param result The DA output from limma
#' @param logFC_threshold The logFC threshold (1)
#' @param p_value_threshold The pvalue threshold (0.05)
#'
#' @return A list with differentially abundant proteins
#' @export
#'
#' @examples
#' \dontrun{
#' Summary_result <- Summarise_DA(DA_results, 1, 0.05)
#' }
Summarise_DA <- function(result, logFC_threshold = 1, p_value_threshold = 0.05){

  up_reg_p <- result %>%
    rownames_to_column(var = "Protein") %>%
    filter(logFC > logFC_threshold & P.Value <= p_value_threshold) %>%
    pull(Protein) %>%
    unique()

  down_reg_p <- result %>%
    rownames_to_column(var = "Protein") %>%
    filter(logFC < -logFC_threshold & P.Value <= p_value_threshold) %>%
    pull(Protein) %>%
    unique()

  up_reg_p_adj <- result %>%
    rownames_to_column(var = "Protein") %>%
    filter(logFC > logFC_threshold & adj.P.Val <= p_value_threshold) %>%
    pull(Protein) %>%
    unique()

  down_reg_p_adj <- result %>%
    rownames_to_column(var = "Protein") %>%
    filter(logFC < -logFC_threshold & adj.P.Val <= p_value_threshold) %>%
    pull(Protein) %>%
    unique()

  print(paste("In total there were",
              (length(up_reg_p)+length(down_reg_p)),
               "significant proteins (p < 0.05)",
              "Following p-value adjustment there are",
              (length(up_reg_p_adj)+length(down_reg_p_adj)),
              "significant proteins."))

  print("Proteins can be accessed inside the list that is produced with this function")

  analysis_results <- list(up_p = up_reg_p,
                           down_p = down_reg_p,
                           up_p_adj = up_reg_p_adj,
                           down_p_adj = down_reg_p_adj)

  return(analysis_results)
}

