#' Function to perform differential abundance analysis using the limma package
#'
#' @param normalised_data Normalised dataframe
#' @param sample_matrix Metadata of experiment
#' @param contrast_formula A charachter to specify how to do the DAA eg. ("disease - control")
#' @param adjust.method Method to adjust p-values ("None", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni")
#'
#' @return Differential expression table
#' @export
#'
#' @examples
#' \dontrun{
#' Differential_Abundance_Analysis(normalised_df,
#'                                 matrix,
#'                                 "old - young",
#'                                 adjust.method = "BH")
#' }
#'
Differential_Abundance_Analysis <- function(normalised_data, sample_matrix, contrast_formula = NULL, adjust.method = "BH") {
  # Convert group to factor
  group_factor <- factor(sample_matrix$group)

  # Design matrix
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- levels(group_factor)

  # Fit linear model
  fit <- limma::lmFit(normalised_data, design)

  # Define contrasts
  if (is.null(contrast_formula)) {
    # Default contrast: second level vs first level
    if (ncol(design) < 2) stop("Need at least two groups for contrast.")
    contrast_formula <- paste0(colnames(design)[2], "-", colnames(design)[1])
  }

  contrast.matrix <- limma::makeContrasts(
    contrasts = contrast_formula,
    levels = design
  )

  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)

  # Get top table for all features
  result <- limma::topTable(fit2, coef = 1, number = Inf, adjust.method = adjust.method)

  return(result)
}
