#' Make a sample matrix
#'
#' @param data Dataset with sample information
#'
#' @return A sample matrix for downstream analysis
#' @export
#'
#' @examples
#' \dontrun{
#' Make_Matrix(data = df)
#' }
Make_Matrix <- function(data){

  sample <- colnames(data)
  group <- vector("character", length(sample))

  cat("Add in the grouping for each sample")

  for (i in seq_along(sample)) {
    prompt <- paste0("[", i, "] ", sample[i], ": ")
    input <- readline(prompt)
    group[i] <- ifelse(nzchar(input), input, sample[i])
  }

  matrix <- data.frame(sample, group)
  return(matrix)

}

