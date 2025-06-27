#' Generate a dummy dataset
#'
#' @param n_proteins Number of proteins
#' @param group1 Name of group 1
#' @param group2 Name of group 2
#' @param n_replicates Number of desired replicates
#' @param missing_frac Percentage of missing values
#'
#' @return A dummy dataset
#' @export
#' @examples
#' \dontrun{
#'   Dummy_Data <- (n_proteins = 500,group1 = "Control",group2 = "Treated",n_replicates = 5,missing_frac = 0.2
#' }
#'
Dummy_Data <- function(
    n_proteins = 500,
    group1 = "Control",
    group2 = "Treated",
    n_replicates = 5,
    missing_frac = 0.2
) {

  set.seed(123)  # for reproducibility

  # Protein names
  protein_names <- paste0("Protein_", seq_len(n_proteins))

  # Baseline log-means (so after exponentiation, get different scales)
  log_means <- rnorm(n_proteins, mean = 20, sd = 1)

  # Choose ~10% proteins to be differentially abundant
  n_diff <- floor(n_proteins / 10)
  diff_indices <- sample(n_proteins, n_diff)

  # Simulate replicates for each group
  simulate_group <- function(is_treated = FALSE) {
    replicate(n_replicates, {
      # For each protein, simulate from log-normal
      log_mu <- log_means
      if (is_treated) {
        log_mu[diff_indices] <- log_mu[diff_indices] + log(4)  # ~4-fold change
      }
      rlnorm(n_proteins, meanlog = log_mu, sdlog = 0.3)
    })
  }

  control_samples <- simulate_group(FALSE)
  treated_samples <- simulate_group(TRUE)

  # Combine
  expr_matrix <- cbind(control_samples, treated_samples)
  colnames(expr_matrix) <- c(
    paste0(group1, "_", seq_len(n_replicates)),
    paste0(group2, "_", seq_len(n_replicates))
  )
  rownames(expr_matrix) <- protein_names

  # Add missing values
  if (missing_frac > 0) {
    n_values <- length(expr_matrix)
    n_missing <- floor(n_values * missing_frac)
    missing_idx <- arrayInd(sample(n_values, n_missing), .dim = dim(expr_matrix))
    expr_matrix[missing_idx] <- NA
  }

  dummy_data <- as.data.frame(expr_matrix)
  return(dummy_data)
}
