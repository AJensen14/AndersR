#' Normalisation comparison function
#'
#' @param filtered_df Filtered dataframe
#' @param sample_matrix Metadata for your experiment
#'
#' @return Coefficient of variation, Density, Principal component analysis, MA and Relative log expression plots
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Test_Normalisation(filtered_df, matrix)
#' }
Test_Normalisation <- function(filtered_df, sample_matrix) {

  # 1. Normalize data --------------------------------------------------------
  median_norm <- medianNormalization(as.matrix(filtered_df)) %>% as.data.frame()
  rownames(median_norm) <- rownames(filtered_df)

  quantile_norm <- performQuantileNormalization(as.matrix(filtered_df)) %>% as.data.frame()
  rownames(quantile_norm) <- rownames(filtered_df)

  log_norm <- log2(filtered_df)
  rownames(log_norm) <- rownames(filtered_df)

  vsn_norm <- performVSNNormalization(as.matrix(filtered_df)) %>% as.data.frame()
  rownames(vsn_norm) <- rownames(filtered_df)

  norm_list <- list(Median = median_norm, Quantile = quantile_norm, Log2 = log_norm, VSN = vsn_norm)

  # 2. CV calculation ---------------------------------------------------------
  calculate_cv <- function(data, meta) {
    data_long <- data %>%
      rownames_to_column("Protein") %>%
      pivot_longer(-Protein, names_to = "sample", values_to = "Intensity") %>%
      left_join(meta, by = "sample")

    data_long %>%
      group_by(Protein, group) %>%
      summarise(
        mean_int = mean(Intensity, na.rm = TRUE),
        sd_int = sd(Intensity, na.rm = TRUE),
        cv = sd_int / mean_int,
        .groups = "drop"
      )
  }

  cv_all <- bind_rows(
    lapply(names(norm_list), function(m) {
      calculate_cv(norm_list[[m]], sample_matrix) %>% mutate(Method = m)
    })
  )

  cv_plot <- ggplot(cv_all, aes(x = Method, y = cv, fill = group)) +
    geom_violin(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Coefficient of Variation Across Replicates", y = "Coefficient of Variation", fill = "Group")

  # 3. Density plots -----------------------------------------------------------
  make_density <- function(data, meta) {
    data %>%
      rownames_to_column("Protein") %>%
      pivot_longer(-Protein, names_to = "sample", values_to = "Intensity") %>%
      left_join(meta, by = "sample")
  }

  dp_all <- bind_rows(
    lapply(names(norm_list), function(m) {
      make_density(norm_list[[m]], sample_matrix) %>% mutate(Method = m)
    })
  )

  density_plot <- ggplot(dp_all, aes(x = Intensity, fill = group)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Method, scales = "free") +
    theme_minimal() +
    labs(y = "Density", fill = "Group", title = "Density Plots Per Normalization Method")

  # 4. PCA plots ---------------------------------------------------------------
  plot_pca <- function(norm_df, method_name, meta) {
    norm_df[is.na(norm_df)] <- 0
    pca <- prcomp(t(norm_df), scale. = TRUE)
    pca_df <- as.data.frame(pca$x)
    pca_df$sample <- rownames(pca_df)
    pca_df <- left_join(pca_df, meta, by = "sample")

    ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
      geom_point(size = 3) +
      labs(title = paste("PCA -", method_name),
           x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
           y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)"),
           color = "Group") +
      theme_minimal() +
      geom_label_repel(aes(label = sample))
  }

  pca_plots <- list(
    Median   = plot_pca(norm_list$Median, "Median", sample_matrix),
    Quantile = plot_pca(norm_list$Quantile, "Quantile", sample_matrix),
    Log2     = plot_pca(norm_list$Log2, "Log2", sample_matrix),
    VSN      = plot_pca(norm_list$VSN, "VSN", sample_matrix)
  )

  pca_grid <- cowplot::plot_grid(plotlist = pca_plots, ncol = 2)

  # 5. MA plot ---------------------------------------------------------------
  plot_ma <- function(norm_df, method_name, meta) {
    group_samples <- split(meta$sample, meta$group)
    if (length(group_samples) < 2) stop("Need at least two groups for MA plot")
    s1 <- group_samples[[1]][1]
    s2 <- group_samples[[2]][1]

    x <- norm_df[[s1]]
    y <- norm_df[[s2]]

    A <- (log2(x) + log2(y)) / 2
    M <- log2(x) - log2(y)

    df_ma <- data.frame(A = A, M = M)

    ggplot(df_ma, aes(x = A, y = M)) +
      geom_point(alpha = 0.4) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      labs(title = paste("MA Plot -", method_name),
           x = "A (mean log2 intensity)", y = "M (log2 fold change)") +
      theme_minimal()
  }

  ma_plots <- list(
    Median = plot_ma(norm_list$Median, "Median", sample_matrix),
    Quantile = plot_ma(norm_list$Quantile, "Quantile", sample_matrix),
    Log2 = plot_ma(norm_list$Log2, "Log2", sample_matrix),
    VSN = plot_ma(norm_list$VSN, "VSN", sample_matrix)
  )

  ma_grid <- cowplot::plot_grid(plotlist = ma_plots, ncol = 2)

  # 6. RLE plot ---------------------------------------------------------------
  plot_rle <- function(norm_df, method_name, meta) {
    medians <- apply(norm_df, 1, median, na.rm = TRUE)
    rle_mat <- sweep(norm_df, 1, medians)

    df_rle <- as.data.frame(rle_mat)
    df_rle$protein <- rownames(df_rle)
    df_long <- tidyr::pivot_longer(df_rle, -protein, names_to = "sample", values_to = "rle")
    df_long <- left_join(df_long, meta, by = "sample")

    ggplot(df_long, aes(x = sample, y = rle, fill = group)) +
      geom_boxplot(outlier.size = 0.2) +
      labs(title = paste("RLE Plot -", method_name),
           y = "Relative Log Expression",
           x = "Sample") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      guides(fill = guide_legend(title = "Group"))
  }

  rle_plots <- list(
    Median = plot_rle(norm_list$Median, "Median", sample_matrix),
    Quantile = plot_rle(norm_list$Quantile, "Quantile", sample_matrix),
    Log2 = plot_rle(norm_list$Log2, "Log2", sample_matrix),
    VSN = plot_rle(norm_list$VSN, "VSN", sample_matrix)
  )

  rle_grid <- cowplot::plot_grid(plotlist = rle_plots, ncol = 2)

  # Return everything in a named list ------------------------------------------
  results <- list(
    CV = cv_plot,
    Density = density_plot,
    PCA = pca_grid,
    MA = ma_grid,
    RLE = rle_grid)

  return(results)

}

