# Have a folder on desktop with three analysis types in
setwd("C:/Users/hlajens2/Desktop/")

df_report <- read_tsv(file = "MS Data/KA_250618_Equine0411_Cartilage_Report.tsv") %>%
  mutate(Naming = case_when(is.na(PG.Genes)~PG.UniProtIds,
                            !is.na(PG.Genes)~PG.Genes),
         Naming = str_remove(Naming, "^;"),
         Naming = unlist(lapply(strsplit(Naming, ";"), function(x)x[1])),
         Naming = make.unique(Naming))

# tidy up data
df_report_samples <- read_tsv(file = "MS Data/KA_250618_Equine0411_Cartilage_Report.tsv") %>%
  select(contains("Quantity")) %>%
  select(!contains("Log"))

# Change NaN to NA
df_report_samples[df_report_samples == "NaN"] <- NA

# Add in the row names
rownames(df_report_samples) <- df_report$Naming

# Now need a way to sort out the column names

df <- Tidy_Spectronaut(file_path = "C:/Users/hlajens2/Desktop/MS Data/KA_250618_Equine0411_Cartilage_Report.tsv")

# Testing column name function
test_vector <- letters[1:10]
new_names <- vector("character", length(test_vector))
for (i in seq_along(test_vector)) {
  prompt <- paste0("[", i, "] ", test_vector[i], ": ")
  input <- readline(prompt)
  new_names[i] <- ifelse(nzchar(input), input, test_vector[i])
}

# check data function
# list --> nprotein per sample, missing values (total)(same graph)
#


# Makes plot if not split by group


# Makes the plot if split by group
plot_group <- df %>%
  melt(variable.name = "sampleID") %>%
  mutate(group = case_when(str_detect(sampleID, "young")~"young",
                           str_detect(sampleID, "old")~"old")) %>%
  group_by(group,sampleID) %>%
  summarise(missing = sum(is.na(value)),
            detected = sum(!is.na(value))) %>%
  melt(id.vars = c('group', 'sampleID')) %>%
  ggplot(aes(x = sampleID, y = value)) +
  geom_bar(stat = "identity", col = "black",
           aes(fill = variable)) +
  labs(x = "Sample ID",
       y = "Number of proteins",
       fill = "Identified",
       title = "Number of proteins missing and detected across samples") +
  facet_wrap(~group, scales = "free_x",labeller = labeller(group = c("young" = str_to_title("young"),
                                                                     "old" = str_to_title("old")))) +
  scale_fill_manual(values = c("#F08080", "#48D1CC"),
                    labels = c("Missing", "Detected")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,
                                   size = 11,
                                   hjust = 1,
                                   vjust = 1.2,
                                   face = "bold"),
        axis.text.y = element_text(size = 11,
                                   face = "bold"),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# work out percent missing total
summary_plot <- df %>%
  melt() %>%
  summarise(Missing = sum(is.na(value)),
            Detected = sum(!is.na(value))) %>%
  mutate(Percent_Missing = Missing/(Missing+Detected)*100,
         Percent_Missing = round(Percent_Missing, 2)) %>%
  melt() %>%
  mutate(variable = factor(variable,
                           levels = c("Detected", "Missing", "Percent_Missing")),
         Position = c(1,1,1),
         value = as.character(value),
         value = case_when(variable == "Percent_Missing"~paste0(value, "%"),
                              variable != "Percent_Missing"~value)) %>%
  ggplot(aes(x = variable, y = Position))+
  geom_point(aes(fill = variable),
             col = "black",
             shape = 21,
             size = 50,
             alpha= 0.3,
             show.legend = F,
             stroke = 3)+
  labs(title = "Summary of protein detection across all samples")+
  geom_text(aes(label = value), size = 10)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 15,
                                   face = "bold",
                                   vjust = 70,
                                   colour = "black"),
        plot.title = element_text(size = 20, hjust = 0.5))+
  scale_fill_manual(values = c("#1E90FF", "#DC143C", "#8A2BE2"))

cowplot::plot_grid(plotlist = list(Missing = plot_group,
                                   Summary = summary_plot),
                   labels = c("A", "B"))


library(dplyr)
library(reshape2)  # or tidyr::pivot_longer
library(stringr)
library(ggplot2)
library(cowplot)

plot_missing_detected <- function(df,
                                  sample_id_var = "sampleID",
                                  group_patterns = c("young", "old", "middle"),
                                  group_labels = NULL,
                                  fill_colors = c("#F08080", "#48D1CC"),
                                  summary_colors = c("#1E90FF", "#DC143C", "#8A2BE2")) {

  # Melt the data
  df_long <- df %>%
    reshape2::melt(variable.name = sample_id_var)

  # Create group_labels if NULL
  if (is.null(group_labels)) {
    group_labels <- setNames(str_to_title(group_patterns), group_patterns)
  }

  # Helper function to assign groups based on patterns
  assign_group <- function(x, patterns) {
    assigned <- rep(NA_character_, length(x))
    for (pattern in patterns) {
      assigned[str_detect(x, pattern) & is.na(assigned)] <- pattern
    }
    assigned
  }

  # Assign groups dynamically
  df_long <- df_long %>%
    mutate(group = assign_group(!!rlang::sym(sample_id_var), group_patterns)) %>%
    filter(!is.na(group))

  # Summary data
  summary_data <- df_long %>%
    summarise(Missing = sum(is.na(value)),
              Detected = sum(!is.na(value))) %>%
    mutate(Percent_Missing = round(Missing/(Missing + Detected) * 100, 2)) %>%
    reshape2::melt() %>%
    mutate(variable = factor(variable, levels = c("Detected", "Missing", "Percent_Missing")),
           Position = 1,
           value = ifelse(variable == "Percent_Missing", paste0(value, "%"), as.character(value)))

  # Plot group
  plot_group <- df_long %>%
    group_by(group, !!rlang::sym(sample_id_var)) %>%
    summarise(missing = sum(is.na(value)),
              detected = sum(!is.na(value)), .groups = "drop") %>%
    reshape2::melt(id.vars = c("group", sample_id_var)) %>%
    ggplot(aes_string(x = sample_id_var, y = "value")) +
    geom_bar(stat = "identity", aes(fill = variable), col = "black") +
    labs(x = "Sample ID",
         y = "Number of proteins",
         fill = "Identified",
         title = "Number of proteins missing and detected across samples") +
    facet_wrap(~group, scales = "free_x", labeller = labeller(group = group_labels)) +
    scale_fill_manual(values = fill_colors, labels = c("Missing", "Detected")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, size = 11, hjust = 1, vjust = 1.2, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          legend.position = "top",
          legend.title = element_text(face = "bold"),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())

  # Summary plot
  summary_plot <- summary_data %>%
    ggplot(aes(x = variable, y = Position)) +
    geom_point(aes(fill = variable), col = "black", shape = 21, size = 50, alpha = 0.3, show.legend = FALSE, stroke = 3) +
    geom_text(aes(label = value), size = 10) +
    labs(title = "Summary of protein detection across all samples") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size = 15, face = "bold", colour = "black"),
          plot.title = element_text(size = 20, hjust = 0.5)) +
    scale_fill_manual(values = summary_colors)

  final_plot <- cowplot::plot_grid(plot_group, summary_plot, labels = c("A", "B"))
  return(final_plot)
}

plot_missing_detected(df,
                      group_patterns = c("young", "old", "middle"),
                      group_labels = c("young" = "Young",
                                       "old" = "Old"),
                      fill_colors = c("#F08080", "#48D1CC"))




# Now i need to do the normalisation testing function
# 1) Variation across replicates for each method (bar plot with stats comparing across normalisation methods)
# 2) density plots per condition for each method
# 3) PCA for each norm method with the grouping variables coloured
# 4) MA plot
# 5) RLE plot

median_norm <- medianNormalization(as.matrix(filtered_df)) %>%
  as.data.frame()
row.names(median_norm) <- rownames(filtered_df)

quantile_norm <- performQuantileNormalization(as.matrix(filtered_df)) %>%
  as.data.frame()
rownames(quantile_norm) <- rownames(filtered_df)

log_norm <- log2(filtered_df)
rownames(log_norm) <- rownames(filtered_df)

vsn_norm <- performVSNNormalization(as.matrix(filtered_df)) %>%
  as.data.frame()
rownames(vsn_norm) <- rownames(filtered_df)

# Now need to link each of these to the four types of analysis
calculate_cv <- function(data, matrix) {
  data_long <- data %>%
    rownames_to_column("Protein") %>%
    pivot_longer(-Protein, names_to = "sample", values_to = "Intensity") %>%
    left_join(matrix, by = "sample")

  cv_df <- data_long %>%
    group_by(Protein, group) %>%
    summarise(
      mean_int = mean(Intensity, na.rm = TRUE),
      sd_int = sd(Intensity, na.rm = TRUE),
      cv = sd_int / mean_int,
      .groups = "drop"
    )

  return(cv_df)
}

cv_median <- calculate_cv(median_norm, matrix) %>% mutate(Method = "Median")
cv_quantile <- calculate_cv(quantile_norm, matrix) %>% mutate(Method = "Quantile")
cv_log <- calculate_cv(log_norm, matrix) %>% mutate(Method = "Log2")
cv_vsn <- calculate_cv(vsn_norm, matrix) %>% mutate(Method = "VSN")

cv_all <- bind_rows(cv_median, cv_quantile, cv_log, cv_vsn)

# Plot
ggplot(cv_all, aes(x = Method, y = cv, fill = Method)) +
  geom_violin(aes(fill = group)) +
  theme_minimal() +
  labs(title = "Coefficient of variation across replicates", y = "Coefficient of Variation")

make_density <- function(data){

  data_long <- data %>%
    rownames_to_column("Protein") %>%
    pivot_longer(-Protein, names_to = "sample", values_to = "Intensity") %>%
    left_join(matrix, by = "sample")

}

dp_med <- make_density(median_norm) %>% mutate(Method = "Median")
dp_log <- make_density(log_norm) %>% mutate(Method = "Log2")
dp_qua <- make_density(quantile_norm) %>% mutate(Method = "Quantile")
dp_vsn <- make_density(vsn_norm) %>% mutate(Method = "VSN")

dp_all <- bind_rows(dp_med, dp_log, dp_qua, dp_vsn)

# Plot
ggplot(dp_all, aes(x = Intensity)) +
  geom_density(aes(fill = group),
               alpha = 0.6) +
  theme_minimal() +
  facet_wrap(~Method)+
  labs(y = "Density", fill = "Condition")

# Now do PCA
plot_pca <- function(norm_df, method_name, sample_metadata) {
  # Transpose so PCA is done on samples
  norm_df[is.na(norm_df)] <- 0
  pca <- prcomp(t(norm_df), scale. = T)
  pca_df <- as.data.frame(pca$x)
  pca_df$sample <- rownames(pca_df)

  # Add group info
  pca_df <- left_join(pca_df, sample_metadata, by = "sample")

  # Plot
  ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = sample)) +
    geom_point(size = 3) +
    labs(
      title = paste("PCA -", method_name),
      x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
      y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)")
    ) +
    theme_minimal()+
    geom_label_repel(aes(label = sample))
}

pca_plots <- list(
  Median   = plot_pca(median_norm, "Median", matrix),
  Quantile = plot_pca(quantile_norm, "Quantile", matrix),
  Log2     = plot_pca(log_norm, "Log2", matrix),
  VSN      = plot_pca(vsn_norm, "VSN", matrix)
)

# To view them one at a time:
plot_grid(plotlist = pca_plots)


# MA plot
plot_ma <- function(norm_df, method_name, metadata) {
  # Get two samples from different groups
  group_samples <- split(metadata$sample, metadata$group)
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
    labs(title = paste("MA Plot -", method_name), x = "A (mean log2 intensity)", y = "M (log2 fold change)") +
    theme_minimal()
}

plot_rle <- function(norm_df, method_name, metadata) {
  medians <- apply(norm_df, 1, median, na.rm = TRUE)
  rle_mat <- sweep(norm_df, 1, medians)  # subtract row medians

  df_rle <- as.data.frame(rle_mat)
  df_rle$protein <- rownames(df_rle)
  df_long <- tidyr::pivot_longer(df_rle, -protein, names_to = "sample", values_to = "rle")
  df_long <- left_join(df_long, metadata, by = "sample")

  ggplot(df_long, aes(x = sample, y = rle, fill = group)) +
    geom_boxplot(outlier.size = 0.2) +
    labs(title = paste("RLE Plot -", method_name), y = "Relative Log Expression") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

ma_plots <- list(
  Median = plot_ma(median_norm, "Median", matrix),
  Quantile = plot_ma(quantile_norm, "Quantile", matrix),
  Log2 = plot_ma(log_norm, "Log2", matrix),
  VSN = plot_ma(vsn_norm, "VSN", matrix)
)

# RLE plots
rle_plots <- list(
  Median = plot_rle(median_norm, "Median", matrix),
  Quantile = plot_rle(quantile_norm, "Quantile", matrix),
  Log2 = plot_rle(log_norm, "Log2", matrix),
  VSN = plot_rle(vsn_norm, "VSN", matrix)
)

cowplot::plot_grid(plotlist = ma_plots)
cowplot::plot_grid(plotlist = rle_plots)

# Testing how to do DA
group_factor <- factor(matrix$group)
design <- model.matrix(~0 + group_factor)
colnames(design) <- levels(group_factor)

# Fit linear model
fit <- lmFit(normalised_df, design)

# Define contrasts (example: compare group2 vs group1)
contrast.matrix <- makeContrasts(
  old_vs_young = old - young,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get top results
result <- topTable(fit2, coef = "old_vs_young", number = Inf, adjust.method = "BH")










