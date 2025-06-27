#' Function to review missingness of data
#'
#' Summarising the amount of missing data within a proteomics dataset
#'
#' @param data The dataframe
#' @param group_patterns The group within your experiment
#' @param group_labels The labels for the figures for each groups
#' @param fill_colors Colours for plot A
#' @param summary_colors Colours for plot B
#'
#' @return A cowplot with two figures. A) Missing proteins per sample. B) Summary of missingness of data
#' @export
#'
#' @examples
#' \dontrun{
#' Check_Data(data = df,
#'            group_patterns = c("young", "old"))
#' }
Check_Data <- function(data,
                       group_patterns = c("young", "old"),
                       group_labels = NULL,
                       fill_colors = c("#F08080", "#48D1CC"),
                       summary_colors = c("#1E90FF", "#DC143C", "#8A2BE2")) {

  # Melt the data
  df_long <- data %>%
    reshape2::melt(variable.name = "sample_id_var")

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
    mutate(group = assign_group(!!rlang::sym("sample_id_var"), group_patterns)) %>%
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
    group_by(group, !!rlang::sym("sample_id_var")) %>%
    summarise(missing = sum(is.na(value)),
              detected = sum(!is.na(value)), .groups = "drop") %>%
    reshape2::melt(id.vars = c("group", "sample_id_var")) %>%
    ggplot(aes_string(x = "sample_id_var", y = "value")) +
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
