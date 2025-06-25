#' Reads and process protein assay from xlsx sheet
#'
#' XLSX file from spectrostar Nano always produces the same sheet and this
#' function attempts to analyse it quickly based of my processing protocol.
#'
#' This protocol hasn't let me down so far.
#'
#' @param file_path Where the file is saved
#' @param n_standard How many replicates of the standards
#' @param standard_cols The number the columns are in
#' @param standard_conc Vector with standard concentrations
#' @param n_samples How many samples in the 96-well plate
#' @param sample_order A vector of sample order. This assumes you went across in triplicates.
#' @param dodgy_standard Added this in cases when one of your standards is a bit rubbish
#'
#' @return A list with all the analysis outputs
#' @export
#'
#' @examples
#' Protein_Assay(file_path = "protein_assay_result.xlsx",
#' n_standard = 8,
#' standard_cols = c(1, 2, 3),
#' standard_conc = c(2000, 1500, 1000, 750, 500, 250, 125, 62.5),
#' n_samples = 24,
#' sample_order = LETTERS[1:24])
#'
Protein_Assay <- function(file_path,
                          n_standard = NULL,
                          standard_cols = NULL,
                          standard_conc = NULL,
                          n_samples = NULL,
                          sample_order = NULL,
                          dodgy_standard = NULL) {

  # make empty list to store things
  Analysis_folder <- list()

  # Default values if arguments are NULL
  if (is.null(n_standard)) n_standard <- 8
  if (is.null(standard_cols)) standard_cols <- c(1, 2, 3)
  if (is.null(standard_conc)) standard_conc <- c(2000, 1500, 1000, 750, 500, 250, 125, 62.5)
  if (is.null(n_samples)) n_samples <- 24
  if (is.null(sample_order)) sample_order <- LETTERS[1:20]
  if (is.null(dodgy_standard)) dodgy_standard <- c(1)

  # Now need code from above
  df <- read_xlsx(path = file_path)

  # filter out non-complete rows
  df_data <- df %>%
    filter(complete.cases(.)) %>%
    column_to_rownames(var = "User: USER")

  # sort out the column names
  colnames(df_data) <- c(1:30)[1:ncol(df_data)]

  # Now make all values 3 decimal places
  df_data_cleaned <- df_data %>%
    mutate(across(, as.numeric)) %>%
    mutate(across(where(is.numeric), ~ round(., 3)))

  # Make heat map of the plate
  heatmap <- df_data_cleaned %>%
    rownames_to_column(var = "row_letter") %>%
    melt() %>%
    ggplot(aes(x = variable, y = (row_letter)))+
    geom_tile(aes(fill = value))+
    scale_fill_gradient2(low = "#4575b4",
                         mid =  "#ffffbf",
                         high = "#d73027")+
    geom_text(aes(label = value))+
    scale_y_discrete(limits = rev)+
    scale_x_discrete(position = "top")+
    theme_minimal()+
    labs(x = "Column number",
         y = "Row letter",
         fill = "Absorbance")+
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 15))

  # Now basically need to take first three cols as standards
  standards <- df_data_cleaned %>%
    dplyr::select(all_of(standard_cols)) %>%
    mutate(avg_abs = rowMeans(across(where(is.numeric))),
           sd = apply(across(where(is.numeric)), 1, sd))

  standards <- standards[1:n_standard, ]

  # Now need to somehow add in the concentration of protein standards
  # easy when coding normally but more difficult for shiny app

  # add in standard data and plot
  final_data <- standards %>%
    mutate(conc = standard_conc) %>%
    filter(!conc %in% dodgy_standard)

  model <- lm(formula = avg_abs~conc, data = final_data) # same

  # extract coefficients
  gradient <- model$coefficients[2] # same
  intercept <- model$coefficients[1] # same

  # plot data
  standards_plot <- final_data %>%
    ggplot(aes(x = conc, y = avg_abs))+
    geom_point()+
    geom_smooth(method = 'lm')+
    geom_errorbar(aes(x = conc, ymin = avg_abs-sd, ymax = avg_abs+sd))+
    annotate(geom = 'text',
             x = max(final_data$conc)/3.2,
             y = max(final_data$avg_abs)/0.8,
             label = paste("Absorbance =",
                           round(gradient, digits = 5),
                           "* Concentration +",
                           round(intercept, digits = 4)))+
    annotate(geom = 'text',
             x = max(final_data$conc)/3.2,
             y = max(final_data$avg_abs)/0.9,
             label = paste("Concentration = Absorbance - ",
                           round(intercept, digits = 4),
                           "/",
                           round(gradient, digits = 5)))+
    theme_minimal()+
    labs(x = "Standard concentration (ug/mL)",
         y = "Average absorbance")+
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 15))+
    xlim(-300,2300)+
    ylim(0,(max(final_data$avg_abs)*1.5))

  # So now I need to average out triplicates from my original assay
  samples <- df_data_cleaned %>%
    dplyr::select(!all_of(standard_cols))

  # Then run a loop to save all the average and sd values into a long dataframe
  count_storage <- 1
  val_1 <- c()
  val_2 <- c()
  val_3 <- c()
  means <- c()
  sd_val <- c()

  for (row in 1:nrow(samples)) {

    for (col in c(1,4,7)) {

      # Retrieve trip values
      value_1 <- as.numeric(samples[row, col])
      value_2 <- as.numeric(samples[row, col+1])
      value_3 <- as.numeric(samples[row, col+2])

      # get summary results
      mean_value <- (value_1+value_2+value_3)/3
      sd_value <-  sd(c(value_1, value_2, value_3))

      # save information
      val_1[count_storage] <- value_1
      val_2[count_storage] <- value_2
      val_3[count_storage] <- value_3
      means[count_storage] <- mean_value
      sd_val[count_storage] <- sd_value

      count_storage <- count_storage+1
    }
  }

  # Store as dataframe
  df_results <- data.frame(val_1,
                           val_2,
                           val_3,
                           means,
                           sd_val)

  # So now we can state how many samples we have
  # in this example we have 20
  actual_samples <- df_results[1:n_samples, ]

  # Then you may want to specify sample order
  samples_results <- actual_samples %>%
    mutate(sample_order = sample_order) %>%
    mutate(Concentration = ((means - intercept)/gradient))

  # Now make a plot with standards and samples with ID labels
  samples_on_stand_plot <- final_data %>%
    ggplot(aes(x = conc, y = avg_abs))+
    geom_point()+
    geom_smooth(method = 'lm')+
    geom_errorbar(aes(x = conc, ymin = avg_abs-sd, ymax = avg_abs+sd))+
    geom_point(inherit.aes = F, data = samples_results, col = "red",
               aes(x = Concentration, y = means))+
    geom_label_repel(inherit.aes = F,
                     data = samples_results,
                     min.segment.length = 0,
                     nudge_y = 0.3,
                     max.overlaps = 20,
                     aes(x = Concentration,
                         y = means,
                         label = sample_order))+
    theme_minimal()+
    labs(x = "Standard concentration (ug/mL)",
         y = "Average absorbance")+
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 15))+
    xlim(-300,2300)+
    ylim(0,(max(final_data$avg_abs)*1.5))

  Analysis_folder$heatmap <- heatmap
  Analysis_folder$standards_plot <- standards_plot
  Analysis_folder$standards_result <- final_data
  Analysis_folder$samples_plot <- samples_on_stand_plot
  Analysis_folder$samples_result <- samples_results

  return(Analysis_folder)

}



