#' Code to tidy up SpectroNaut output for analysis
#'
#' This function reads a Spectronaut report (TSV file) and returns a cleaned
#' dataframe of protein quantities.
#'
#' @param file_path The location your file is saved
#'
#' @return A tidy dataframe from SpectroNaut
#' @export
#'
#' @examples
#' \dontrun{
#' tidy_df <- Tidy_Spectronaut("data/spectronaut_output.tsv")
#' head(tidy_df)
#' }
Tidy_Spectronaut <- function(file_path){

  df_report <- read_tsv(file = file_path) %>%
    mutate(Naming = case_when(is.na(PG.Genes)~PG.UniProtIds,
                              !is.na(PG.Genes)~PG.Genes),
           Naming = str_remove(Naming, "^;"),
           Naming = unlist(lapply(strsplit(Naming, ";"), function(x)x[1])),
           Naming = make.unique(Naming),
           Naming = case_when(Naming == ""~PG.UniProtIds,
                              Naming != ""~Naming))

  # tidy up data
  df_report_samples <- read_tsv(file = file_path) %>%
    dplyr::select(contains("Quantity")) %>%
    dplyr::select(!contains("Log"))

  # Change NaN to NA
  df_report_samples[df_report_samples == "NaN"] <- NA

  # Add in the row names
  rownames(df_report_samples) <- df_report$Naming

  return(df_report_samples)

}
