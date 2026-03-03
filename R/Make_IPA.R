#' Function to make an IPA suitable file
#'
#' @param de_file Differential expression data generated with the AndersR package
#' @param tsv_file The TSV file containing all data (Make sure the gene names are the rownames)
#'
#' @returns An IPA suitable dataset
#' @export
#'
#' @examples
#' \dontrun{
#' IPA_data <- Make_IPA("file_path_de.csv", "file_path_all_data.tsv")
#' }
Make_IPA <- function(de_file, tsv_file){

  # Read in datasets
  daa <- read.csv(file = de_file)
  overall <- Tidy_Spectronaut(file_path =  tsv_file)

  # find out what is lost in DAA results
  lost <- setdiff(row.names(overall), daa$X)

  # Make a dataframe similar to the daa file
  empty_data <- data.frame(X = lost,
                           logFC = 0,
                           AveExpr = NA,
                           t = NA,
                           P.Value = 1,
                           adj.P.Val = 1,
                           B = NA)
  # bind rows
  total_data <- rbind(daa,empty_data)

  # return the data
  return(total_data)

}
