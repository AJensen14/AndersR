#' Function to run all the Gene Ontology analysis at once
#'
#' @param group_up The group that corresponds to a positive fold change
#' @param group_down The group that corresponds to a negative fold change
#' @param DA_summary The DA summary object produced using the AndersR package
#' @param save_file The location to save images
#' @param name Name to help with saving files
#'
#' @returns Multiple ggplots from clusterProfiler
#' @export
#'
#' @examples
#' \dontrun{
#' GO_from_DA_all("Old", "Young", DA_summary, "C:/file_path/plots/", "Experiment_1")
#' }
GO_from_DA_all <- function(group_up, group_down, DA_summary, save_file, name = NULL){

  if (is.null(name)) name <- "GO"

  # run loop
  for (i in c("MF", "BP", "CC")) {

    GO_from_DA(DA_summary = DA_summary, ontology = i, use_adj_p = FALSE, group_up = group_up, group_down = group_down, top_n = 15)

    ggsave(filename = paste0(name,"_",i, "dotplot.png"),
           path = paste0(save_file),
           width = 14,
           height = 11,
           dpi = 600)
  }

}
