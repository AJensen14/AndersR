#' Title A function to load required packages for AndersR
#'
#' @return Nothing just loads packages. Will come up with a message if some of the packages are not installed.
#' @export
#'
#' @examples
#' \dontrun{
#' Load_Packages()
#' }
Load_Packages <- function(){

  pkgs <- c(
    "usethis", "dplyr", "ggplot2", "reshape2", "readr", "tibble", "stringr",
    "cowplot", "NormalyzerDE", "SummarizedExperiment", "tidyverse", "rstatix",
    "ggpubr", "ggrepel", "limma", "scales", "colorspace", "clusterProfiler",
    "org.Hs.eg.db"
  )

  not_installed <- pkgs[!pkgs %in% rownames(installed.packages())]

  if (length(not_installed) > 0) {
    stop("The following required package(s) are not installed:\n",
         paste0(" - ", not_installed, collapse = "\n"),
         call. = FALSE)
  }

  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))

}
