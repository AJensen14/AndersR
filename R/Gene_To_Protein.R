#' Function to convert gene names into protein names
#'
#' @param DA_result The result from DA analysis
#' @param species The species to look in. First letter of the genus then the species
#'
#' @return A new DA data with protein description included
#' @export
#'
#' @examples
#' \dontrun{
#' DA_results <- Differential_Abundance_Analysis(normalised_df,
#'                                               matrix,
#'                                              "old - young",
#'                                               adjust.method = "BH")
#' Updated_DA <- Gene_To_Protein(DA_results, "hsapiens")
#'
#' }
Gene_To_Protein <- function(DA_result, species = "hsapiens"){

  ensembl <- useMart("ensembl", dataset = paste0(species, "_gene_ensembl"))
  result <- getBM(
    attributes = c(
      "external_gene_name",
      "uniprotswissprot",
      "uniprotsptrembl",
      "description"),
    filters = "external_gene_name",
    values = rownames(DA_results),
    mart = ensembl)

  result$uniprot <- ifelse(
    result$uniprotswissprot != "",
    result$uniprotswissprot,
    result$uniprotsptrembl
  )

  result_clean <- result %>%
    group_by(external_gene_name) %>%
    summarise(UniProt_ID = paste(unique(na.omit(uniprot)), collapse = "; "),
              Protein_Description = paste(unique(na.omit(description)), collapse = "; "))

  output <- data.frame(Gene = rownames(DA_results),
                       stringsAsFactors = FALSE)

  # Merge, preserving all input genes
  output <- left_join(
    output,
    result_clean,
    by = c("Gene" = "external_gene_name"))

  output$UniProt_ID[is.na(output$UniProt_ID)] <- "not_found"
  output$Protein_Description[is.na(output$Protein_Description)] <- "not_found"

  DA_new <- DA_result %>%
    mutate(Prot_name = output$Protein_Description,
           UniProt_ID = output$UniProt_ID)

  return(DA_new)

}
