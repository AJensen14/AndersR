# Testing functions (Muscle)

# Tidy spectronaut
df <- Tidy_Spectronaut(file_path = "C:/Users/hlajens2/Desktop/MS Data/KA_250618_Equine0411_Muscle_Report.tsv")

# tidy column names
df <- Tidy_Column_Names(data = df)

# Add in grouping
df <- Add_Grouping(data = df)

# check data
data_check <- Check_Data(data = df,group_patterns = c("young", "old"),
           group_labels = c("young" = "Young",
                            "old" = "Old"))

# Make matrix
matrix <- Make_Matrix(data = df)
table(matrix$group)

# Filter proteins
filtered_df <- Filter_Proteins(data = df,
                               matrix = matrix,
                               min_frac = 0.7)

# test normalisation
Test_Normalisation(filtered_df, matrix)

# Normalise
normalised_df <- Normalise_Data(filtered_df = filtered_df, method = "VSN")

# Do differential expression analysis
DA_results <- Differential_Abundance_Analysis(normalised_df, matrix, "old - young", "BH")

# Make volcano
volc_plot <- Volcano_Plot(DA_results, 0.05, 0.58)+
  labs(title = "Volcano plot (Red = Higher in old, Blue = Higher in young)")

# Summarise DA result
results <- Summarise_DA(result = DA_results, 0.58, 0.05)
results$down_p

# All working
# To do
# Add in some kind of volcano legend also fix title
# interactive plotly function
# Test on all 3 MS dataset

tissue <- "Muscle"
write.csv(x = filtered_df,
          file = paste0("E:/Lab work/Progress/Omics_results/filtered_data_",
                        tissue,
                        ".csv"))
# Save DA results
write.csv(x = DA_results,
          file = paste0("E:/Lab work/Progress/Omics_results/DA_Results_",
                        tissue,
                        ".csv"))

# save data checking
ggsave(plot = data_check,
       filename = paste0("data_check_",
                         tissue,
                         ".png"),
       width = 16, height = 12,
       dpi = 600, device = "png",
       path = "E:/Lab work/Progress/Omics_results/")

# Save volcano plot
ggsave(plot = volc_plot,
       filename = paste0("volc_plot_",
                         tissue,
                         ".png"),
       width = 16, height = 16,
       dpi = 800, device = "png",
       path = "E:/Lab work/Progress/Omics_results/")


# trying to make gene to protein function
map_genes_to_proteins <- function(genes, species = "hsapiens") {
  # Connect to Ensembl
  ensembl <- useMart("ensembl", dataset = paste0(species, "_gene_ensembl"))

  # Retrieve mappings
  result <- getBM(
    attributes = c(
      "external_gene_name",
      "uniprotswissprot",
      "uniprotsptrembl",
      "description"
    ),
    filters = "external_gene_name",
    values = genes,
    mart = ensembl
  )

  # Create a combined UniProt ID column
  result$uniprot <- ifelse(
    result$uniprotswissprot != "",
    result$uniprotswissprot,
    result$uniprotsptrembl
  )

  # Keep only relevant columns
  result <- result[, c("external_gene_name", "uniprot", "description")]

  # For each gene, keep the first match

  result_clean <- result %>%
    group_by(external_gene_name) %>%
    summarise(
      UniProt_ID = first(na.omit(uniprot)),
      Protein_Description = first(na.omit(description))
    )

  # Prepare full output with all input genes
  output <- data.frame(
    Gene = genes,
    UniProt_ID = NA_character_,
    Protein_Description = NA_character_,
    stringsAsFactors = FALSE
  )

  # Fill in matches
  output <- left_join(output, result_clean, by = c("Gene" = "external_gene_name"))

  # If no match, set "not_found"
  output$UniProt_ID <- ifelse(is.na(output$UniProt_ID), "not_found", output$UniProt_ID)
  output$Protein_Description <- ifelse(is.na(output$Protein_Description), "not_found", output$Protein_Description)

  return(output)
}

genes = c("TIMP3", "ACTN3", "ACAT1", "H1-10", "BUMH")
species = "hsapiens"

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
