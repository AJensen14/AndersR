library(AndersR)

dummy <- Dummy_Data(n_proteins = 1000, group1 = "control", "treated", n_replicates = 5, missing_frac = 0.1)

# check data
data_check <- Check_Data(data = dummy,group_patterns = c("control", "treated"),
                         group_labels = c("control" = "Con",
                                          "treated" = "Treat"))
data_check

# Make matrix
matrix <- Make_Matrix(data = dummy)
table(matrix$group)

# checked filtered
Check_Filtering(data = dummy, matrix = matrix, n_samples = 5)

# filter out proteins
filtered_df <- Filter_Proteins(data = dummy,
                               matrix = matrix,
                               min_frac = 0.75)

Test_Normalisation(filtered_df = filtered_df, sample_matrix = matrix)

# do norm
norm_data <- Normalise_Data(filtered_df = filtered_df, method = "quantile")

# Do DE analysis
?Differential_Abundance_Analysis
DE_result <- Differential_Abundance_Analysis(norm_data,
                                             matrix,
                                             contrast_formula = "treated - control",
                                             adjust.method = "BH")

# Make volc
Volcano_Plot(DE_result, pval_threshold = 0.05, logFC_threshold = 1)

# boxplot_protein
Boxplot_Proteins(data = norm_data, DA_result = DE_result,
                 matrix = matrix, proteins = c("Protein_426", "Protein_894"),
                 group1 = "control", group2 = "treated", p_type = "raw")

# summarise
result <- Summarise_DA(result = DE_result, logFC_threshold = 1, p_value_threshold = 0.05)
result$up_p_adj
result$down_p_adj
Boxplot_Proteins(data = norm_data, DA_result = DE_result,
                 matrix = matrix, proteins = result$down_p_adj[1:10],
                 group1 = "control", group2 = "treated", p_type = "adj")



