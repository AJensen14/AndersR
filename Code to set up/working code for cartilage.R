# Testing functions (cartilage)

# Tidy spectronaut
df <- Tidy_Spectronaut(file_path = "C:/Users/hlajens2/Desktop/MS Data/KA_250618_Equine0411_Cartilage_Report.tsv")

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
Test_Normalisation(filtered_df, matrix) # Quantile and Log2 (old 50 as an outlier)

# Normalise
normalised_df <- Normalise_Data(filtered_df = filtered_df, method = "VSN")

# Do differential expression analysis
DA_results <- Differential_Abundance_Analysis(normalised_df, matrix, "old - young", "BH")

# Make volcano
volc_plot <- Volcano_Plot(DA_results, 0.05, 1)+
  labs(title = "Volcano plot (Red = Higher in old, Blue = Higher in young)")

# Summarise DA result
results <- Summarise_DA(result = DA_results, 1, 0.05)

# Save filtered_df
tissue <- "Cartilage"
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


