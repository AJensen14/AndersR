#' Plot a PCA from a normalised proteomics dataset
#'
#' Simple function to plot the results from a principal component analysis
#'
#' @param normalised_df A normalised data frame generated with AndersR
#' @param matrix Sample matrix generated with AndersR
#' @param labels A vector of labels for your samples
#'
#' @returns A ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' Plot_PCA(normalised_df, matrix, c("Old", "Young"))
#' }
Plot_PCA <- function(normalised_df, matrix, labels){

  tran_data <- t(normalised_df)
  tran_data[is.na(tran_data)] <- 0

  PCA <- prcomp(tran_data, scale = T)
  PCA_df <- PCA$x %>%
    as.data.frame() %>%
    mutate(group = matrix$group)

  pca_variance <- summary(PCA)
  pc1_variance <- pca_variance$importance[2, 1] * 100  # % variance for PC1
  pc2_variance <- pca_variance$importance[2, 2] * 100  # % variance for PC2

  PCA_df %>%
    ggplot(aes(x = PC1, y = PC2))+
    theme_minimal(base_size = 15)+
    theme(axis.title = element_text(size = 22),
          panel.grid.major = element_line(colour ="grey"),
          panel.background = element_rect(colour = "white"),
          axis.text = element_text(size = 12, colour = "black"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 13, colour = "black"))+
    geom_point(aes(fill = group),
               size = 5, shape = 21)+
    labs(x = paste("PC1 -", round(pc1_variance, 2), "% variance"),
         y = paste("PC2 -", round(pc2_variance, 2), "% variance"),
         fill = "")+
    scale_fill_manual(values = c("#edae49", "#7209b7"), labels = labels)

}
