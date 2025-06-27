# Learning how to make functions
library(usethis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(readr)
library(tibble)
library(stringr)
library(cowplot)
library(NormalyzerDE)
library(SummarizedExperiment)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(limma)
library(scales)
library(colorspace)

# Adding in the protein assay function
use_r("Protein_Assay")
# This opens a new tab and makes the file in the correct place
# Then you write the code for the function
# Then you click inside the function and press code --> insert roxygen species
# Then fill in the information
# Then in the CONSOLE
# run devtools::document() to save the function
# run devtools::load_all() to load the functions
# run ?function_name to inspect the details on the function

# Next thing to do is to push to github
# firstly i am going to make a readme file to explain the package
usethis::use_readme_rmd()
# Then write in the read me file using devtools::readme()
usethis::create_github_token()
gitcreds::gitcreds_set()
# Then run usethis::use_github()

# To add more functions --> use_r() makes the function as before

# Then do the roxygen skeleton

# document it

# then commit in the git tab in the top right

# make a reference file for each function with examples


# Making new functions
usethis::use_r("Tidy_Spectronaut")
usethis::use_r("Tidy_Column_Names")
usethis::use_r("Add_Grouping")
usethis::use_r("Check_Data")
usethis::use_r("Make_Matrix")
usethis::use_r("Filter_Proteins")
usethis::use_r("Test_Normalisation")
usethis::use_r("Normalise_Data")
usethis::use_r("Differential_Abundance_Analysis")
usethis::use_r("Volcano_Plot")
usethis::use_r("Summarise_DA")

# This is a good start but need to clean up
# Quantile norm loses gene names - Fixed
# Function to map protein names to gene names - done
usethis::use_r("Gene_To_Protein")

# Function to plot effect of increasing stringency of missingness
usethis::use_r("Check_Filtering")

# Summarise DA function needs better names in the output list - done
# Then function to boxplot a protein
usethis::use_r("Boxplot_Proteins")
# function to make dummy data
usethis::use_r("Dummy_Data")

# Then push to github
# Then write in the read me file using devtools::readme()
usethis::create_github_token()
gitcreds::gitcreds_set()
# Then run usethis::use_github()
usethis::use_github()
# Re-purpose shiny app with these functions
# Re-develop app
# Send mandy sample metadata
# Go through examples with mandy


