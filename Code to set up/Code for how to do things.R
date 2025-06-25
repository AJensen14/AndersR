# Learning how to make functions
library(usethis)

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
