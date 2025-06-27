# code to start making an R package
# packages
library(usethis)

# setwd
setwd("E:/Anders R package/")

# Running this makes a new r project called andersR in the set directory
create_package("AndersR")

# Checking this file theres loads of stuff in it
# So now I should go to that 

# test installing
devtools::install_github("AJensen14/AndersR")
library(AndersR)

# Working