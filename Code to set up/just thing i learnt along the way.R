# Remove all functions
rm(list = ls()[sapply(ls(), function(x) is.function(get(x)))])
