% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tidy_Column_Names.R
\name{Tidy_Column_Names}
\alias{Tidy_Column_Names}
\title{Change column names to something more informative}
\usage{
Tidy_Column_Names(data)
}
\arguments{
\item{data}{The dataset to change columns on}
}
\value{
A dataframe with new columns
}
\description{
Function to help adjust column names in a more interactive way.
}
\examples{
#' \dontrun{
colnames(df)
tidy_df <- Tidy_Column_Names(df)
colnames(df)
}
}
