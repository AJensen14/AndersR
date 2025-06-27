#' Normalisation function
#'
#' This function utilises the NormalyserDE package to help with normalisation
#'
#' @param filtered_df A filtered dataframe
#' @param method The method of normalisation (Median, Log2, Quantile or VSN)
#'
#' @return A normalised dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' Normalise_Data(filtered_df, method = "Median")
#' }
Normalise_Data <- function(filtered_df, method){

  matrix_to_norm <- as.matrix(filtered_df)

  method <- tolower(method)
  method <- str_to_title(method)

  if (method == "Median"){

    result <- medianNormalization(matrix_to_norm)
    result <- as.data.frame(result)
    rownames(result) <- rownames(filtered_df)

  } else if (method == "Log2"){

    result <- log2(matrix_to_norm)
    result <- as.data.frame(result)
    rownames(result) <- rownames(filtered_df)

  } else if (method == "Quantile"){

    result <- performQuantileNormalization(matrix_to_norm)
    result <- as.data.frame(result)
    rownames(result) <- rownames(filtered_df)

  } else if (method == "Vsn"){

    result <- performVSNNormalization(matrix_to_norm)
    result <- as.data.frame(result)
    rownames(result) <- rownames(filtered_df)

  } else {

    print("Please select one of the following: Median, Log2, Quantile or VSN")
    print("If you want to do another method this package can't handle it and I aint coding it")

  }

  return(result)

}

