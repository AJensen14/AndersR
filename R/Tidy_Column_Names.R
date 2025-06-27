#' Change column names to something more informative
#'
#' Function to help adjust column names in a more interactive way.
#'
#' @param data The dataset to change columns on
#'
#' @return A dataframe with new columns
#' @export
#'
#' @examples
#' #' \dontrun{
#' colnames(df)
#' tidy_df <- Tidy_Column_Names(df)
#' colnames(df)
#' }
Tidy_Column_Names <- function(data){

  old_col_names <- colnames(data)
  new_col_names <- vector("character", length(old_col_names))

  cat("Rename columns (press Enter to keep current name, but you can't retype the same name):\n")

  for (i in seq_along(old_col_names)) {
    repeat {
      prompt <- paste0("[", i, "] ", old_col_names[i], ": ")
      input <- readline(prompt)

      if (!nzchar(input)) {
        # User pressed enter: keep old name
        new_col_names[i] <- old_col_names[i]
        break
      } else if (input == old_col_names[i]) {
        cat("You typed the same name. Please enter a new name or press Enter to keep the current one.\n")
      } else {
        new_col_names[i] <- input
        break
      }
    }
  }

  colnames(data) <- new_col_names
  return(data)
}

