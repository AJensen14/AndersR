#' Add in grouping information to column names
#'
#' Extra column naming function in case the original names are very long
#'
#' @param data Dataframe to change
#'
#' @return Dataframe with new columns
#' @export
#'
#' @examples
#' \dontrun{
#' Add_Grouping(df)
#' }
Add_Grouping <- function(data) {

  old_names <- colnames(data)
  new_names <- vector("character", length(old_names))

  cat("Add in group labels for each column (press Enter to skip):\n\n")

  for (i in seq_along(old_names)) {
    repeat {
      prompt <- paste0("[", i, "] ", old_names[i], " group: ")
      input <- readline(prompt)

      # Prevent duplicate name
      if (input == "") {
        new_names[i] <- old_names[i]
        break
      } else if (paste0(old_names[i], "_", input) == old_names[i]) {
        cat("Resulting name would be the same as the original. Please choose a different group label or press Enter to skip.\n")
      } else {
        new_names[i] <- paste0(old_names[i], "_", input)
        break
      }
    }
  }

  colnames(data) <- new_names
  return(data)
}

