# Function needs to return a dataframe containing week info
# you input target km and

Weekly_Runs <- function(target_km){

  days_of_week <- c("Monday",
                    "Tuesday",
                    "Wednesday",
                    "Thursday",
                    "Friday",
                    "Saturday",
                    "Sunday")

  distances <- numeric(length(days_of_week))

  for (i in seq_along(days_of_week)) {
    input <- readline(
      paste("Enter distance for", days_of_week[i], "(km): ")
    )
    # Convert to numeric
    distances[i] <- as.numeric(input)
  }

  week_data <- data.frame(
    day = days_of_week,
    distance = distances
  )

  total_km <- sum(week_data$distance)

  if (total_km == target_km){

    print("That plan works for your target")

  } else if (total_km > target_km){

    print(paste("That plan means you are running", total_km-target_km,
                "too many kilometers."))

  } else if (total_km < target_km){

    print(paste("That plan means you are running", target_km-total_km,
                "too few kilometers."))

  }

  return(week_data)

}

