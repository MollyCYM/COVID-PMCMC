# Create a sequence of dates from January 30th, 2020 to April 28th, 2020
library(latex2exp)
start_date <- as.Date("2020-01-30")
end_date <- as.Date("2020-12-02")
date_sequence <- seq.Date(start_date, end_date, by = "days")

# Create a function to generate the dummy values
dummy_function <- function(time_step) {
  date_numeric <- as.numeric(time_step)
  if (time_step >= as.Date("2020-01-30") && time_step <= as.Date("2020-03-23")) {
    return(0)
  } else if (time_step >= as.Date("2020-03-24") && time_step <= as.Date("2020-06-15")) {
    return(1)
  } else if (time_step >= as.Date("2020-06-16") && time_step <= as.Date("2020-11-05")) {
    return(0)
  } 
  else if (time_step >= as.Date("2020-11-06") && time_step <= as.Date("2020-12-02")) {
    return(1)
  } 
  else {
    return(NA)  # Handle values outside the specified range
  }
}

# Apply the function to generate the dummy values for all dates
dummy_values <- sapply(date_sequence, dummy_function)

# Create a vector of colors based on the dummy values (1 = red, 0 = blue)
line_colors <- ifelse(dummy_values == 1, "red", "blue")

# Create a sequence of weekly dates
weekly_date_sequence <- seq.Date(start_date, end_date, by = "weeks")

# Initialize the plot
plot(date_sequence, dummy_values, type = "n", ylim = c(-0.1, 1.1), 
     xlab = "Weekly Dates", ylab = TeX(r"(Lockdown dummy $L_{t}$)"), xaxt = "n",
     main =TeX(r"(England lockdown policy implementation roadmap for $L_{t}$)"))


current_color <- line_colors[1]
start_date <- date_sequence[1]
for (i in 2:length(date_sequence)) {
  if (line_colors[i] != current_color || i == length(date_sequence)) {
    if (current_color == "red") {
      segments(start_date, dummy_values[i - 1], date_sequence[i - 1], dummy_values[i - 1], col = "red")
    } else {
      segments(start_date, dummy_values[i - 1], date_sequence[i - 1], dummy_values[i - 1], col = "blue")
    }
    current_color <- line_colors[i]
    start_date <- date_sequence[i]
  }
}



# Customize the x-axis labels to show week dates
date_labels <- format(weekly_date_sequence, "%b %d")
axis(1, at = weekly_date_sequence, labels = date_labels, cex.axis = 0.7, las = 2)

# Add solid vertical red lines at specified dates
abline(v = as.Date(c("2020-03-24", "2020-11-06")), col = "red")

# Add solid vertical blue line at specified date
abline(v = as.Date("2020-06-16"), col = "blue")

# Add grid lines for clarity
grid()

text(x = date_sequence[24], y = 0.06, "Initial outbreak period", col = "blue",cex = 0.6)
text(x = date_sequence[98], y = 0.9, "First National Lockdown", col = "red",cex = 0.6)
text(x = date_sequence[203], y = 0.06, "Lifting and eased period", col = "blue",cex = 0.6)
text(x = date_sequence[277], y = 0.9, "Second National Lockdown", col = "red",cex = 0.6)