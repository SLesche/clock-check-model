set.seed(42)

# Simulation parameters
n_participants <- 50   # Number of simulated people
max_time <- 300        # 5 minutes in seconds
dt <- 1                # Time step (1 second)
k <- 1.2             # Uncertainty growth rate
lambda <- 0.05        # Cost penalty for error
check_cost <- 0.1      # Small cost for checking the clock

# Function to simulate a participant's behavior
simulate_participant <- function() {
  time <- seq(0, max_time, by = dt)  # Time steps
  sigma_sq <- k * time               # Growing variance
  checks <- c()                      # Store time points of checking
  
  for (t in time) {
    # Compute expected free energy components
    info_gain <- 0.5 * log(k * t + 1e-6)  # Avoid log(0)
    expected_cost <- lambda * ((t - 300)^2 + k * t)
    
    # Decision rule: check if uncertainty exceeds cost
    if (info_gain > expected_cost + check_cost) {
      checks <- c(checks, t)  # Store time of check
    }
  }
  
  return(checks)
}

# Simulate multiple participants
all_checks <- lapply(1:n_participants, function(x) simulate_participant())

# Convert to a data frame for plotting
library(ggplot2)
library(dplyr)
check_df <- do.call(rbind, lapply(1:n_participants, function(i) {
  if (length(all_checks[[i]]) > 0) {
    data.frame(participant = i, check_time = all_checks[[i]])
  } else {
    data.frame(participant = i, check_time = NA)
  }
}))

# Plot the results
ggplot(check_df, aes(x = check_time)) +
  geom_histogram(binwidth = 10, fill = "blue", alpha = 0.6, color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Time Checks",
       x = "Time (seconds)",
       y = "Frequency of Clock Checks")
  