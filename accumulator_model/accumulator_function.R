check_time <- function(target_time, threshold, k, dstep) {
  max_iters <- ceiling(target_time / dstep)
  actual_time <- numeric(max_iters)
  est_time <- numeric(max_iters)
  
  cur_est_time <- 0
  cur_time <- 0
  i <- 0
  
  while ((target_time * threshold > cur_est_time) & cur_time < target_time) {
    i <- i + 1
    cur_time <- cur_time + dstep
    check_times[i] <- cur_time
    cur_est_time <- rnorm(1, cur_est_time + dstep, sqrt(k) * cur_time)
    est_time[i] <- cur_est_time
  }
  
  check_times <- check_times[1:i]
  est_time <- est_time[1:i]
  
  list(check_times = check_times, 
       est_time = est_time, 
       checked = ifelse(max(check_times) < target_time, 1, 0),
       params = list(
         target_time = target_time,
         threshold = threshold, 
         k = k, 
         dstep = dstep
       ))
}
block_checks <- function(target_time, threshold, rel_k, dstep) {
  # Calculate parameter k based on relative k and target time
  k <- rel_k / target_time
  
  # Predefine maximum iterations
  max_iters <- ceiling(target_time / dstep)
  
  # Initialize result arrays
  actual_time <- numeric(max_iters)
  clock_check <- numeric(max_iters)
  est_time <- numeric(max_iters)
  
  # Initialize variables for tracking time
  cur_est_time <- 0
  cur_time <- 0
  i <- 0
  time_to_target <- target_time - cur_time
  estimated_time_since_last_cc <- 0
  time_since_last_cc <- 0
  time_of_last_cc <- 0
  
  # Main simulation loop
  while ((time_to_target * threshold > estimated_time_since_last_cc) & cur_time < target_time) {
    # Inner loop for time updates
    while ((time_to_target * threshold > estimated_time_since_last_cc) & cur_time < target_time) {
      i <- i + 1
      cur_time <- cur_time + dstep
      time_since_last_cc <- time_since_last_cc + dstep
      
      # Record actual time and estimate time
      actual_time[i] <- cur_time
      cur_est_time <- rnorm(1, cur_est_time + dstep, sqrt(k) * time_since_last_cc)
      est_time[i] <- cur_est_time
      estimated_time_since_last_cc <- cur_est_time - time_of_last_cc
    }
    
    # Update parameters for clock check
    time_to_target <- target_time - cur_time
    clock_check[i] <- 1
    time_since_last_cc <- 0
    time_of_last_cc <- cur_time
    cur_est_time <- cur_time
    estimated_time_since_last_cc <- 0
  }
  
  # Trim excess elements
  actual_time <- actual_time[1:i]
  est_time <- est_time[1:i]
  clock_check <- clock_check[1:i]
  
  # Return results
  list(
    actual_time = actual_time,
    est_time = est_time,
    clock_check = clock_check,
    checked = as.integer(max(actual_time) < target_time),
    params = list(
      target_time = target_time,
      threshold = threshold,
      k = k,
      dstep = dstep
    )
  )
}

simulate_clock_checks <- function(n_simulations, target_time, threshold, rel_k, dstep) {
  checked <- numeric(n_simulations)
  clock_check_times <- c()
  
  for (sim in 1:n_simulations) {
    # Run the simulation using the block_checks function
    result <- block_checks(target_time, threshold, rel_k, dstep)
    
    # Extract the times where clock checks occurred
    clock_check_times <- c(clock_check_times, result$actual_time[which(result$clock_check == 1)])
    
    checked[sim] <- result$checked
  }
  
  print(mean())
  # Return the times of clock checks as a vector
  return(clock_check_times)
}


plot_time_est <- function(result) {
  # Set up the plot
  plot(result$actual_time, result$est_time, xlim = c(0, result$params$target_time), 
       type = "l", col = "blue", lwd = 2, 
       xlab = "Actual Time", ylab = "Estimated Time", 
       main = "Actual vs Estimated Time with Clock Checks")
  
  # Add vertical lines for clock checks
  abline(v = result$actual_time[which(result$clock_check == 1)], col = "red", lty = 2, lwd = 1)
  
  # Add a reference line (y=x) for perfect estimation
  abline(a = 0, b = 1, col = "gray", lwd = 2, lty = 2)
}

# Run 100 simulations and track the `cur_time` when `est_time` is max
run_simulations <- function(n_simulations, target_time, threshold, k, dstep) {
  max_times <- numeric(n_simulations)
  checked <- numeric(n_simulations)
  
  for (sim in 1:n_simulations) {
    result <- check_time(target_time, threshold, k, dstep)
    check_times <- result$check_times
    est_times <- result$est_times
    
    # Find the index of the maximum `est_time`
    max_index <- which.max(est_times)
    max_times[sim] <- check_times[max_index]
    
    checked[sim] <- result$checked
  }
  
  print(mean(checked))
  print(mean(max_times))
  return(max_times[which(checked == 1)])
}

# Parameters
target_time <- 300
threshold <- 0.75
k <- 0.001
dstep <- 1
n_simulations <- 1000

plot_time_est(check_time(target_time, threshold, k, dstep))

# Run the simulations
max_times <- run_simulations(n_simulations, target_time, threshold, k, dstep)
print(paste0("Percentage of checks: ", length(max_times)/n_simulations))

# Summary of results
hist(max_times, breaks = 100)
mean(max_times)
target_time * threshold
