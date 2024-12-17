check_time <- function(target_time, threshold, k, dstep) {
  max_iters <- ceiling((target_time - threshold) / dstep)
  check_times <- numeric(max_iters)
  est_times <- numeric(max_iters)
  
  cur_est_time <- 0
  cur_time <- 0
  i <- 0
  
  while ((target_time - threshold > cur_est_time) & cur_time < target_time) {
    i <- i + 1
    cur_time <- cur_time + dstep
    check_times[i] <- cur_time
    cur_est_time <- rnorm(1, cur_est_time + dstep, k * cur_time)
    est_times[i] <- cur_est_time
  }
  
  check_times <- check_times[1:i]
  est_times <- est_times[1:i]
  
  list(check_times = check_times, 
       est_times = est_times, 
       checked = ifelse(max(check_times) < target_time, 1, 0),
       params = list(
         target_time = target_time,
         threshold = threshold, 
         k = k, 
         dstep = dstep
       ))
}

plot_time_est <- function(result){
  plot(result$check_times, result$est_times, xlim = c(0, result$params$target_time), type = "l")
  abline(h = result$params$target_time - result$params$threshold, col = "red", lwd = 2, lty = 2)
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
  return(max_times[which(checked == 1)])
}

# Parameters
target_time <- 300
threshold <- 100
k <- 0.02
dstep <- 1
n_simulations <- 1000

plot_time_est(check_time(target_time, threshold, k, dstep))

# Run the simulations
max_times <- run_simulations(n_simulations, target_time, threshold, k, dstep)

# Summary of results
hist(max_times, breaks = 50)

