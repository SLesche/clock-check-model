plot_likelihood <- function(k, g, threshold, a, m, t_min = 0.01, t_max = 1) {
  # Helper function: Normal CDF
  normal_cdf <- function(mean, x, sd) {
    pnorm(x, mean = mean, sd = sd)
  }
  
  # Hazard rate function
  hazard_rate <- function(x, k, g, threshold, a) {
    cdf_val <- normal_cdf(1, x, k * x)  # target_time is 1
    g + (1 - g) / (1 + exp(-a * (1 - cdf_val) - threshold))
  }
  
  # Survival function
  survival_function <- function(x, k, g, threshold, a, m, t_min) {
    delta_x <- (x - t_min) / m
    t_values <- seq(t_min, x, length.out = m + 1)
    
    # Compute integral of hazard rates
    integral <- 0
    for (i in seq_len(m)) {
      hazard_left <- hazard_rate(t_values[i], k, g, threshold, a)
      hazard_right <- hazard_rate(t_values[i + 1], k, g, threshold, a)
      integral <- integral + (hazard_left + hazard_right) * delta_x / 2
    }
    
    return(1 - integral)
  }
  
  # Compute likelihood for each time point
  time_points <- seq(t_min, t_max, length.out = 500)  # Fine-grained time steps
  likelihoods <- sapply(time_points, function(x) {
    h_x <- hazard_rate(x, k, g, threshold, a)
    S_x <- survival_function(x, k, g, threshold, a, m, t_min)
    log(h_x * S_x)
  })
  
  # Plot the likelihood
  plot(
    time_points, likelihoods, type = "l", col = "blue", lwd = 2,
    xlab = "Time (t)", ylab = "Likelihood", 
    main = "Likelihood of Response Time"
  )
  grid()
}

# Example usage
plot_likelihood(
  k = 1000,       # Scaling parameter
  g = 0.1,       # Baseline hazard
  threshold = 0.0001,  # Threshold parameter
  a = 10,        # Steepness parameter
  m = 100        # Precision of integral estimation
)
