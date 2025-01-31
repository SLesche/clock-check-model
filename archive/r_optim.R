clock_check_lik <- function(x, t_target, k, g, threshold, a, m) {
  # Step size for the trapezoidal rule
  delta_x <- (x - 0.01) / m  # This is the width of each small interval
  
  # Precompute normal CDF values for the range
  normal_cdf_left <- numeric(m)
  normal_cdf_right <- numeric(m)
  
  # Calculate normal CDF values for each segment
  for (i in 1:m) {
    t_left <- 0.01 + (i - 1) * delta_x
    t_right <- 0.01 + i * delta_x
    
    # Store normal CDF values in arrays
    normal_cdf_left[i] <- pnorm(t_target, mean = t_left, sd = k * t_left)
    normal_cdf_right[i] <- pnorm(t_target, mean = t_right, sd = k * t_right)
  }
  
  # Compute the integral using the trapezoidal rule
  integral <- 0
  for (i in 1:(m - 1)) {
    hazard_left <- g + (1 - g) / (1 + exp(-a * (1 - normal_cdf_left[i]) - threshold))
    hazard_right <- g + (1 - g) / (1 + exp(-a * (1 - normal_cdf_left[i + 1]) - threshold))
    
    # Apply trapezoidal rule: (f(t_left) + f(t_right)) * delta_x / 2
    integral <- integral + (hazard_left + hazard_right) * delta_x / 2
  }
  
  # Compute the survival probability and action probability
  survival_prob <- 1 - max(integral / m, 0)  # Ensure no negative integral
  action_prob <- g + (1 - g) / (1 + exp(-a * (1 - pnorm(t_target, mean = x, sd = k * x)) - threshold))
  
  # Return the log-likelihood
  return(log(max(survival_prob * action_prob, 1e-10)))  # Avoid log(0)
}


log_likelihood <- function(params, known_t_to_target, clock_check_time, m) {
  k <- params[1]
  g <- params[2]
  threshold <- params[3]
  a <- params[4]
  
  total_log_likelihood <- 0
  for (n in 1:length(clock_check_time)) {
    total_log_likelihood <- total_log_likelihood + 
      clock_check_lik(clock_check_time[n], known_t_to_target[n], k, g, threshold, a, m)
  }
  
  return(-total_log_likelihood)  # Return negative log-likelihood for minimization
}

log_likelihood(c(1, 0.6, 0.3, 10), known_t_to_target, clock_check_time, 100)

data <- read.csv("archive/diffusion_data.csv")

fit_data <- data %>% 
  filter(known_t_to_target != 0, time_since_last_cc != 0) %>% 
  filter(known_t_to_target / block_duration == 1)

clock_check_time = fit_data$time_since_last_cc
known_t_to_target = fit_data$known_t_to_target

m <- 100  # Precision of integral estimation

# Initial parameter guess
initial_params <- c(k = 1, g = 0.1, threshold = 0.1, a = 10)


# Run the optimization
result <- optim(
  par = initial_params, 
  fn = log_likelihood, 
  known_t_to_target = known_t_to_target, 
  clock_check_time = clock_check_time, 
  m = m, 
  method = "BFGS"  # Use a suitable optimization method
)

# Optimized parameters
optimized_params <- result$par
optimized_params

