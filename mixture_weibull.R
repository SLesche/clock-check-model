library(rstan)
library(dplyr)

data <- read.csv("archive/diffusion_data.csv")

clean_data <- data %>% 
  # filter(
  #   block_duration == known_t_to_target
  # ) %>%
  mutate(
    r = time_since_last_cc / known_t_to_target
  ) %>% 
  filter(r > 0.01) # make sure that the data is not too close to 0 (this cause issues in integration)

# Set up the data
# Replace with your actual data
stan_data <- list(
  N = nrow(clean_data),  # Number of events per participant
  clock_check_time = clean_data$r,
  known_t_to_target = clean_data$known_t_to_target,
  group = clean_data$participant,
  G = length(unique(clean_data$participant))
)

# Fit the model
options(mc.cores = 4)
fit <- stan(
  file = "weibull_model_mixed_timedep.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
  # control=list(adapt_delta=0.99, stepsize = 0.01, max_treedepth =15)
)

# saveRDS(fit, "stanfit_weibull_model.rds")
# Check convergence with traceplots
traceplot(fit)

# Get a summary of the posterior distributions
print(fit)

# Extract posterior samples
posterior_samples <- extract(fit)

# Function to generate predicted distributions from posterior samples
generate_ppc <- function(fit, times) {
  # Extract posterior samples
  posterior_samples <- extract(fit) # Convert to a data frame if necessary
  
  # Get parameter names (assuming standard naming conventions)
  pi_samples <- posterior_samples$pi
  alpha1_samples <- posterior_samples$alpha1
  sigma1_samples <- posterior_samples$sigma1
  alpha2_samples <- posterior_samples$alpha2
  sigma2_samples <- posterior_samples$sigma2
  
  # Number of posterior samples
  n_samples <- length(pi_samples)
  
  # Initialize a matrix to store predictive densities
  pred_matrix <- matrix(0, nrow = n_samples, ncol = length(times))
  
  # Loop over posterior samples to compute predictive densities
  for (i in 1:n_samples) {
    pi <- pi_samples[i]
    alpha1 <- alpha1_samples[i]
    sigma1 <- sigma1_samples[i]
    alpha2 <- alpha2_samples[i]
    sigma2 <- sigma2_samples[i]
    
    # Mixture Weibull densities for each time point
    f1 <- dweibull(times, shape = alpha1, scale = sigma1) # Constant hazard Weibull
    f2 <- dweibull(times, shape = alpha2, scale = sigma2) # Flexible hazard Weibull
    
    # Mixture density
    pred_matrix[i, ] <- pi * f1 + (1 - pi) * f2
  }
  
  # Summarize predictive densities
  pred_mean <- colMeans(pred_matrix)
  pred_ci <- apply(pred_matrix, 2, quantile, probs = c(0.025, 0.975)) # 95% CI
  
  # Return results as a list
  return(list(
    times = times,
    pred_mean = pred_mean,
    pred_ci_lower = pred_ci[1, ],
    pred_ci_upper = pred_ci[2, ]
  ))
}

# Example usage
# Simulate a fitted model (replace with your actual fitted model)
# Define time points for predictions
times <- seq(0.01, 1.5, length.out = 1000)

# Generate posterior predictive distribution
ppc_results <- generate_ppc(fit, times)

# Calculate density of observed data
observed_density <- density(clean_data$r)

# Plot posterior predictive results
plot(ppc_results$times, ppc_results$pred_mean, type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Density", main = "Posterior Predictive Check with Observed Data")
lines(ppc_results$times, ppc_results$pred_ci_lower, col = "red", lty = 2)
lines(ppc_results$times, ppc_results$pred_ci_upper, col = "red", lty = 2)

# Overlay observed data density
lines(observed_density$x, observed_density$y, col = "black", lwd = 2, lty = 3)

# Add legend
legend("topright", legend = c("Mean Predictive Density", "95% CI", "Observed Data"),
       col = c("blue", "red", "black"), lty = c(1, 2, 3), lwd = 2)
