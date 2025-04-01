library(rstan)
library(dplyr)

data <- read.csv("archive/diffusion_data.csv")

clean_data <- data %>% 
  filter(
    block_duration == known_t_to_target
  ) %>%
  mutate(
    r = time_since_last_cc / known_t_to_target
  ) %>% 
  filter(r > 0.01) # make sure that the data is not too close to 0 (this cause issues in integration)

# Set up the data
# Replace with your actual data
stan_data <- list(
  N = nrow(clean_data),  # Number of events per participant
  clock_check_time = clean_data$r,
  group = clean_data$participant,
  G = length(unique(clean_data$participant))
)

# Fit the model
options(mc.cores = 4)
fit <- stan(
  file = "weibull_model.stan",
  data = stan_data,
  iter = 4000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 1000, # Number of warmup iterations
  # control=list(adapt_delta=0.99, stepsize = 0.01, max_treedepth =15)
)

# saveRDS(fit, "stanfit_weibull_model.rds")
# Check convergence with traceplots
traceplot(fit)

# Get a summary of the posterior distributions
print(fit)

# Simulate data using the 3PL Weibull PDF
simulate_3pl_weibull <- function(a, b, c, N) {
  # Generate random samples from the modified Weibull distribution
  t_samples <- numeric(N)
  
  for (i in 1:N) {
    # Inverse CDF sampling approach or rejection sampling
    u <- runif(1) # Uniform random variable
    t <- 0        # Initial guess for time
    
    # Use a numerical approach to solve for t satisfying F(t) = u
    while (TRUE) {
      f_t <- a * (b + c * t) * t^(b - 1) * exp(c * t) * exp(-a * t^b * exp(c * t)) # PDF
      F_t <- 1 - exp(-a * t^b * exp(c * t))                                       # CDF
      
      if (abs(F_t - u) < 1e-6) break
      t <- t + 0.001 # Increment time to converge
    }
    
    t_samples[i] <- t
  }
  
  return(t_samples)
}

# Load your fitted Stan model (assumed already run)
# Example: fit <- stan_model("path/to/your/model.stan")
# Extract posterior samples
posterior_samples <- extract(fit)

# Posterior predictive sampling
N <- length(clean_data$r) # Number of observations
n_draws <- 500                # Number of posterior draws for PPC

ppc_results <- list()
for (i in 1:n_draws) {
  # Draw parameters from posterior
  a <- posterior_samples$a[i]
  b <- posterior_samples$b[i]
  c <- posterior_samples$c[i]
  
  # Generate predictive samples
  ppc_results[[i]] <- simulate_3pl_weibull(a, b, c, N)
}

# Compute predictive mean and credible intervals
ppc_samples <- do.call(rbind, ppc_results)
ppc_mean <- apply(ppc_samples, 2, mean)
ppc_ci_lower <- apply(ppc_samples, 2, quantile, probs = 0.025)
ppc_ci_upper <- apply(ppc_samples, 2, quantile, probs = 0.975)

# Plot PPC results
ggplot() +
  geom_histogram(aes(x = clock_check_time, y = ..density..), bins = 30, fill = "gray", alpha = 0.5) +
  geom_line(aes(x = clock_check_time, y = ppc_mean), color = "blue", size = 1) +
  geom_ribbon(aes(x = clock_check_time, ymin = ppc_ci_lower, ymax = ppc_ci_upper), fill = "red", alpha = 0.3) +
  labs(title = "Posterior Predictive Check", x = "Time", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::percent)
# Function to generate predicted distributions from posterior samples
generate_ppc <- function(fit, times) {
  # Extract posterior samples
  posterior_samples <- extract(fit) # Convert to a data frame if necessary
  
  # Get parameter names (assuming standard naming conventions)
  pi_samples <- posterior_samples$pi
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
    sigma1 <- sigma1_samples[i]
    alpha2 <- alpha2_samples[i]
    sigma2 <- sigma2_samples[i]
    
    # Mixture Weibull densities for each time point
    f1 <- dweibull(times, shape = 1, scale = sigma1) # Constant hazard Weibull
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
