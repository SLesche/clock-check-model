library(dplyr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

data <- read.csv("archive/diffusion_data.csv")

fit_data <- data %>% 
  filter(known_t_to_target != 0, time_since_last_cc != 0)


stan_data <- list(
  N = nrow(fit_data),
  clock_check_time = fit_data$time_since_last_cc,
  known_t_to_target = fit_data$known_t_to_target,
  m = 100
)

# Fit the model
fit <- stan(
  file = "activation_model.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
  # init = list(
  #   list(
  #     # Priors
  #     mu_k = rnorm(1, 1, 0.5),
  #     sigma_k = rgamma(1, 1, 1),
  #     mu_sigma_0 = rgamma(1, 1, 1),
  #     sigma_sigma_0 = rgamma(1, 1, 1),
  #     mu_raw_theta = runif(1, 0, 0.8),
  #     sigma_raw_theta = runif(1, 0, 0.2),
  #     mu_sigma_err = rbeta(1, 2, 6),
  #     sigma_sigma_err = rbeta(1, 2, 10)
  #   )
  # )
)

# R-hat statistics
print(fit)

# Trace plots for parameters
traceplot(fit)

# Extract posterior samples from the fitted model
posterior_samples <- rstan::extract(fit)


# Posterior samples for sigma_0, k, and threshold
g_samples <- posterior_samples$g
k_samples <- posterior_samples$k
threshold_samples <- posterior_samples$threshold

# Visualize the posteriors of parameters
ggplot(data.frame(g_samples), aes(x = g_samples)) + 
  geom_density() + ggtitle("Posterior of sigma_0")

ggplot(data.frame(k_samples), aes(x = k_samples)) + 
  geom_density() + ggtitle("Posterior of k")

ggplot(data.frame(threshold_samples), aes(x = threshold_samples)) + 
  geom_density() + ggtitle("Posterior of threshold")

# threshold_growth_samples <- posterior_samples$theta_rate
t_target_samples <- data$known_t_to_target / data$block_duration  # the target times for each observation

# Number of posterior samples
N_simulations <- 100

# Initialize a matrix to store the simulated response times
simulated_data <- matrix(NA, nrow = N_simulations, ncol = length(t_target_samples))

# Loop through posterior samples and generate new response times
for (i in 1:N_simulations) {
  for (j in 1:length(t_target_samples)) {
    # Get the posterior samples for the current iteration
    sigma_0_sample <- sigma_0_samples[i]
    k_sample <- k_samples[i]
    threshold_sample <- threshold_samples[i]
    t_target_sample <- t_target_samples[j]
    
    # Get the predicted response time
    simulated_data[i, j] <- get_t_predicted(t_target_sample, sigma_0_sample, k_sample, threshold_sample)
  }
}

# Plot the observed data and simulated data
library(ggplot2)

# Convert simulated data to a long format for ggplot
simulated_long <- data.frame(
  time = as.vector(simulated_data),
  target_time = data$known_t_to_target / data$block_duration,
  type = "Simulated"
)

# Add the observed data to the dataframe
observed_data <- data.frame(
  time = data$time_since_last_cc / data$block_duration,
  target_time = data$known_t_to_target / data$block_duration,
  type = "Observed"
)

# Combine observed and simulated data
combined_data <- rbind(
  simulated_long,
  observed_data
)

# Plot the comparison
combined_data %>% 
  # filter(target_time == 1) %>%
  ggplot(aes(x = time)) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity", fill = "blue") +
  geom_density(aes(color = type), size = 1) +
  labs(title = "Posterior Predictive Check: Observed vs Simulated",
       x = "Time",
       y = "Density") +
  scale_color_manual(values = c("Observed" = "red", "Simulated" = "blue")) +
  theme_minimal()