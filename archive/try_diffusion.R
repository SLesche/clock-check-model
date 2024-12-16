library(rstan)
library(dplyr)

get_t_predicted <- function(t_target, sigma_0, k, threshold){
  t_pred = t_target / (1 + qnorm(1 - threshold) * sigma_0 * k)
  
  return(t_pred)
}
get_pred_ratio <- function(sigma, k, threshold){
  ratio = 1 / (qnorm(1-threshold) * sigma * k + 1)
  return(ratio)
}

simulate_blockwise_behavior <- function(t_total, sigma_0, k, threshold){
  check_times = c()
  t_start = 0
  delta = 30
  while (t_start < t_total - delta) {
    t_pred = get_t_predicted(t_total - t_start, sigma_0, k, threshold)
    t_check = t_start + t_pred
    check_times = c(check_times, t_check)
    
    t_start = t_check
  }
  
  return(check_times)
}

data <- read.csv("diffusion_data.csv")

stan_data <- list(
  Nsubj = length(unique(data$participant)),  # Number of participants
  Nblocks = max(unique(data$block_num)) + 1,
  N = nrow(data),  # Number of events per participant
  person_id = data$participant,
  P = length(unique(data$participant)),
  block_id = data$block_num + 1,
  observed_time = data$time_since_last_cc / data$block_duration,
  known_t_to_target = data$known_t_to_target / data$block_duration
)

# Fit the model
fit <- stan(
  file = "diffusion_model_cc_hierarch.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 1,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
)

# R-hat statistics
print(fit)

# Trace plots for parameters
traceplot(fit)

# Extract posterior samples from the fitted model
posterior_samples <- extract(fit)

# Posterior samples for sigma_0, k, and threshold
sigma_0_samples <- posterior_samples$sigma_0
k_samples <- posterior_samples$k
threshold_samples <- posterior_samples$theta
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

