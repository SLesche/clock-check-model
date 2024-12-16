library(rstan)

data <- read.csv("cleaned_data.csv")

# Set up the data
# Replace with your actual data
stan_data <- list(
  Nsubj = length(unique(data$id)),  # Number of participants
  N = nrow(data),  # Number of events per participant
  id = data$id,
  t_since_last_check = data$t_since_last_check / 300,
  known_t_to_target = data$known_t_to_target / 300,
  check = data$check
)

# Fit the model
fit <- stan(
  file = "simple_linear.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 1,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
)

# Check convergence with traceplots
traceplot(fit)

# Get a summary of the posterior distributions
print(fit)

# Posterior predictive checks (PPC)
# Extract the simulated data from the posterior
y_rep <- extract(fit)$y_rep

# Compare the observed and simulated data (e.g., plot or summary)
# Plot the observed data vs. the predicted values
observed_check <- data$check
table(observed_check + y_rep)

# Alternatively, you can use the ppc_dens_overlay function for density overlays
ppc_dens_overlay(observed_check, y_rep)

# Other diagnostics: 
# 1. Check for Rhat values and effective sample size to ensure convergence
summary(fit)$summary[, c("Rhat", "n_eff")]

# 2. Posterior parameter diagnostics (e.g., mu_sigma_0, mu_k)
mu_sigma_0_samples <- extract(fit)$mu_sigma_0
mu_k_samples <- extract(fit)$mu_k

# Visualize the posteriors of parameters
ggplot(data.frame(mu_sigma_0_samples), aes(x = mu_sigma_0_samples)) + 
  geom_density() + ggtitle("Posterior of mu_sigma_0")

ggplot(data.frame(mu_k_samples), aes(x = mu_k_samples)) + 
  geom_density() + ggtitle("Posterior of mu_k")
