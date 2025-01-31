library(dplyr)
library(rstan)
library(model.matrix)

data <- read.csv("weibull_data.csv")

clean_data <- data %>% 
  # filter(
  #   block_duration == known_t_to_target
  # ) %>%
  mutate(
    r_check = time_since_last_cc / known_t_to_target,
    r_to_target = known_t_to_target / block_duration,
  ) %>% # filter(r > 1) %>% View()
  filter(r_check < 1.2) %>% 
  filter(cens == 0)

stan_data <- list(
  N = nrow(clean_data),  # Number of events per participant
  clock_check_time = clean_data$r_check,
  K_k = 2,
  K_lambda = 2,
  X_k = model.matrix( ~ r_to_target, data = clean_data),
  X_lambda = model.matrix(~ r_to_target, data = clean_data)
)

# Fit the model
options(mc.cores = parallel::detectCores())
fit <- stan(
  file = "simple_weibull_timedep.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
  # control = list(adapt_delta = 0.95, max_treedepth = 15, stepsize = 0.01)
)

traceplot(fit)
print(fit)

# Assuming 'fit' is the result from sampling the Stan model
# Extract the posterior samples for alpha_subject and sigma_subject
posterior_samples <- extract(fit)
library(rstan)
library(loo)
library(bayesplot)

# Function to generate posterior predictive samples from Weibull regression
generate_pp_samples <- function(stan_fit, stan_data, num_samples = 1000) {
  # Extract posterior draws for k and lambda
  posterior_samples <- rstan::extract(stan_fit)
  
  # Compute k and lambda for each posterior draw
  k_pred <- exp(stan_data$X_k %*% t(posterior_samples$k))  # Matrix multiplication
  lambda_pred <- exp(stan_data$X_lambda %*% t(posterior_samples$lambda))
  
  # Generate posterior predictive samples
  y_rep <- matrix(nrow = nrow(k_pred), ncol = num_samples)
  for (i in 1:num_samples) {
    y_rep[, i] <- rweibull(nrow(k_pred), shape = k_pred[, i], scale = lambda_pred[, i])
  }
  
  return(y_rep)
}

# Function for model comparison using WAIC & LOO
compute_model_criteria <- function(stan_fit) {
  log_lik <- extract_log_lik(stan_fit)  # Extract log-likelihood
  loo_result <- loo(log_lik)
  waic_result <- waic(log_lik)
  
  return(list(loo = loo_result, waic = waic_result))
}

# Function to visualize posterior predictive checks
plot_pp_checks <- function(y_obs, y_rep) {
  ppc_dens_overlay(y_obs, y_rep[1:100, ])  # Overlay density plot
}

y_rep <- generate_pp_samples(fit, stan_data, 500)

plot_pp_checks(stan_data$clock_check_time, t(y_rep))

# Example usage:
# fit <- stan("your_model.stan", data = stan_data, iter = 2000, chains = 4)
# y_rep <- generate_pp_samples(fit, stan_data)
# model_criteria <- compute_model_criteria(fit)
# plot_pp_checks(stan_data$clock_check_time, y_rep)
