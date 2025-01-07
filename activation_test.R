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