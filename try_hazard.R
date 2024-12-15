library(rstan)
library(dplyr)

data <- read.csv("hazard_data.csv")

stan_data <- list(
  Nsubj = length(unique(data$participant)),  # Number of participants
  Nblocks = max(unique(data$block_num)) + 1,
  N = nrow(data),  # Number of events per participant
  id = data$participant,
  block_id = data$block_num + 1,
  check_time = data$time_since_start,
  block_length = data$block_duration
)

# Fit the model
fit <- stan(
  file = "hazard_model_hierarch.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 1,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
)
