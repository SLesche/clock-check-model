library(dplyr)
library(rstan)
library(model.matrix)
library(survival)

data <- read.csv("weibull_data.csv")

clean_data <- data %>%
  # filter(
  #   block_duration == known_t_to_target
  # ) %>%
  mutate(
    r_check = time_since_last_cc / known_t_to_target,
    r_to_target = known_t_to_target / 300,
  ) %>% # filter(r > 1) %>% View()
  mutate(
    censor_reason = ifelse(cens == 1 & accessed_pm == 1, "waited_for_pm", 
                           ifelse(cens == 1 & accessed_pm == 0, "clock_ran_out", "no_censor"))
  ) %>% 
  # mutate(
  #   r_check = ifelse(cens == 1, 1, r_check)
  # ) %>% 
  # filter(r_check < 5, r_check > 0) %>%
  mutate(
    event = ifelse(censor_reason == "clock_ran_out", 1 - cens, 1),
  ) %>% 
  ungroup()

data_censored <- clean_data %>% 
  filter(event == 0)
data_uncensored <- clean_data %>% 
  filter(event == 1)

hist(data_uncensored$time_since_last_cc, breaks = 50)
hist(data_censored$time_since_last_cc, breaks = 50)

stan_data <- list(
  N_censored = nrow(data_censored),
  N_uncensored = nrow(data_uncensored),
  censored_times = data_censored$time_since_last_cc,
  uncensored_times = data_uncensored$time_since_last_cc,
  K_k = 2,
  K_lambda = 2,
  X_k_censored = model.matrix( ~ known_t_to_target, data = data_censored),
  X_lambda_censored = model.matrix(~ known_t_to_target, data = data_censored),
  X_k_uncensored = model.matrix( ~ known_t_to_target, data = data_uncensored),
  X_lambda_uncensored = model.matrix(~ known_t_to_target, data = data_uncensored)
)

# Fit the model
options(mc.cores = parallel::detectCores())
fit <- stan(
  file = "simple_weibull_censored.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
  # control = list(adapt_delta = 0.95, max_treedepth = 15, stepsize = 0.01)
)

traceplot(fit)
print(fit)

