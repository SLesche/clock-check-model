library(dplyr)
library(rstan)
library(model.matrix)
library(survival)
library(ggsurvfit)

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
  filter(r_check < 2, r_check > 0) %>%
  mutate(is_first_guess = ifelse(cc_time == 0, 1, 0)) %>% 
  mutate(
    event = ifelse(censor_reason == "clock_ran_out", 1 - cens, 1),
  ) %>% 
  filter(!is.na(event)) %>% 
  # filter(r_to_target > 0.95) %>% 
  filter(
    accessed_pm == 1
  ) %>%
  ungroup()
# 
# data_censored <- clean_data %>% 
#   filter(event == 0)
# data_uncensored <- clean_data %>% 
#   filter(event == 1)
# 
# hist(data_uncensored$time_since_last_cc, breaks = 50)
# hist(data_censored$time_since_last_cc, breaks = 50)
# hist(data_uncensored$r_check, breaks = 50)
# hist(data_censored$r_check, breaks = 50)
hist(clean_data$r_check, breaks = 50)

stan_data <- list(
  N = nrow(clean_data),
  times = clean_data$r_check,
  event = clean_data$event,
  K_k = 1,
  K_lambda = 1,
  K_lambda2 = 1,
  K_mixture = 1,
  X_k = model.matrix( ~1, data = clean_data),
  X_lambda = model.matrix(~ 1, data = clean_data),
  X_lambda2 = model.matrix(~ 1, data = clean_data),
  X_mixture = model.matrix(~ 1, data = clean_data)
)

# Fit the model
options(mc.cores = parallel::detectCores())
fit <- stan(
  file = "simple_weibull_censored_mixed_timedep.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
  # control = list(adapt_delta = 0.95, max_treedepth = 15, stepsize = 0.01)
)

traceplot(fit)
print(fit)

# Kaplan-Meier Plot
survfit2(Surv(r_check, event) ~ is_first_guess, data = clean_data %>%  filter(r_check < 1.1)) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )


