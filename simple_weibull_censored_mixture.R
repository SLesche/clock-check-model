library(dplyr)
library(rstan)
library(survival)
library(ggsurvfit)


simulate_mixture_weibull <- function(n, g, k1, lambda1, lambda2){
  ifelse(runif(n, 0, 1) <= g, rweibull(n, 1, lambda2), rweibull(n, k1, lambda1))
}

plot_pred_vs_actual <- function(pred_surv, actual_surv, event_status, data) {
  # Create survival objects
  pred_fit <- survfit(Surv(pred_surv, event_status) ~ 1)
  actual_fit <- survfit(Surv(actual_surv, event_status) ~ 1, data = data)
  
  # Convert to data frames for ggplot
  pred_df <- broom::tidy(pred_fit) %>% mutate(Type = "Predicted")
  actual_df <- broom::tidy(actual_fit) %>% mutate(Type = "Actual")
  
  # Combine the data
  survival_df <- bind_rows(pred_df, actual_df)
  
  # Plot both curves
  ggplot(survival_df, aes(x = time, y = estimate, color = Type)) +
    geom_step() +
    labs(
      x = "Days",
      y = "Overall survival probability",
      title = "Predicted vs. Actual Survival Curves"
    ) +
    theme_minimal()
}

simulated_data <- data.frame(
  times = simulate_mixture_weibull(10000, 0.1, 2.45, 0.5, 0.2),
  event = 1
)

stan_data <- list(
  N = nrow(simulated_data),
  times = simulated_data$times,
  event = simulated_data$event,
  K_k = 1,
  K_lambda = 1,
  K_lambda2 = 1,
  K_mixture = 1,
  X_k = model.matrix( ~1, data = simulated_data),
  X_lambda = model.matrix(~ 1, data = simulated_data),
  X_lambda2 = model.matrix(~ 1, data = simulated_data),
  X_mixture = model.matrix(~ 1, data = simulated_data)
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

# Example usage (assuming clean_data exists)
plot_pred_vs_actual(
  pred_surv = get_mixture_prediction(nrow(clean_data), 0.2, 2.45, 0.65, 1), 
  actual_surv = clean_data$r_check,
  event_status = clean_data$event,
  data = clean_data
)


survfit2(Surv(get_mixture_prediction(1000, 0.1, 2.45, 0.65, 1)) ~ 1) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )



# Kaplan-Meier Plot
survfit2(Surv(r_check, event) ~ 1, data = clean_data) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )


