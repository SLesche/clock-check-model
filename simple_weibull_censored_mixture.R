library(dplyr)
library(rstan)
library(survival)
library(ggsurvfit)

inv_logit <- function(x){
  exp(x) / (1 + exp(x))
}

simulate_mixture_weibull <- function(n, g, k1, lambda1, lambda2){
  ifelse(runif(n, 0, 1) <= g, rweibull(n, 1, lambda2), rweibull(n, k1, lambda1))
}

simulate_from_predicted <- function(data, fit, ndraws){
  result = matrix(NA, nrow = nrow(data), ncol = ndraws)
  
  posteriors = extract(fit, pars = c("k_pred", "lambda_pred", "lambda2_pred", "mixture_prob"))
  
  sample_ids = sample(1:nrow(posteriors$k_pred), ndraws)
  posterior_samples_k = posteriors$k_pred[sample_ids, ]
  posterior_samples_lambda = posteriors$lambda_pred[sample_ids, ]
  posterior_samples_lambda2 = posteriors$lambda2_pred[sample_ids, ]
  posterior_samples_mixture = posteriors$mixture_prob[sample_ids, ]
  
  for (idraw in 1:ndraws){
    result[, idraw] = simulate_mixture_weibull(
      nrow(data), 
      posterior_samples_mixture[idraw, ],
      posterior_samples_k[idraw, ],
      posterior_samples_lambda[idraw, ],
      posterior_samples_lambda2[idraw, ]
      )
  }
  
  return(result)
}

plot_pred_vs_actual <- function(pred_surv, actual_surv, event_status, data) {
  # Ensure pred_surv is a matrix
  if (!is.matrix(pred_surv)) {
    pred_surv <- matrix(pred_surv, ncol = 1)  # Convert vector to single-column matrix
  }
  
  # Create actual survival fit
  actual_fit <- survfit(Surv(actual_surv, event_status) ~ 1, data = data)
  actual_df <- broom::tidy(actual_fit) %>% mutate(Type = "Actual", idraw = 0)
  
  # Create predicted survival fits for each column of pred_surv
  pred_dfs <- lapply(1:ncol(pred_surv), function(i) {
    pred_fit <- survfit(Surv(pred_surv[, i], event_status) ~ 1)
    broom::tidy(pred_fit) %>% mutate(Type = "Predicted", idraw = i)  # Assign idraw for grouping
  })
  
  # Combine all predicted survival curves
  pred_df <- bind_rows(pred_dfs)
  
  # Combine actual and predicted data
  survival_df <- bind_rows(pred_df, actual_df)
  
  # Plot survival curves
  ggplot(survival_df, aes(x = time, y = estimate, group = idraw)) +
    geom_step(data = pred_df, aes(color = "Predicted"), alpha = 0.4, size = 0.5) +  # Blue, low alpha
    geom_step(data = actual_df, aes(color = "Actual"), size = 1.2) +  # Red, thicker line
    scale_color_manual(values = c("Actual" = "#FF474C", "Predicted" = "lightblue")) +
    labs(
      x = "Relative Check Time",
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
  X_k = model.matrix( ~ 1, data = simulated_data),
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

# returns k = exp(0.9), lambda = exp(-0.7), lambda2 = exp(-1.53), g = exp(-2.19) / (1 + exp(-2.19))
plot_pred_vs_actual(
  simulate_mixture_weibull(nrow(simulated_data), exp(-2.34) / (1 + exp(-2.34)), exp(0.89), exp(-0.7), exp(-1.70)),
  simulated_data$times,
  simulated_data$event,
  simulated_data
)

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
  K_k = 2,
  K_lambda = 2,
  K_lambda2 = 2,
  K_mixture = 2,
  X_k = model.matrix( ~ is_first_guess, data = clean_data),
  X_lambda = model.matrix(~ is_first_guess, data = clean_data),
  X_lambda2 = model.matrix(~ is_first_guess, data = clean_data),
  X_mixture = model.matrix(~ is_first_guess, data = clean_data)
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

traceplot(fit, pars = c("k[1]", "lambda[1]", "lambda2[1]", "mixture_rate[1]"))

bayesplot::mcmc_areas(
  fit,
  pars = c("k[1]", "k[2]", "lambda[1]", "lambda[2]"),
  prob = 0.9
)

posteriors <- extract(fit, pars = c("k[1]", "lambda[1]", "lambda2[1]", "mixture_rate[1]"))

predicted_data <- simulate_from_predicted(clean_data, fit, 100)

plot_pred_vs_actual(
  predicted_data,
  clean_data$r_check,
  clean_data$event,
  clean_data
)


## ---- Hierarchical
stan_data <- list(
  N = nrow(clean_data),
  times = clean_data$r_check,
  event = clean_data$event,
  K_k = 1,
  K_lambda = 1,
  K_lambda2 = 1,
  K_mixture = 1,
  S = length(unique(clean_data$subject_id)),
  subject_id = clean_data$subject_id,
  X_k = model.matrix( ~ 1, data = clean_data),
  X_lambda = model.matrix(~ 1, data = clean_data),
  X_lambda2 = model.matrix(~ 1, data = clean_data),
  X_mixture = model.matrix(~ 1, data = clean_data)
)

# Fit the model
options(mc.cores = parallel::detectCores())
hierarch_fit <- stan(
  file = "hierarch_weibull_censored_mixed_timedep.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
  # control = list(adapt_delta = 0.95, max_treedepth = 15, stepsize = 0.01)
)

traceplot(hierarch_fit, pars = c("k[1]", "lambda[1]", "lambda2[1]", "mixture_rate[1]"))

bayesplot::mcmc_areas(
  hierarch_fit,
  pars = c("k[1]", "k[2]", "k[3]"),
  prob = 0.9
)

posteriors <- extract(hierarch_fit, pars = c("k[1]", "lambda[1]", "lambda2[1]", "mixture_rate[1]"))

predicted_data <- simulate_from_predicted(clean_data, hierarch_fit, 100)

plot_pred_vs_actual(
  predicted_data,
  clean_data$r_check,
  clean_data$event,
  clean_data
)






