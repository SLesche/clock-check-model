library(dplyr)
library(rstan)
library(survival)
library(ggsurvfit)

# Functions ----
inv_logit <- function(x){
  exp(x) / (1 + exp(x))
}

simulate_mixture_weibull <- function(n, g, k1, lambda1, lambda2 = 1){
  ifelse(runif(n, 0, 1) <= g, rweibull(n, 1, lambda2), rweibull(n, k1, lambda1))
}

simulate_from_predicted <- function(data, fit, ndraws){
  result = matrix(NA, nrow = nrow(data), ncol = ndraws)
  
  posteriors = extract(fit, pars = c("k_pred", "lambda_pred", "mixture_prob"))
  
  sample_ids = sample(1:nrow(posteriors$k_pred), ndraws)
  posterior_samples_k = posteriors$k_pred[sample_ids, ]
  posterior_samples_lambda = posteriors$lambda_pred[sample_ids, ]
  # posterior_samples_lambda2 = 1
  posterior_samples_mixture = posteriors$mixture_prob[sample_ids, ]
  
  for (idraw in 1:ndraws){
    result[, idraw] = simulate_mixture_weibull(
      nrow(data), 
      posterior_samples_mixture[idraw, ],
      posterior_samples_k[idraw, ],
      posterior_samples_lambda[idraw, ],
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

# Data ----
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
  mutate(
    subject_id = dense_rank(participant)
  ) %>% 
  # filter(event == 1) %>% 
  filter(!(is_first_guess == 1 & r_check < 0.1)) %>% 
  # filter(is_first_guess == 0) %>%
  ungroup() 


# Plot survival curves
surv_data <- clean_data %>% 
  mutate(
    times_to_target = case_when(
      known_t_to_target > 250 ~ "long",
      known_t_to_target <= 250 & known_t_to_target > 100 ~ "semi-long",
      known_t_to_target <= 100 & known_t_to_target > 50 ~ "mid",
      known_t_to_target <= 50 ~ "short"
    )
  )

ggsurvfit(survfit(Surv(r_check, event) ~ times_to_target, data = surv_data)) 


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

pred_matrix <- model.matrix( ~ 1, data = clean_data)

stan_data <- list(
  N = nrow(clean_data),
  times = clean_data$r_check,
  event = clean_data$event,
  K_k = ncol(pred_matrix),
  K_lambda = ncol(pred_matrix),
  K_mixture = ncol(pred_matrix),
  X_k = pred_matrix,
  X_lambda = pred_matrix,
  X_mixture = pred_matrix
)

# Fit the model
options(mc.cores = parallel::detectCores())
fit <- stan(
  file = "very_simple_weibull_censored_mixed_timedep.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
  # control = list(adapt_delta = 0.95, max_treedepth = 15, stepsize = 0.01)
)

traceplot(fit, pars = c("k", "lambda", "mixture_rate"))

bayesplot::mcmc_areas(
  fit,
  pars = c("k[1]", "lambda[1]", "mixture_rate[1]"),
  prob = 0.9
)

inv_logit(-4)

posteriors <- extract(fit, pars = c("mixture_prob"))

lm(colMeans(posteriors$mixture_prob) ~ clean_data$r_to_target) %>% summary()

predicted_data <- simulate_from_predicted(clean_data, fit, 100)

plot_pred_vs_actual(
  predicted_data,
  clean_data$r_check,
  clean_data$event,
  clean_data
)+ 
  xlim(0, 2)

## ---- Hierarchical

pred_matrix <- model.matrix( ~ r_to_target, data = clean_data)

stan_data <- list(
  N = nrow(clean_data),
  times = clean_data$r_check,
  event = clean_data$event,
  K_k = ncol(pred_matrix),
  K_lambda = ncol(pred_matrix),
  # K_lambda2 = 1,
  K_mixture = ncol(pred_matrix),
  J = length(unique(clean_data$subject_id)),
  group = clean_data$subject_id,
  X_k = pred_matrix,
  X_lambda = pred_matrix,
  # X_lambda2 = model.matrix(~ is_first_guess, data = clean_data),
  X_mixture = pred_matrix
)

# Fit the model
options(mc.cores = parallel::detectCores())
hierarch_fit <- stan(
  file = "very_simple_hierarch_weibull_censored_mixed_timedep.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
  # control = list(adapt_delta = 0.95, max_treedepth = 15, stepsize = 0.01)
)

subject_id <- unique(clean_data$subject_id)
subject_k <- colMeans(as.data.frame(extract(hierarch_fit, "k_raw")))
# subject_k_eff <- colMeans(as.data.frame(extract(hierarch_fit, "k_raw")))

subject_lambda <- colMeans(as.data.frame(extract(hierarch_fit, "lambda_raw")))
# subject_lambda_eff <- colMeans(as.data.frame(extract(hierarch_fit, "lambda_raw")))
subject_g <- colMeans(as.data.frame(extract(hierarch_fit, "mixture_rate_raw")))
# subject_g_eff <- colMeans(as.data.frame(extract(hierarch_fit, "mixture_rate_raw")))

cc_metrics <- data.frame(
  subject_id = subject_id,
  k = subject_k[1:length(subject_id)],
  lambda = subject_lambda[1:length(subject_id)],
  g = subject_g[1:length(subject_id)],
  k_effect = subject_k[(length(subject_id) + 1):(2*length(subject_id))],
  lambda_effect = subject_lambda[(length(subject_id) + 1):(2*length(subject_id))],
  g_effect = subject_g[(length(subject_id) + 1):(2*length(subject_id))]
)

# traceplot(hierarch_fit, pars = c("k_raw", "lambda_raw", "mixture_rate_raw"))

# bayesplot::mcmc_areas(
#   hierarch_fit,
#   pars = c("mu_mixture[1]", "mu_mixture[2]"),
#   prob = 0.9
# )

# posteriors <- extract(hierarch_fit, pars = c(""))

predicted_data <- simulate_from_predicted(clean_data, hierarch_fit, 100)

plot_pred_vs_actual(
  predicted_data,
  clean_data$r_check,
  clean_data$event,
  clean_data
)+
  xlim(0, 2)

# Other metrics ----
mean_times <- clean_data %>% 
  group_by(subject_id) %>% 
  summarize(
    mean_time_between_ccs = mean(time_since_last_cc, na.rm = TRUE),
    mean_r_check = mean(r_check, na.rm = TRUE),
    n_checks = sum(cens == 0),
    pm_acc = mean(unique(pm_acc)),
    pm_count = unique(pm_count)
  )

# Mean time between clock checks

# Compute strategic clock checking
strat_check <- clean_data %>% 
  mutate(
    is_last_quarter = ifelse(time_since_start > 300 * 3/4, 1, 0)
  ) %>% 
  count(subject_id, is_last_quarter) %>% 
  group_by(subject_id) %>% 
  mutate(
    strat_check = round(100*n / sum(n), 2)
  ) %>% 
  rename(
    "n_checks_last_quarter" = n
  ) %>% 
  filter(is_last_quarter == 1)

cc_compare <- cc_metrics %>% 
  left_join(., mean_times) %>% 
  left_join(., strat_check) 

cc_compare %>% 
  select(-subject_id, -is_last_quarter) %>% 
  cor() %>% View()
