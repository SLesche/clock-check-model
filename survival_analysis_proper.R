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
    r_to_target = known_t_to_target / block_duration,
  ) %>% # filter(r > 1) %>% View()
  mutate(event = 1 - cens) %>% 
  # filter(r_check < 1.2, r_check > 0) %>%
  # filter(accessed_pm == 1) %>% 
  mutate(subject_id = dense_rank(subject_id))

# Get mean time between clock_checks
cc_metrics <- clean_data %>% 
  group_by(subject_id, block_num, accessed_pm) %>% 
  summarize(
    mean_time_between_ccs = mean(time_since_last_cc),
    time_last_cc = max(time_since_start),
    n_ccs = n(),
    strat_cc = sum(time_since_start > 300 * 0.75) / n()
  )

cc_metrics %>% 
  glm(accessed_pm ~ mean_time_between_ccs + time_last_cc + n_ccs + strat_cc, data = ., family =binomial(link = "logit")) %>% 
  summary()

# Kaplan-Meier Plot
survfit2(Surv(r_check, event) ~ time_till_pm, data = clean_data %>% filter(r_check > 0, r_check < 1.5) %>% mutate(time_till_pm = paste(floor(known_t_to_target / 60)))) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )
