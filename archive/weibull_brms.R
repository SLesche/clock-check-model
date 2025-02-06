library(dplyr)
library(brms)

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
  mutate(
    event = ifelse(censor_reason == "clock_ran_out", 1 - cens, 1),
  ) %>% 
  # filter(
  #   accessed_pm == 1
  # ) %>%
  ungroup()

data_censored <- clean_data %>% 
  filter(event == 0)
data_uncensored <- clean_data %>% 
  filter(event == 1)

hist(data_uncensored$r_check, breaks = 50)
hist(data_censored$r_check, breaks = 50)

weibull_model <- brm(
  bf(r_check | cens(event) ~ r_to_target + (r_to_target | subject_id)),  # Weibull ATF model
  family = weibull(),
  data = clean_data,
  chains = 4, iter = 2000, warmup = 1000, cores = 4,
  # prior = c(
  #   prior(normal(0, 2), class = "b"),  # Prior for regression coefficient
  #   prior(normal(0, 2), class = "Intercept") # Prior for Intercept
  #   # prior(normal(1, 1), class = "shape")
  #   )  # Prior for shape parameter
)

summary(weibull_model)

pp_check(weibull_model, ndraws = 50)

posterior_samples(weibull_model) %>% head()

plot(conditional_effects(weibull_model, "x"))
