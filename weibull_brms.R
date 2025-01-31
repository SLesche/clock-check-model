library(dplyr)
library(brms)

data <- read.csv("archive/diffusion_data.csv")

clean_data <- data %>% 
  # filter(
  #   block_duration == known_t_to_target
  # ) %>%
  mutate(
    r = time_since_last_cc / known_t_to_target
  ) %>% 
  filter(r > 0.01, known_t_to_target > 0) %>% 
  arrange(participant) %>% # make sure that the data is not too close to 0 (this cause issues in integration)
  mutate(subject_id = dense_rank(participant)) %>% 
  mutate(known_t_to_target)

weibull_model <- brm(
  bf(time_since_last_cc ~ 1 + known_t_to_target, shape ~ 1 + (1 | participant)),  # Weibull AFT model
  family = weibull(),
  data = clean_data,
  chains = 4, iter = 2000, warmup = 1000, cores = 4,
  prior = c(
    prior(normal(0, 2), class = "b"),  # Prior for regression coefficient
    prior(normal(0, 2), class = "Intercept") # Prior for Intercept
    # prior(normal(1, 1), class = "shape")
    )  # Prior for shape parameter
)

summary(weibull_model)

pp_check(weibull_model)

posterior_samples(weibull_model) %>% head()

plot(conditional_effects(weibull_model, "x"))
