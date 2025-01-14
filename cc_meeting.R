library(dplyr)
library(purrr)

prob_target <- function(t_target, t_step, k){
  return(1 - pnorm(t_target, t_step, k*t_step))
}

prob_action <- function(prob_target, c, g, a = 10){
  return(g + (1 - g)/(1 + exp(-a*(log(prob_target / (1 - prob_target)) + c))))
}

plot_action_prob <- function(t_target, nsteps, g, c, k, a){
  t_steps <- seq(0, t_target, t_target/nsteps)
  probability_action <- prob_action(prob_target(t_target, t_steps, k), c * t_target, g, a)
  # prob_no_action_yet <- prob_no_action_until_time(t_target, k, g, threshold, a, nsteps)
  plot(t_steps, probability_action)
}
# Hazard function
hazard_function <- function(t, t_target, k, g, a, c) {
  z <- (1 - pnorm(t_target, t, k * t)) / pnorm(t_target, t, k * t)
  action_prob <- g + (1 - g) / (1 + exp(-a * c) * z^-a)
  return(action_prob)
}

# Survival function
survival_function <- function(t, t_target, k, g, a, c) {
  # Integrate the hazard function from 0 to t
  cumulative_hazard <- integrate(function(u) hazard_function(u, t_target, k, g, a, c), 
                                 lower = 0, upper = t)$value
  return(exp(-cumulative_hazard))
}

# PDF of action times
pdf_action_time <- function(t, t_target, k, g, a, c) {
  lambda_t <- hazard_function(t, t_target, k, g, a, c)
  S_t <- survival_function(t, t_target, k, g, a, c)
  return(lambda_t * S_t)
}

# Parameters
t_target <- 1    # Example target time
k <- 2    # Example scale parameter
g <- 0.1     # Example baseline probability
a <- 0.1   # Example slope parameter
c <- 9# Example scaling factor
n_steps <- 100
t_steps <- seq(0, t_target, t_target/n_steps)
plot(t_steps, purrr::map_dbl(t_steps, ~hazard_function(., t_target, k, g, a,c)))
plot(t_steps, purrr::map_dbl(t_steps, ~survival_function(., t_target, k, g, a, c)))
plot(t_steps, purrr::map_dbl(t_steps, ~pdf_action_time(., t_target, k, g, a, c)))


data <- read.csv("archive/diffusion_data.csv")

fit_data <- data %>% 
  filter(known_t_to_target != 0, time_since_last_cc != 0)

# Plot of first clock_check
fit_data %>% filter(known_t_to_target == block_duration) %>% pull(time_since_start) %>% hist()

times <- seq(0, max(fit_data$block_duration), 0.5)

surv_per <- purrr::map_dbl(times, \(x) sum(first_checks$time_since_start > x) / length(first_checks$time_since_start))

plot(times, surv_per)

checks <- fit_data %>% 
  mutate(
    rate_not_checked = time_since_last_cc / known_t_to_target
  ) 

checks %>%filter(known_t_to_target < 300) %>%  pull(rate_not_checked) %>% hist(breaks = 50)
