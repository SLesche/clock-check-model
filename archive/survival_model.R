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

