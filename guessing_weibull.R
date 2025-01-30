hazard_function <- function(t, b, k, g){
  g + (1 - g)*b*k*t^(k-1)
}

cumulative_hazard <- function(t, b, k, g){
  g*t + (1 - g)  *b*t^k
}

pdf <- function(t, b, k, g){
  hazard_function(t, b, k, g) * exp(-cumulative_hazard(t, b, k, g))
}

times <- seq(0, 1, 0.01)

plot(times, hazard_function(times, 3, 3, 0))
plot(times, cumulative_hazard(times, 3, 3, 0))
plot(times, pdf(times, 3, 2, 0.65))


data <- read.csv("archive/diffusion_data.csv")

fit_data <- data %>% 
  filter(known_t_to_target != 0, time_since_last_cc != 0) %>% 
  filter(
  block_duration == known_t_to_target
  ) %>%
  mutate(
    r = time_since_last_cc / known_t_to_target
  ) %>% 
  filter(r > 0)

stan_data <- list(
  N = nrow(fit_data),  # Number of events per participant
  clock_check_time = fit_data$r,
  known_t_to_target = fit_data$known_t_to_target
)

# Fit the model
options(mc.cores = parallel::detectCores())
fit <- stan(
  file = "weibull_guessing.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
)

traceplot(fit)
print(fit)

pdf <- function(t, b, k, g){
  hazard_function(t, b, k, g) * exp(-cumulative_hazard(t, b, k, g))
}

plot(times, pdf(times, 0.52, 1.95, 0))
lines(density(fit_data$r))
