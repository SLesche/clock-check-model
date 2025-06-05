library(survival)
library(ggplot2)
library(dplyr)
library(rstan)
library(brms)

data <- read.csv("archive/diffusion_data.csv")

fit_data <- data %>% 
  filter(known_t_to_target != 0, time_since_last_cc != 0) %>% 
  # filter(
  # block_duration == known_t_to_target
  # ) %>%
  mutate(
    r = time_since_last_cc / known_t_to_target
  ) %>% 
  filter(r > 0)

stan_data <- list(
  N = nrow(fit_data),  # Number of events per participant
  clock_check_time = fit_data$r
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

samples <- extract(fit)

# Example posterior samples (from Stan model)
lambda_0_post <- samples$lambda_0  # Posterior samples of λ₀
eta_post <- samples$eta            # Posterior samples of η
beta_post <- samples$beta          # Posterior samples of β

# Compute predicted CDF for a grid of time points
time_grid <- seq(0, max(failure_times), length.out = 100)
predicted_cdfs <- sapply(1:length(lambda_0_post), function(i) {
  lambda_0 <- lambda_0_post[i]
  eta <- eta_post[i]
  beta <- beta_post[i]
  M_t <- lambda_0 * time_grid + (eta / beta) * time_grid^beta
  1 - exp(-M_t)  # Convert cumulative intensity to predicted CDF
})

# Compute the mean predicted CDF (posterior predictive mean)
predicted_cdf_mean <- rowMeans(predicted_cdfs)

# Observed failure times (sorted)
empirical_cdf <- ecdf(failure_times)  # Creates an empirical CDF function

# Plot empirical CDF
plot(empirical_cdf, main = "Empirical vs. Predicted CDF",
     xlab = "Time", ylab = "CDF", col = "black", lwd = 2)

# Add predicted CDF
lines(time_grid, predicted_cdf_mean, col = "blue", lwd = 2)
legend("bottomright", legend = c("Empirical CDF", "Predicted CDF"),
       col = c("black", "blue"), lwd = 2)



plot(density(fit_data$r))
summary(lm(r ~ known_t_to_target, data = fit_data))

failure_times <- fit_data$r

n <- length(failure_times)
ranked_times <- sort(failure_times)
cumulative_prob <- (1:n - 0.3) / (n + 0.4)  #  Empirical cumulative distribution function
# cumulative_prob <- (1:n) / (n + 1)  # Median rank
log_time <- log(ranked_times)
log_minus_log <- log(-log(1 - cumulative_prob))

plot_data <- data.frame(
  log_time = log_time,
  log_minus_log = log_minus_log
)

ggplot(plot_data, aes(x = log_time, y = log_minus_log)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    title = "Weibull Probability Plot",
    x = "Log(Failure Time)",
    y = "Log(-Log(1 - Cumulative Probability))"
  ) +
  theme_minimal()
