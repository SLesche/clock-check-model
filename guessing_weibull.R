library(ggplot2)
library(dplyr)
library(rstan)
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
  # filter(
  # block_duration == known_t_to_target
  # ) %>%
  mutate(
    r = time_since_last_cc / known_t_to_target
  ) %>% 
  filter(r > 0)

results <- data.frame(
  subject = unique(fit_data$participant),
  g = numeric(length(unique(fit_data$participant))),
  b = numeric(length(unique(fit_data$participant))),
  k = numeric(length(unique(fit_data$participant)))
)

for (isubject in 1:nrow(results)){
  cc_data = fit_data %>% filter(participant == results[isubject, "subject"])
  fit = stan(
    file = "weibull_guessing.stan",
    data = list(N = nrow(cc_data), clock_check_time = cc_data$r),
    iter = 2000,  # Number of iterations
    chains = 4,   # Number of MCMC chains
    warmup = 500, # Number of warmup iterations
    # control = list(adapt_delta = 0.95, max_treedepth = 15, stepsize = 0.01)
  )
  
  g = mean(extract(fit)$g)
  b = mean(extract(fit)$b)
  k = mean(extract(fit)$k)
  
  results[isubject, "g"] = g
  results[isubject, "k"] = k
  results[isubject, "b"] = b
}

predicted_pdf <- matrix(NA, nrow = nrow(results), ncol = length(times))
for (isubject in 1:nrow(results)){
  subj_pdf = pdf(times, results[isubject, "b"], results[isubject, "k"], results[isubject, "g"])
  
  predicted_pdf[isubject, ] = subj_pdf
}

mean_pdf <- colMeans(predicted_pdf, na.rm = TRUE)
plot(times, mean_pdf)
lines(density(fit_data$r))

full_data <- read.csv("data_clean.csv") 

data <- results %>% left_join(., full_data %>% select(subject = participant, pm_count) %>% distinct(subject, pm_count))
data$cc_count <- clean_data %>% count(subject_id) %>% pull(n)

cor(data[, 2:5], use = "pairwise.complete.obs")

full_data %>% 
  distinct(participant, pm_count)

stan_data <- list(
  N = nrow(fit_data),  # Number of events per participant
  clock_check_time = fit_data$r,
  known_t_to_target = (fit_data$known_t_to_target / fit_data$block_duration),
  J = length(unique(fit_data$participant)),
  subject = fit_data$participant
)

# Fit the model
options(mc.cores = parallel::detectCores())
fit <- stan(
  file = "weibull_guessing_hierarch.stan",
  data = stan_data,
  iter = 4000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 1000, # Number of warmup iterations
  # control = list(adapt_delta = 0.95, max_treedepth = 15, stepsize = 0.01)
)

traceplot(fit)
print(fit)

pdf <- function(t, b, k, g){
  hazard_function(t, b, k, g) * exp(-cumulative_hazard(t, b, k, g))
}

plot(times, pdf(times, 7.15, 3.78, 0.2))
lines(density(fit_data$r))

# Extract posterior samples
posterior_samples <- extract(fit)

