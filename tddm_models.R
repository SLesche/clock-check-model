library(rstan)

# Simulated data
N <- 1000
t_target <- runif(N, 0.5, 100)
theta <- 0.8
sigma <- 0.3
z <- theta * t_target
mu <- z
lambda <- z^2 / sigma^2
rt <- statmod::rinvgauss(N, mean = mu, shape = lambda)

stan_data <- list(
  N = N,
  rt = rt,
  t_target = t_target
)

# Fit the model
fit <- stan(
  file = "tddm_model.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  cores = 4
)

print(fit)


library(survival)
library(ggplot2)
library(dplyr)
library(rstan)
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


stan_data <- list(
  N = nrow(clean_data),  # Number of events per participant
  t_target = clean_data$known_t_to_target,
  rt = clean_data$time_since_last_cc,
  J = length(unique(clean_data$subject_id)),
  subj = clean_data$subject_id
)

# Fit the model
options(mc.cores = parallel::detectCores())

# Fit the model
fit <- stan(
  file = "mixture_model_ddm_hierarch.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  cores = 4
)
print(fit)


ppc_ddm_ig <- function(fit, observed_rt, t_target, ndraws = 100) {
  # Extract posterior draws
  posterior <- rstan::extract(fit)
  n <- length(observed_rt)
  draws <- sample(1:length(posterior$theta), ndraws)
  
  sim_data <- matrix(NA, nrow = ndraws, ncol = n)
  
  for (i in 1:ndraws) {
    theta_i <- posterior$theta[draws[i]]
    sigma_i <- posterior$sigma[draws[i]]
    z_i <- theta_i * t_target
    mu_i <- z_i
    lambda_i <- z_i^2 / sigma_i^2
    
    sim_data[i, ] <- statmod::rinvgauss(n, mean = mu_i, shape = lambda_i)
  }
  
  # Plot posterior predictive densities
  bayesplot::ppc_dens_overlay(y = observed_rt, yrep = sim_data[1:20, ])
}

ppc_ddm_mixture <- function(fit, observed_rt, t_target, ndraws = 100) {
  library(statmod)   # for rinvgauss
  library(bayesplot)
  
  # Extract posterior draws
  posterior <- rstan::extract(fit)
  N <- length(observed_rt)
  draws <- sample(1:length(posterior$theta), ndraws)
  
  # Prepare storage for simulated data
  sim_data <- matrix(NA, nrow = ndraws, ncol = N)
  
  # Loop over posterior draws
  for (i in seq_len(ndraws)) {
    theta_i <- posterior$theta[draws[i]]
    sigma_i <- posterior$sigma[draws[i]]
    g_i     <- posterior$g[draws[i]]
    
    z_i     <- theta_i * t_target
    mu_i    <- z_i
    lambda_i <- z_i^2 / sigma_i^2
    
    # Draw guesses or DDM samples
    is_guess <- rbinom(N, size = 1, prob = g_i)
    ddm_rt   <- statmod::rinvgauss(N, mean = mu_i, shape = lambda_i)
    guess_rt <- runif(N, min = 0, max = t_target)
    
    sim_data[i, ] <- ifelse(is_guess == 1, guess_rt, ddm_rt)
  }
  
  # Plot posterior predictive densities
  bayesplot::ppc_dens_overlay(y = observed_rt, yrep = sim_data[1:min(50, ndraws), ])
}

ppc_ddm_mixture(fit,
           clean_data$time_since_last_cc,
           t_target = clean_data$known_t_to_target,
           200)


subject_id <- unique(clean_data$subject_id)
subject_theta <- colMeans(as.data.frame(extract(fit, "theta")))
subject_sigma <- colMeans(as.data.frame(extract(fit, "sigma")))
subject_g <- colMeans(as.data.frame(extract(fit, "g")))
cc_metrics <- data.frame(
  subject_id = subject_id,
  theta = subject_theta,
  sigma = subject_sigma,
  g = subject_g
)

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
