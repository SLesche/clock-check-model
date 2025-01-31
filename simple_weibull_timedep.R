library(dplyr)
library(rstan)
library(model.matrix)

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
  mutate(subject_id = dense_rank(participant))

stan_data <- list(
  N = nrow(clean_data),  # Number of events per participant
  clock_check_time = clean_data$r,
  K_k = 2,
  K_lambda = 1,
  X_k = model.matrix( ~ known_t_to_target, data = clean_data),
  X_lambda = model.matrix(~1, data = clean_data)
)

# Fit the model
options(mc.cores = parallel::detectCores())
fit <- stan(
  file = "simple_weibull_timedep.stan",
  data = stan_data,
  iter = 2000,  # Number of iterations
  chains = 4,   # Number of MCMC chains
  warmup = 500, # Number of warmup iterations
  # control = list(adapt_delta = 0.95, max_treedepth = 15, stepsize = 0.01)
)

traceplot(fit)
print(fit)

# Assuming 'fit' is the result from sampling the Stan model
# Extract the posterior samples for alpha_subject and sigma_subject
posterior_samples <- extract(fit)

num_samples <- length(posterior_samples$alpha_subject)

# Time points to evaluate the PDF
time_points <- seq(0, max(clean_data$r), length.out = 100)

# Matrix to store the PDF values for each time point and posterior sample
pdf_matrix <- matrix(0, nrow = num_samples, ncol = length(time_points))
# Function to compute the Weibull PDF
weibull_pdf <- function(t, alpha, sigma) {
  return((alpha / sigma) * (t / sigma)^(alpha - 1) * exp(-(t / sigma)^alpha))
}
# Loop over the posterior samples and compute the Weibull PDF
for (i in 1:num_samples) {
  pdf_matrix[i, ] <- weibull_pdf(time_points, posterior_samples$alpha_subject[i], posterior_samples$sigma_subject[i])
}
expected_pdf <- apply(pdf_matrix, 2, mean)

plot(time_points, expected_pdf)
lines(density(clean_data$r)$x, density(clean_data$r)$y)

# Extract subject-level parameters
alpha_subject_samples <- posterior_samples$alpha_subject  # Subject-level alpha parameters
sigma_subject_samples <- posterior_samples$sigma_subject  # Subject-level sigma parameters

# Assuming you have the number of subjects, S
S <- length(unique(clean_data$subject_id))  # Number of subjects (assuming subject_id is available)

# Create a dataframe to store the summary statistics
subject_ids <- 1:S  # You can replace this with the actual subject_id if available
alpha_subject_mean <- apply(alpha_subject_samples, 2, mean)  # Mean of alpha_subject for each subject
sigma_subject_mean <- apply(sigma_subject_samples, 2, mean)  # Mean of sigma_subject for each subject

# Create a dataframe with subject IDs, mean alpha, and mean sigma
df_subject_params <- data.frame(
  subject_id = subject_ids,
  alpha_subject = alpha_subject_mean,
  sigma_subject = sigma_subject_mean
)

# Optionally, add 95% credible intervals for each subject
alpha_subject_ci <- apply(alpha_subject_samples, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
sigma_subject_ci <- apply(sigma_subject_samples, 2, function(x) quantile(x, probs = c(0.025, 0.975)))

# Adding credible intervals to the dataframe
df_subject_params$alpha_subject_lower <- alpha_subject_ci[1, ]
df_subject_params$alpha_subject_upper <- alpha_subject_ci[2, ]
df_subject_params$sigma_subject_lower <- sigma_subject_ci[1, ]
df_subject_params$sigma_subject_upper <- sigma_subject_ci[2, ]

# View the dataframe
print(df_subject_params)

full_data <- read.csv("data_clean.csv") 

results <- df_subject_params %>% left_join(., full_data %>% select(subject_id = participant, pm_count) %>% distinct(subject_id, pm_count))

results$cc_count <- clean_data %>% count(subject_id) %>% pull(n)
results$b = results$sigma_subject ^ -results$alpha_subject

cor(results[, c(2, 3, 8, 9, 10)], use = "pairwise.complete.obs")
