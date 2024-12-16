n_participants = 100
mu_sigma_0 = 0.5
sigma_sigma_0 = 0.3

mu_k = 1
sigma_k = 0.4

participant_sigma = rnorm(1, mu_sigma_0, sigma_sigma_0)
participant_k = rnorm(1, mu_k, sigma_k)

total_time = 300
t_to_target = 100 / total_time
time_since_last_cc = 50 / total_time
sigma_est = participant_sigma[1] * exp(participant_k[1]*time_since_last_cc)

prob_check = 1 - pnorm(t_to_target, time_since_last_cc, sigma_est)
