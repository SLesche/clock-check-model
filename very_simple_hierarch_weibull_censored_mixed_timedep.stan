data {
  int<lower=1> N;  // Number of observations (censored + uncensored)
  int<lower=1> J;  // Number of groups/individuals
  int<lower=1> K_k; 
  int<lower=1> K_lambda;
  int<lower=1> K_mixture;

  vector<lower=0>[N] times;  // Event or censoring times
  int<lower=0,upper=1> event[N];  // 1 = event observed, 0 = censored

  matrix[N, K_k] X_k;  
  matrix[N, K_lambda] X_lambda;    
  matrix[N, K_mixture] X_mixture;

  int<lower=1,upper=J> group[N];  // Group or individual index for hierarchical structure
}

parameters {
  // Population-level parameters
  vector[K_k] mu_k;
  vector[K_lambda] mu_lambda;
  vector[K_mixture] mu_mixture;

  // Group-level (hierarchical) parameters
  vector[K_k] k_raw[J];
  vector[K_lambda] lambda_raw[J];
  vector[K_mixture] mixture_rate_raw[J];

  // Hyperparameters (standard deviations)
  vector<lower=0>[K_k] sigma_k;
  vector<lower=0>[K_lambda] sigma_lambda;
  vector<lower=0>[K_mixture] sigma_mixture;
}

transformed parameters {
  vector<lower=0>[N] k_pred;
  vector<lower=0>[N] lambda_pred;
  vector<lower=0, upper=1>[N] mixture_prob;

  for (n in 1:N) {
    int g = group[n];
    k_pred[n] = exp(X_k[n] * (mu_k + sigma_k .* k_raw[g]));
    lambda_pred[n] = exp(X_lambda[n] * (mu_lambda + sigma_lambda .* lambda_raw[g]));
    mixture_prob[n] = inv_logit(X_mixture[n] * (mu_mixture + sigma_mixture .* mixture_rate_raw[g]));
  }
}

model {
  // Priors on hyperparameters
  mu_k ~ normal(0, 1);
  mu_lambda ~ normal(0, 1);
  mu_mixture ~ normal(0, 1);

  sigma_k ~ normal(0, 1);
  sigma_lambda ~ normal(0, 1);
  sigma_mixture ~ normal(0, 1);

  // Hierarchical priors
  for (j in 1:J) {
    k_raw[j] ~ normal(0, 1);
    lambda_raw[j] ~ normal(0, 1);
    mixture_rate_raw[j] ~ normal(0, 1);
  }

  // Likelihood
  for (n in 1:N) {
    if (event[n] == 1) {
      target += log_mix(mixture_prob[n], 
        weibull_lpdf(times[n] | 1, 1),
        weibull_lpdf(times[n] | k_pred[n], lambda_pred[n])
      );
    } else {
      target += log_mix(mixture_prob[n], 
        weibull_lccdf(times[n] | 1, 1),
        weibull_lccdf(times[n] | k_pred[n], lambda_pred[n])
      );
    }
  }
}
