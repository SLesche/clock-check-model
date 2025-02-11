data {
  int<lower=1> N;  // Total number of observations (censored + uncensored)
  vector<lower=0>[N] times;  // Event or censoring times
  int<lower=0,upper=1> event[N];  // 1 = event observed, 0 = censored

  int<lower=1> K_k; 
  int<lower=1> K_lambda;
  int<lower=1> K_mixture;

  matrix[N, K_k] X_k;  
  matrix[N, K_lambda] X_lambda;    
  matrix[N, K_mixture] X_mixture;
  
  int<lower=1> S;  // Number of subjects
  int<lower=1, upper=S> subject_id[N];  // Subject identifiers
}

parameters {
  // Global fixed effects (population-level means)
  vector[K_k] k_global;            
  vector[K_lambda] lambda_global;  
  vector[K_mixture] mixture_rate_global; 

  // Random intercepts for each subject
  vector[S] k_sub;            
  vector[S] lambda_sub;  
  vector[S] mixture_rate_sub; 
  
  // Hyperpriors for the random effects (variance across subjects)
  real<lower=0> sigma_k;  
  real<lower=0> sigma_lambda;
  real<lower=0> sigma_mixture_rate;
}

transformed parameters {
  vector<lower=0>[N] k_pred;
  vector<lower=0>[N] lambda_pred;
  vector<lower=0, upper=1>[N] mixture_prob;

  for (n in 1:N) {
    int s = subject_id[n];  // Get the subject index for observation n
    
    // Regression + subject-specific random intercept
    k_pred[n] = exp(X_k[n] * k_global + k_sub[s]);  
    lambda_pred[n] = exp(X_lambda[n] * lambda_global + lambda_sub[s]);  
    mixture_prob[n] = inv_logit(X_mixture[n] * mixture_rate_global + mixture_rate_sub[s]);  
  }
}

model {
  // Priors for fixed effects (global means)
  k_global ~ normal(0, 1);
  lambda_global ~ normal(0, 1);
  mixture_rate_global ~ normal(0, 1);

  // Priors for hyperparameters (variance of random effects)
  sigma_k ~ gamma(2, 8);
  sigma_lambda ~ gamma(2, 8);
  sigma_mixture_rate ~ gamma(2, 8);

  // Hierarchical structure: subject-level intercepts drawn from global regression mean
  k_sub ~ normal(0, sigma_k);
  lambda_sub ~ normal(0, sigma_lambda);
  mixture_rate_sub ~ normal(0, sigma_mixture_rate);

  // Likelihood function
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
