data {
  int<lower=1> N;  // Total number of observations (censored + uncensored)
  vector<lower=0>[N] times;  // Event or censoring times
  int<lower=0,upper=1> event[N];  // 1 = event observed, 0 = censored

  int<lower=1> K_k; 
  int<lower=1> K_lambda;
  int<lower=1> K_lambda2;  
  int<lower=1> K_mixture;

  matrix[N, K_k] X_k;  
  matrix[N, K_lambda] X_lambda;    
  matrix[N, K_lambda2] X_lambda2;    
  matrix[N, K_mixture] X_mixture;

  int<lower=1> S;  // Number of subjects
  int<lower=1, upper=S> subject_id[N];  // Subject identifiers
}

parameters {
  // Fixed effects
  vector[K_k] k_global;            
  vector[K_lambda] lambda_global;  
  vector[K_lambda2] lambda2_global;  
  vector[K_mixture] mixture_rate_global; 

  // Random effects (subject-specific)
  vector[K_k] k[S];            
  vector[K_lambda] lambda[S];  
  vector[K_lambda2] lambda2[S];  
  vector[K_mixture] mixture_rate[S]; 
  
  // Hyperpriors for the random effects
  real<lower=0> sigma_k;  
  real<lower=0> sigma_lambda;
  real<lower=0> sigma_lambda2;
  real<lower=0> sigma_mixture_rate;
}

transformed parameters {
  vector<lower=0>[N] k_pred;            
  vector<lower=0>[N] lambda_pred;  
  vector<lower=0>[N] lambda2_pred;  
  vector<lower=0, upper=1>[N] mixture_prob; 

  for (n in 1:N) {
    k_pred[n] = exp(X_k[n] * (k_global + k[subject_id[n]]));  // Subject-specific k
    lambda_pred[n] = exp(X_lambda[n] * (lambda_global + lambda[subject_id[n]]));  // Subject-specific lambda
    lambda2_pred[n] = exp(X_lambda2[n] * (lambda2_global + lambda2[subject_id[n]]));  // Subject-specific lambda2
    mixture_prob[n] = inv_logit(X_mixture[n] * (mixture_rate_global + mixture_rate[subject_id[n]]));  // Subject-specific mixture rate
  }
}

model {
  // Priors for fixed effects
  k_global ~ normal(0, 1);
  lambda_global ~ normal(0, 1);
  lambda2_global ~ normal(0, 1);
  mixture_rate_global ~ normal(0, 1);

  // Priors for random effects (subject-specific parameters)
  for (s in 1:S) {
    k[s] ~ normal(0, sigma_k);            
    lambda[s] ~ normal(0, sigma_lambda);  
    lambda2[s] ~ normal(0, sigma_lambda2);  
    mixture_rate[s] ~ normal(0, sigma_mixture_rate); 
  }

  // Hyperpriors for the random effect standard deviations
  sigma_k ~ gamma(2, 4);
  sigma_lambda ~ gamma(2, 4);
  sigma_lambda2 ~ gamma(2, 4);
  sigma_mixture_rate ~ gamma(2, 4);

  // Likelihood
  for (n in 1:N) {
    if (event[n] == 1) {
      target += log_mix(mixture_prob[n], 
        weibull_lpdf(times[n] | 1, lambda2_pred[n]),
        weibull_lpdf(times[n] | k_pred[n], lambda_pred[n])
      );
    } else {
      target += log_mix(mixture_prob[n], 
        weibull_lccdf(times[n] | 1, lambda2_pred[n]),
        weibull_lccdf(times[n] | k_pred[n], lambda_pred[n])
      );
    }
  }
}
