data {
  int<lower=1> N;                  // Number of observations
  vector<lower=0>[N] clock_check_time; // Observation times
  int<lower=1> K_k; // Number of predictors
  int<lower=1> K_lambda; // for
  matrix[N, K_k] X_k; // Predictor Matrix for k
  matrix[N, K_lambda] X_lambda; // Predictor Matrix for lambda
}

parameters {
  vector[K_k] k;            // Regression coefficients for k
  vector[K_lambda] lambda;            // Regression coefficients for lambda
}

transformed parameters {
  vector<lower=0>[N] k_pred;  // Computed alpha for each observation
  vector<lower=0>[N] lambda_pred;  // Computed sigma for each observation

  for (n in 1:N) {
    k_pred[n] = exp(X_k[n] * k);
    lambda_pred[n] = exp(X_lambda[n] * lambda);
  }
}

model {
  // Priors
  k ~ normal(0, 1);         
  lambda ~ normal(0, 1);
  
  // Likelihood
  for (n in 1:N) {
    target += weibull_lpdf(clock_check_time[n] | k_pred[n], lambda_pred[n]);
  }
}
