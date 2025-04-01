data {
  int<lower=1> N;
  vector<lower=0>[N] clock_check_time;
  int<lower=1> K_k;
  int<lower=1> K_lambda;
  matrix[N, K_k] X_k;
  matrix[N, K_lambda] X_lambda;
}

parameters {
  vector[K_k] k;
  vector[K_lambda] lambda;
}

transformed parameters {
  vector<lower=0>[N] k_pred = exp(X_k * k);
  vector<lower=0>[N] lambda_pred = exp(X_lambda * lambda);
}

model {
  k ~ normal(0, 1);
  lambda ~ normal(0, 1);
  
  // Vectorized likelihood
  target += weibull_lpdf(clock_check_time | k_pred, lambda_pred);
}
