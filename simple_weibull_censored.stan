data {
  int<lower=1> N_uncensored;  // Number of uncensored events
  int<lower=1> N_censored;    // Number of censored events

  vector<lower=0>[N_uncensored] uncensored_times;  // Times for event = 1
  vector<lower=0>[N_censored] censored_times;      // Times for event = 0

  int<lower=1> K_k; 
  int<lower=1> K_lambda;
  matrix[N_uncensored, K_k] X_k_uncensored;  // Predictor matrix for uncensored times
  matrix[N_censored, K_k] X_k_censored;      // Predictor matrix for censored times

  matrix[N_uncensored, K_lambda] X_lambda_uncensored;  
  matrix[N_censored, K_lambda] X_lambda_censored;    
}

parameters {
  vector[K_k] k;            // Regression coefficients for k
  vector[K_lambda] lambda;  // Regression coefficients for lambda
}

transformed parameters {
  vector<lower=0>[N_uncensored] k_pred_uncensored;  
  vector<lower=0>[N_uncensored] lambda_pred_uncensored;

  vector<lower=0>[N_censored] k_pred_censored;  
  vector<lower=0>[N_censored] lambda_pred_censored;

  k_pred_uncensored = exp(X_k_uncensored * k);
  lambda_pred_uncensored = exp(X_lambda_uncensored * lambda);

  k_pred_censored = exp(X_k_censored * k);
  lambda_pred_censored = exp(X_lambda_censored * lambda);
}


model {
  // Priors
  k ~ normal(0, 1);
  lambda ~ normal(0, 1);

  // Likelihood for uncensored
  target += weibull_lpdf(uncensored_times | k_pred_uncensored, lambda_pred_uncensored);

  // Likelihood for censored
  target += weibull_lccdf(censored_times | k_pred_censored, lambda_pred_censored);
}

// generated quantities {
//   // Posterior predictive samples for uncensored data
//   vector[N_uncensored] uncensored_pred;
//   for (n in 1:N_uncensored) {
//     uncensored_pred[n] = weibull_rng(k_pred_uncensored[n], lambda_pred_uncensored[n]);
//   }
// 
//   // vector[N_uncensored] log_lik_uncensored;
//   // log_lik_uncensored = weibull_lpdf(uncensored_times | k_pred_uncensored, lambda_pred_uncensored);
// 
//   // Posterior predictive samples for censored data
//   vector[N_censored] censored_pred;
//   for (n in 1:N_censored) {
//     censored_pred[n] = weibull_rng(k_pred_censored[n], lambda_pred_censored[n]);
//   }
// 
//   // vector[N_censored] log_lik_censored;
//   // log_lik_censored = weibull_lccdf(censored_times | k_pred_censored, lambda_pred_censored);
// 
// }
