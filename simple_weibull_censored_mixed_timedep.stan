data {
  int<lower=1> N_uncensored;  // Number of uncensored events
  int<lower=1> N_censored;    // Number of censored events

  vector<lower=0>[N_uncensored] uncensored_times;  // Times for event = 1
  vector<lower=0>[N_censored] censored_times;      // Times for event = 0

  int<lower=1> K_k; 
  int<lower=1> K_lambda;
  // int<lower=1> K_k2;          // Additional predictor dimensions for second Weibull
  int<lower=1> K_lambda2;     // Additional predictor dimensions for second Weibull
  int<lower=> K_mixture;
  matrix[N_uncensored, K_k] X_k_uncensored;  // Predictor matrix for uncensored times
  matrix[N_censored, K_k] X_k_censored;      // Predictor matrix for censored times

  matrix[N_uncensored, K_lambda] X_lambda_uncensored;  
  matrix[N_censored, K_lambda] X_lambda_censored;    

  // matrix[N_uncensored, K_k2] X_k2_uncensored;  // Predictor matrix for second Weibull uncensored
  // matrix[N_censored, K_k2] X_k2_censored;      // Predictor matrix for second Weibull censored

  matrix[N_uncensored, K_lambda2] X_lambda2_uncensored;  
  matrix[N_censored, K_lambda2] X_lambda2_censored;    
  
  matrix[N_uncensored, K_mixture] X_mixture_uncensored;
  matrix[N_censored, K_mixture] X_mixture_censored;
}

parameters {
  vector[K_k] k;            // Regression coefficients for k
  vector[K_lambda] lambda;  // Regression coefficients for lambda
  
  // vector[K_k2] k2;          // Regression coefficients for k2 (second Weibull)
  vector[K_lambda2] lambda2;  // Regression coefficients for lambda2 (second Weibull)
  
  vector[K_mixture] mixture_rate; // Regression for mixture probability
}

transformed parameters {
  vector<lower=0>[N_uncensored] k_pred_uncensored;  
  vector<lower=0>[N_uncensored] lambda_pred_uncensored;

  vector<lower=0>[N_censored] k_pred_censored;  
  vector<lower=0>[N_censored] lambda_pred_censored;

  // vector<lower=0>[N_uncensored] k2_pred_uncensored;  
  vector<lower=0>[N_uncensored] lambda2_pred_uncensored;

  // vector<lower=0>[N_censored] k2_pred_censored;  
  vector<lower=0>[N_censored] lambda2_pred_censored;
  
  vector<lower=0>[N_uncensored] mixture_rate_uncensored;
  vector<lower=0>[N_censored] mixture_rate_censored;

  // Predicted values for both Weibull distributions
  k_pred_uncensored = exp(X_k_uncensored * k);
  lambda_pred_uncensored = exp(X_lambda_uncensored * lambda);

  k_pred_censored = exp(X_k_censored * k);
  lambda_pred_censored = exp(X_lambda_censored * lambda);

  // k2_pred_uncensored = exp(X_k2_uncensored * k2);
  lambda2_pred_uncensored = exp(X_lambda2_uncensored * lambda2);

  // k2_pred_censored = exp(X_k2_censored * k2);
  lambda2_pred_censored = exp(X_lambda2_censored * lambda2);
  
  mixture_rate_uncensored = inv_logit(X_mixture_uncensored * mixture_rate);
  mixture_rate_censored = inv_logit(X_mixture_censored * mixture_rate);
}

model {
  // Priors
  k ~ normal(0, 1);
  lambda ~ normal(0, 1);
  
  // k2 ~ normal(0, 1);  // Prior for second Weibull coefficients
  lambda2 ~ normal(0, 1);  // Prior for second Weibull coefficients

  // Mixture prior
  mixture_rate ~ normal(0, 1);  // Beta prior for the mixture rate (between 0 and 1)

  // Likelihood for uncensored data
  target += log_mix(mixture_rate, 
    weibull_lpdf(uncensored_times | k_pred_uncensored, lambda_pred_uncensored),
    weibull_lpdf(uncensored_times | 1, lambda2_pred_uncensored)
  );

  // Likelihood for censored data
  target += log_mix(mixture_rate, 
    weibull_lccdf(censored_times | k_pred_censored, lambda_pred_censored),
    weibull_lccdf(censored_times | 1, lambda2_pred_censored)
  );
}

// generated quantities {
//   // Posterior predictive samples for uncensored data
//   vector[N_uncensored] uncensored_pred;
//   for (n in 1:N_uncensored) {
//     uncensored_pred[n] = weibull_rng(k_pred_uncensored[n], lambda_pred_uncensored[n]);
//   }
// 
//   // Posterior predictive samples for censored data
//   vector[N_censored] censored_pred;
//   for (n in 1:N_censored) {
//     censored_pred[n] = weibull_rng(k_pred_censored[n], lambda_pred_censored[n]);
//   }
// }
