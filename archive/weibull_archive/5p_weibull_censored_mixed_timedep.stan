data {
  int<lower=1> N;
  vector<lower=0>[N] times;
  int<lower=0,upper=1> event[N];

  int<lower=1> K_k;
  int<lower=1> K_lambda;
  int<lower=1> K_lambda2;
  int<lower=1> K_mixture;

  matrix[N, K_k] X_k;  
  matrix[N, K_lambda] X_lambda;    
  matrix[N, K_lambda2] X_lambda2;
  matrix[N, K_mixture] X_mixture;
}

parameters {
  vector[K_k] k;   // Coefficients for adult mortality shape
  vector[K_k] k2;  // Coefficients for infant mortality shape
  vector[K_lambda] lambda;  
  vector[K_lambda2] lambda2;
  vector[K_mixture] mixture_rate;
}

transformed parameters {
  vector<lower=1>[N] k_pred = 1 + exp(X_k * k);  // Adult mortality (k > 1)
  vector<lower=0, upper=1>[N] k2_pred = inv_logit(X_k * k2);  // Infant mortality (k < 1)

  vector<lower=0>[N] lambda_pred = exp(X_lambda * lambda);
  vector<lower=0>[N] lambda2_pred = exp(X_lambda2 * lambda2);
  vector<lower=0, upper=1>[N] mixture_prob = inv_logit(X_mixture * mixture_rate);
}

model {
  k ~ normal(0, 1);  // Centered around 1 after transformation
  k2 ~ normal(0, 1);  // Constrained via inv_logit
  lambda ~ normal(0, 1);
  lambda2 ~ normal(0, 1);
  mixture_rate ~ normal(0, 1);

  for (n in 1:N) {
    if (event[n] == 1) {
      target += log_mix(mixture_prob[n], 
        weibull_lpdf(times[n] | k2_pred[n], lambda2_pred[n]),
        weibull_lpdf(times[n] | k_pred[n], lambda_pred[n])
      );
    } else {
      target += log_mix(mixture_prob[n], 
        weibull_lccdf(times[n] | k2_pred[n], lambda2_pred[n]),
        weibull_lccdf(times[n] | k_pred[n], lambda_pred[n])
      );
    }
  }
}
