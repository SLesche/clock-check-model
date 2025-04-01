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
}

parameters {
  vector[K_k] k;            
  vector[K_lambda] lambda;  
  vector[K_lambda2] lambda2;
  vector[K_mixture] mixture_rate; 
}

transformed parameters {
  vector<lower=0>[N] k_pred = exp(X_k * k);
  vector<lower=0>[N] lambda_pred = exp(X_lambda * lambda);
  vector<lower=0>[N] lambda2_pred = exp(X_lambda2 * lambda2);
  vector<lower=0, upper=1>[N] mixture_prob = inv_logit(X_mixture * mixture_rate);
}

model {
  k ~ normal(0, 1);
  lambda ~ normal(0, 1);
  lambda2 ~ normal(0, 1);
  mixture_rate ~ normal(0, 1);

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
