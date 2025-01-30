functions {
  real nhpp_weibull_mixture_lpdf(int[] event_count, real[] event_time, real alpha, real sigma, real lambda_g) {
    real log_prob = 0;
    for (i in 1:num_elements(event_time)) {
      real lambda_t = (alpha / sigma) * pow(event_time[i] / sigma, alpha - 1); // Weibull intensity function
      real f_nhpp = lambda_t * exp(-lambda_t * event_time[i]);
      real f_guess = lambda_g * exp(-lambda_g * event_time[i]); // Constant intensity for random guessing
      log_prob += log(f_nhpp + f_guess);
    }
    return log_prob;
  }
}

data {
  int<lower=1> N;                  // Number of observations
  real<lower=0> event_time[N];      // Clock check times
  int<lower=0> event_count[N];      // Number of events per interval
}

parameters {
  real<lower=0> alpha;              // Shape parameter for Weibull intensity
  real<lower=0> sigma;              // Scale parameter for Weibull intensity
  real<lower=0> lambda_g;           // Constant rate for random guessing
}

model {
  // Priors
  alpha ~ gamma(2, 2);
  sigma ~ normal(1, 1);
  lambda_g ~ normal(1, 1);

  // Likelihood
  target += nhpp_weibull_mixture_lpdf(event_count, event_time, alpha, sigma, lambda_g);
}