functions {
  real mixture_weibull_lpdf(real[] t, real pi, real alpha1, real sigma1, real alpha2, real sigma2) {
    real log_prob = 0;
    for (i in 1:num_elements(t)) {
      // Density of Weibull component 1 (constant hazard: alpha1 ~ 1)
      real f1 = exp(weibull_lpdf(t[i] | alpha1, sigma1));
      
      // Density of Weibull component 2 (flexible hazard: alpha2, sigma2 free)
      real f2 = exp(weibull_lpdf(t[i] | alpha2, sigma2));
      
      // Mixture probability
      log_prob += log(pi * f1 + (1 - pi) * f2);
    }
    return log_prob;
  }
}

data {
  int<lower=1> N;                  // Number of observations
  real<lower=0> clock_check_time[N]; // Observation times
}

parameters {
  real<lower=0, upper=1> pi;       // Mixture weight
  real<lower=0> alpha1;            // Shape parameter for constant hazard Weibull
  real<lower=0> sigma1;            // Scale parameter for constant hazard Weibull
  real<lower=0> alpha2;            // Shape parameter for flexible hazard Weibull
  real<lower=0> sigma2;            // Scale parameter for flexible hazard Weibull
}

model {
  // Priors
  pi ~ beta(2, 2);                 // Mixture weight
  alpha1 ~ normal(1, 0.01;         // Shape for constant hazard Weibull (close to 1)
  sigma1 ~ normal(1, 1);           // Scale for constant hazard Weibull
  alpha2 ~ gamma(2, 2);            // Shape for flexible hazard Weibull
  sigma2 ~ normal(1, 1);           // Scale for flexible hazard Weibull

  // Likelihood
  target += mixture_weibull_lpdf(clock_check_time | pi, alpha1, sigma1, alpha2, sigma2);
}
