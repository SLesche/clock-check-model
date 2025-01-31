functions {
  real mixture_weibull_lpdf(real[] t, real pi, real sigma1, real alpha2, real[] sigma2) {
    real log_prob = 0;
    for (i in 1:num_elements(t)) {
      // Density of Weibull component 1 (constant hazard: alpha1 ~ 1)
      real f1 = exp(weibull_lpdf(t[i] | 1, sigma1));
      
      // Density of Weibull component 2 (flexible hazard: alpha2, sigma2[i])
      real f2 = exp(weibull_lpdf(t[i] | alpha2, sigma2[i]));
      
      // Mixture probability
      log_prob += log(pi * f1 + (1 - pi) * f2);
    }
    return log_prob;
  }
}

data {
  int<lower=1> N;                  // Number of observations
  real<lower=0> clock_check_time[N]; // Observation times
  real<lower=0> known_t_to_target[N]; // Time till target known by individuals
}

parameters {
  real<lower=0, upper=1> pi;       // Mixture weight
  real<lower=0> sigma1;            // Scale parameter for constant hazard Weibull
  real<lower=0> alpha2;            // Shape parameter for flexible hazard Weibull
  real beta0;                      // Intercept for sigma2 model
  real beta1;                      // Slope for sigma2 model
}

transformed parameters {
  real sigma2[N];                  // Scale parameter for flexible hazard Weibull (varying)
  for (i in 1:N) {
    sigma2[i] = exp(beta0 + beta1 * known_t_to_target[i]);
  }
}

model {
  // Priors
  pi ~ beta(2, 2);
  sigma1 ~ normal(1, 1);
  alpha2 ~ gamma(2, 2);
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, 1);

  // Likelihood
  target += mixture_weibull_lpdf(clock_check_time | pi, sigma1, alpha2, sigma2);
}
