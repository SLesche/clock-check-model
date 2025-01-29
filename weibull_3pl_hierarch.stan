functions {
  // Define the PDF for the 3PL modified Weibull distribution
  real modified_weibull_lpdf(real t, real a, real b, real c) {
    real term1 = log(a) + log(b + c * t) + (b - 1) * log(t) + c * t;
    real term2 = -a * pow(t, b) * exp(c * t);
    return term1 + term2;
  }
}

data {
  int<lower=1> N;                   // Number of observations
  real<lower=0> clock_check_time[N]; // Observation times
  int<lower=1> G;                   // Number of groups
  int<lower=1, upper=G> group[N];   // Group assignment for each observation
}

parameters {
  // Hyperparameters for hierarchical priors
  real<lower=0> mu_a;               // Mean of scale parameter (a)
  real<lower=0> mu_b;               // Mean of shape parameter (b)
  real mu_c;                        // Mean of growth/decay parameter (c)
  
  real<lower=0> sigma_a;            // Standard deviation of scale parameter
  real<lower=0> sigma_b;            // Standard deviation of shape parameter
  real<lower=0> sigma_c;            // Standard deviation of growth/decay parameter

  // Group-level parameters
  real<lower=0> a[G];               // Scale parameter for each group
  real<lower=0> b[G];               // Shape parameter for each group
  real<lower=0> c[G];                        // Growth/decay parameter for each group
}

model {
  // Priors for hyperparameters
  mu_a ~ gamma(2, 2);
  mu_b ~ gamma(2, 2);
  mu_c ~ gamma(2, 2);

  sigma_a ~ gamma(2, 2);
  sigma_b ~ gamma(2, 2);
  sigma_c ~ gamma(2, 2);

  // Priors for group-level parameters
  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  c ~ normal(mu_c, sigma_c);

  // Likelihood
  for (n in 1:N) {
    target += modified_weibull_lpdf(clock_check_time[n] | a[group[n]], b[group[n]], c[group[n]]);
  }
}
