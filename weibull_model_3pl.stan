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
}

parameters {
  real<lower=0> a;                  // Scale parameter
  real<lower=0> b;                  // Shape parameter
  real c;                           // Growth/decay parameter (can be negative or positive)
}

model {
  // Priors
  a ~ gamma(2, 2);                  // Prior for scale parameter
  b ~ gamma(2, 2);                  // Prior for shape parameter
  c ~ normal(0, 1);                 // Prior for growth/decay parameter

  // Likelihood
  for (n in 1:N) {
    target += modified_weibull_lpdf(clock_check_time[n] | a, b, c);
  }
}
