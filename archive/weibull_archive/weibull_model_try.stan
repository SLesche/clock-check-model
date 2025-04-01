data {
  int<lower=1> N;         // Number of failures
  real<lower=0> times[N]; // Failure times (sorted)
}

parameters {
  real<lower=0> lambda_0; // Constant hazard rate
  real<lower=0> eta;      // Scale parameter for power-law component
  real<lower=1> beta;     // Shape parameter (> 1 for increasing hazard rate)
}

model {
  // Priors
  lambda_0 ~ normal(0, 10);     // Prior for constant hazard rate
  eta ~ normal(0, 10);          // Prior for scale parameter
  beta ~ normal(2, 1);          // Prior for shape parameter

  // Log-likelihood
  for (i in 1:N) {
    target += log(lambda_0 + eta * beta * pow(times[i], beta - 1));
  }
  target += -lambda_0 * times[N] + eta * pow(times[N], beta) / beta;
}
