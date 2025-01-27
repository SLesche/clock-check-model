functions {
  // Hazard function integrand for integrate_1d
  real hazard_integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    real k = theta[1];
    real g = theta[2];
    real c = theta[3];
    
    if (x < 1e-2){
      return g;
    }
    
    return g + (1 - g) / (1 + exp((1 - x) / (k * x) - c));
  }

  // Likelihood function combining hazard and survival
  real clock_check_lik(real hazard, real cum_hazard) {
    return hazard * exp(-cum_hazard);
  }
}

data {
  int<lower=1> N;                     // Number of observations
  real<lower=0> clock_check_time[N];  // Observation times
}

transformed data {
  array[0] real x_r; // No real data constants
  array[0] int x_i;  // No integer data constants
}

parameters {
  real<lower=0> k;
  real<lower=0, upper=1> g;
  real c;
}

model {
  // Priors
  k ~ gamma(2, 4);
  g ~ beta(2, 6);
  c ~ cauchy(0, 10);

  // Likelihood
  for (n in 1:N) {
    // Regularize clock_check_time to avoid issues at 0
    real t_safe = fmax(clock_check_time[n], 1e-10);
    real hazard = g + (1 - g) / (1 + exp((1 - t_safe) / (k * t_safe) - c));
    real cum_hazard = integrate_1d(
      hazard_integrand,
      0.0, // Lower bound offset to avoid zero
      t_safe,       // Upper bound
      {k, g, c},
      x_r,
      x_i
    );
    target += log(clock_check_lik(hazard, cum_hazard));
  }
}
