functions {
  // Hazard function integrand for integrate_1d
  real hazard_integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    real k = theta[1];
    real g = theta[2];
    real c = theta[3];
    
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
  real<lower=0> eta; // Precision of the timing intervall
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

transformed parameters {
  real<lower=0, upper=1> integration_stop;

  // Calculate the integration stop from 0 to the intervall where timing becomes relevant
  integration_stop = 1 / (k*(log((1-g)/eta - 1) + c) +1);
}

model {
  // Priors
  k ~ gamma(2, 4);
  g ~ beta(2, 6);
  c ~ cauchy(0, 10);

  // Loop through observations
  for (n in 1:N) {
    // Regularize clock_check_time to avoid issues at 0
    real t_safe = fmax(clock_check_time[n], 1e-10);
    real hazard = g + (1 - g) / (1 + exp((1 - t_safe) / (k * t_safe) - c));

    real cum_hazard_large = 0.0;
    // Split integration logic
    if (t_safe > integration_stop) {
      cum_hazard_large = integrate_1d(
        hazard_integrand,
        integration_stop,  // Lower bound offset to avoid zero
        t_safe,            // Upper bound
        {k, g, c},
        x_r,
        x_i,
        1e-10
      );
    }
    
    real cum_hazard_small = g * integration_stop;
    
    // Log-likelihood contribution
    target += log(clock_check_lik(hazard, cum_hazard_small + cum_hazard_large));
  }
}
