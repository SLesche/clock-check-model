functions {
  real clock_check_lik(real x, real t_target, real k, real g, real threshold, real a, int m) {
    real integral = 0.0;
    real delta_x = (x - 0.01) / m;  // Step size for the trapezoidal rule

    // Precompute normal CDF values for the range
    real normal_cdf_left[m];
    real normal_cdf_right[m];
    
    for (i in 1:m) {
      real t_left = 0.01 + (i - 1) * delta_x;
      real t_right = 0.01 + i * delta_x;

      // Store normal CDF values in arrays
      normal_cdf_left[i] = normal_cdf(t_target, t_left, k * t_left);
      normal_cdf_right[i] = normal_cdf(t_target, t_right, k * t_right);
    }

    // Loop over M intervals to compute the integral using the trapezoidal rule
    for (i in 1:m) {
      real hazard_left = g + (1 - g) / (1 + exp(-a * (1 - normal_cdf_left[i]) - threshold));
      real hazard_right = g + (1 - g) / (1 + exp(-a * (1 - normal_cdf_right[i]) - threshold));

      // Add the area of this trapezoid to the integral
      integral += (hazard_left + hazard_right) * delta_x / 2.0;
    }

    // Compute the survival probability and action probability
    real survival_prob = 1 - fmax(integral, 0.0);  // Avoid negative values
    real action_prob = g + (1 - g) / (1 + exp(-a * (1 - normal_cdf(t_target, x, k * x)) - threshold));

    // Return the log-likelihood
    return log(fmax(survival_prob * action_prob, 1e-10));  // Avoid log(0)
  }
}

data {
  int<lower=1> N;
  real<lower=0> known_t_to_target[N];
  real<lower=0> clock_check_time[N];
  int<lower=1> m; // Precision of integral estimation
}

parameters {
  real<lower=0> k;
  real<lower=0,upper=1> g;
  real<lower=0,upper=1> threshold;
  real<lower=0> a;
}

model {
  // Priors
  k ~ normal(1, 0.5);
  g ~ beta(2, 6);
  threshold ~ beta(2, 6);
  a ~ normal(10, 1);

  // Likelihood of observed data
  for (n in 1:N) {
    target += clock_check_lik(clock_check_time[n], known_t_to_target[n], k, g, threshold, a, m);
  }
}
