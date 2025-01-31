functions {
  real hazard_function(real x, real k, real g, real c){
    return(g + (1 - g) / (1 + exp((1 - x)/kx - c)))
  }
  
  // Define the hazard function
  real hazard_integrand(real x[], real xc, real[] theta, real[] x_r, int[] x_i) {
    // Extract parameters
    real k = theta[1];
    real g = theta[2];
    real c = theta[3];
    return g + (1 - g) / (1 + exp((1 - x) / (k * x) - c));
  }
  
  // Define the cumulative hazard function using integrate_1d
  real cumulative_hazard(real x, real k, real g, real c) {
    // Define the integration limits
    real a = 0; // Lower bound of integration
    real b = x; // Upper bound of integration
    
    // Define theta to pass parameters to the hazard_function
    real theta[3] = {k, g, c};
    
    // Call integrate_1d to compute the integral of hazard_function from a to b
    return integrate_1d(hazard_integrand, a, b, theta, rep_array(0.0, 0), rep_array(0, 0), 1e-6);
  }
  

  real clock_check_lik(real x, real k, real g, real c) {
    real inst_hazard = hazard_function(x, k, g, c);
    real survival = exp(-cumulative_hazard(x, k, g, c));
    
    return(inst_hazard*survival)
  }
}

data {
  int<lower=1> N;
  real<lower=0> known_t_to_target[N];
  real<lower=0> clock_check_time[N];
}

parameters {
  real<lower=0> k;
  real<lower=0,upper=1> g;
  real c;
}

model {
  // Priors
  k ~ gamma(2, 4);
  g ~ beta(2, 6);
  c ~ cauchy(0, 10);

  // Likelihood of observed data
  for (n in 1:N) {
    target += clock_check_lik(clock_check_time[n], k, g, c);
  }
}
