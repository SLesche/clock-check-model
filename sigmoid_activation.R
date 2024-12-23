# Load necessary library
library(pracma)  # For numerical integration

# Define the functions
prob_target <- function(t_target, t_step, k){
  return(1 - pnorm(t_target, t_step, k*t_step))
}

prob_action <- function(prob_target, g, threshold, a = 10){
  return(g + (1 - g)/(1 + exp(-a*(prob_target - threshold))))
}

compute_pdf <- function(t_steps, t_target, k, g, threshold, a, dstep){
  n <- length(t_steps)
  f <- numeric(n)
  cumulative_prob <- numeric(n)
  
  for (i in 1:n) {
    # Compute prob_target for this t_step
    p_target <- prob_target(t_target, t_steps[i], k)
    
    # Compute prob_action
    p_action <- prob_action(p_target, g, threshold, a)
    
    # Compute the integral of f up to t_step[i]
    if (i > 1) {
      cumulative_prob[i] <- trapz(t_steps[1:(i - 1)], f[1:(i- 1)])
    }
    
    # Compute f(t_step)
    f[i] <- (1 - cumulative_prob[i]) * p_action
  }
  
  return(f)
}

# Example usage
t_target <- 0.2                    # Target time
t_steps <- seq(0, t_target, by = 0.01)  # Time steps from 0.1 to 10
k <- 1                      # Scaling factor for prob_target
g <- 0.1                  # Base probability for prob_action
threshold <- 0.5       # Threshold for prob_action
a <- 10                            # Sharpness of the sigmoid in prob_action

# Compute the scaled PDF
pdf_values <- compute_pdf(t_steps, t_target, k, g, threshold, a)

# Plot the PDF
plot(t_steps, pdf_values, type = "l", col = "blue", lwd = 2,
     xlab = "t_step", ylab = "Probability Density",
     main = "Scaled PDF of Action Timing")

