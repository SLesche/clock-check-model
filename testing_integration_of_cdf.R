# Load required library
library(stats)

# Standard normal PDF
standard_normal_pdf <- function(x) {
  (1 / sqrt(2 * pi)) * exp(-x^2 / 2)
}

# Standard normal CDF
standard_normal_cdf <- function(x) {
  pnorm(x)  # Built-in function in R for the standard normal CDF
}

# Analytical solution for the integral of the CDF
integral_of_cdf_analytical <- function(x) {
  x * pnorm(x) + dnorm(x)  # x * CDF(x) + PDF(x)
}

# Numerical integration of the CDF from -Inf to x
integral_of_cdf_numerical <- function(x) {
  integrate(function(t) standard_normal_cdf(t), lower = -Inf, upper = x)$value
}

# Test the functions
x_values <- seq(-3, 3, by = 0.01)  # Values of x for testing

# Calculate results
results <- data.frame(
  x = x_values,
  Analytical = sapply(x_values, integral_of_cdf_analytical),
  Numerical = sapply(x_values, integral_of_cdf_numerical)
)

# Print the results
print(results)

# Check the difference between analytical and numerical results
results$Difference <- abs(results$Analytical - results$Numerical)
print(results)

# Plot the results for visual comparison
plot(results$x, results$Analytical, type = "l", col = "blue", lwd = 2,
     ylab = "Integral of CDF", xlab = "x", main = "Integral of Standard Normal CDF")
lines(results$x, results$Numerical, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Analytical", "Numerical"), col = c("blue", "red"),
       lty = c(1, 2), lwd = 2)
