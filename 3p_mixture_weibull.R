# Load necessary library
library(ggplot2)

# 3PL Hazard Function
hazard_3pl <- function(t, a, b, c) {
  c + (1 - c) / (1 + exp(-a * (t - b)))
}

# Cumulative Hazard Function
cum_hazard_3pl <- function(t, a, b, c) {
  ct <- c * t
  logistic_int <- (1 / a) * log(1 + exp(a * (t - b))) - (1 / a) * log(1 + exp(-a * b))
  return(ct + (1 - c) * logistic_int)
}

# Survival Function (S(t))
survival_3pl <- function(t, a, b, c) {
  H_t <- cum_hazard_3pl(t, a, b, c)
  return(exp(-H_t))
}

# PDF (Failure Probability Function, F(t) = 1 - S(t))
pdf_3pl <- function(t, a, b, c) {
  return(survival_3pl(t, a, b, c) * hazard_3pl(t, a, b, c))
}

# Define parameter values
a <- 2000    # Steepness of the hazard increase
b <- 0    # Threshold time where failure rate increases rapidly
c <- 0.4    # Baseline failure rate (guessing component)

# Generate time points
t_vals <- seq(0, 1, length.out = 100)

# Compute function values
haz_vals <- sapply(t_vals, hazard_3pl, a, b, c)
cum_haz_vals <- sapply(t_vals, cum_hazard_3pl, a, b, c)
surv_vals <- sapply(t_vals, survival_3pl, a, b, c)
pdf_vals <- sapply(t_vals, pdf_3pl, a, b, c)

# Create a data frame for plotting
df <- data.frame(
  t = t_vals,
  Hazard = haz_vals,
  Cumulative_Hazard = cum_haz_vals,
  Survival = surv_vals,
  PDF = pdf_vals
)

# Plot the functions
ggplot(df, aes(x = t)) +
  geom_line(aes(y = Hazard, color = "Hazard Rate")) +
  # geom_line(aes(y = Cumulative_Hazard, color = "Cumulative Hazard")) +
  geom_line(aes(y = Survival, color = "Survival Function")) +
  geom_line(aes(y = PDF, color = "PDF (Failure Probability)")) +
  labs(
    title = "3PL-Based Hazard Model",
    x = "Time",
    y = "Function Value",
    color = "Functions"
  ) +
  theme_minimal()

