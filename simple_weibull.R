library(dplyr)
library(survival)

data <- read.csv("archive/diffusion_data.csv")

clean_data <- data %>% 
  # filter(
  #   block_duration == known_t_to_target
  # ) %>%
  mutate(
    r = time_since_last_cc / known_t_to_target
  ) %>% 
  filter(r > 0.01, known_t_to_target > 0) # make sure that the data is not too close to 0 (this cause issues in integration)

weibull_model <- survreg(Surv(clean_data$r) ~ known_t_to_target, data = clean_data, dist = "weibull")

summary(weibull_model)
library(ggplot2)

# Assuming 'weibull_model' is already fitted and 'data' contains your actual data
# Get the predicted values from the Weibull regression model
predicted_values <- predict(weibull_model, type = "response")

hist(predicted_values)
# Add predicted values to the data
data$predicted_r <- predicted_values

# Plot the data and the fitted model
ggplot(data, aes(x = known_t, y = r)) +
  geom_point(color = "blue", size = 3) +  # Original data points
  geom_line(aes(x = known_t, y = predicted_r), color = "red", size = 1) +  # Fitted Weibull curve
  labs(title = "Weibull Regression Fit", x = "Known_t", y = "r") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the plot title