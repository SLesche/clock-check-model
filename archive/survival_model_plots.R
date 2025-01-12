library(tidyverse)

data <- read.csv("data_clean.csv")

# Plot of first clock_check
first_checks <- data %>% 
  group_by(participant, block_num) %>% 
  filter(clock_check == 1, time_since_start != 0, time_to_end != 0) %>% 
  slice_head(n = 1)

first_checks %>% pull(time_since_start) %>% hist()

times <- seq(0, max(first_checks$block_duration), 0.5)

surv_per <- purrr::map_dbl(times, \(x) sum(first_checks$time_since_start > x) / length(first_checks$time_since_start))

plot(times, surv_per)

checks <- data %>% 
  group_by(participant, block_num) %>% 
  filter(clock_check == 1, time_since_start != 0, time_to_end != 0, time_since_last_cc != 0)

checks <- checks %>% 
  mutate(
    rate_not_checked = time_since_last_cc / known_t_to_target
  ) 

checks %>%filter(known_t_to_target < 50) %>%  pull(rate_not_checked) %>% hist(breaks = 100)

model <- lme4::lmer(rate_not_checked ~ known_t_to_target + (1|participant), data = checks)
summary(model)

library(ggplot2)

# Plot the data with regression lines per participant
ggplot(checks, aes(x = known_t_to_target, y = rate_not_checked, color = factor(participant))) +
  geom_point() +  # Plot the raw data points
  geom_smooth(method = "lm", aes(group = participant), se = FALSE, linetype = "solid") +  # Add regression lines for each participant
  theme_minimal() +
  labs(x = "Known Time to Target", y = "Rate Not Checked", color = "Participant") +
  theme(legend.position = "none")  # Optionally hide the legend

