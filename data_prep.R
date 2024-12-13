library(tidyverse)
library(data.table)

# get data 
data <- fread("data_experiment2_uncleaned.csv")

data_clean <- data %>% 
  dplyr::select(participant, block_num, misc_key.started, space.rt, 
                misc_key.keys, section_id) %>% 
  filter(section_id == "M") %>%
  dplyr::rename(time = misc_key.started) %>%
  separate_rows(space.rt, convert = TRUE) %>% # sometimes more space.rt in one row
  mutate(
    clock_check = ifelse(!is.na(space.rt), 1, 0)
  ) %>% 
  group_by(participant, block_num) %>% 
  mutate(
    start = first(time),
    end = last(time),
    block_duration = end - start,
    time_since_start = time - start
  ) %>% 
  ungroup()

data_clean <- data_clean %>% 
  mutate(
    time_to_end = block_duration - time_since_start
  ) %>% 
  mutate(
    cc_time = ifelse(clock_check == 1, time_since_start, NA)
  ) %>% 
  fill(cc_time, .direction = "down") %>% 
  mutate(
    time_since_last_cc = time_since_start - cc_time
  ) %>% 
  mutate(
    time_since_last_cc = ifelse(is.na(time_since_last_cc), time_since_start, time_since_last_cc)
  ) %>% 
  mutate(time_since_last_cc = ifelse(time_since_last_cc == 0, lag(time_since_last_cc), time_since_last_cc))
