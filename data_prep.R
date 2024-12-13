library(tidyverse)
library(data.table)

# get data 
data <- fread("data_experiment2_uncleaned.csv")

# create df for pm per block 
df_pm <- df_full %>%
  filter(section_id == "M") %>%
  dplyr::rename(time = misc_key.started) %>% 
  mutate(block_num = case_when(section_id == "M" ~ block_num + 1,
                               TRUE ~ block_num)) %>%
  dplyr::select(participant, block_num, accessed_pm, start_time, time) %>% 
  # accessed pm is written in first line of next block and needs to be moved up
  mutate_at(c("accessed_pm"), list(lead), n = 1 ) %>% 
  group_by(participant, block_num) %>%
  # make new variable with last entry for each block, needed for acc
  mutate(time_last = last(time[!is.na(time)])) %>% 
  # calculate absolute value from 5 min for each block 
  mutate(pm_acc = ifelse(block_num == 1 & accessed_pm == 1, 
                         abs(last(time_last) - start_time - 300),
                  ifelse(block_num == 2 & accessed_pm == 1, 
                          abs(last(time_last) - start_time - 600), 
                  ifelse(block_num == 3 & accessed_pm == 1, 
                          abs(last(time_last) - start_time - 900),
                  ifelse(block_num == 4 & accessed_pm == 1, 
                          abs(last(time_last) - start_time - 1200),
                  ifelse(block_num == 5 & accessed_pm == 1, 
                          abs(last(time_last) - start_time - 1500), 
                  ifelse(block_num == 6 & accessed_pm == 1, 
                          abs(last(time_last) - start_time - 1800), 
                  ifelse(block_num == 7 & accessed_pm == 1, 
                          abs(last(time_last) - start_time - 2100), 
                  ifelse(block_num == 8 & accessed_pm == 1, 
                          abs(last(time_last) - start_time - 2400), NA))))))))) %>%
  # drop empty lines (double rows for blocks otherwise)
  drop_na(accessed_pm) %>% 
  unique() %>% 
  # make pm_acc NA if time is > 30 sec (response interval)
  mutate(pm_acc = case_when(pm_acc > 30 ~ NA_real_, TRUE ~ pm_acc)) %>% 
  # make accessed_pm 0 if pm_acc is NA 
  mutate(accessed_pm = case_when(is.na(pm_acc) ~ 0, TRUE ~ accessed_pm)) %>% 
  group_by(participant) %>% 
  # make pm_count depended on pm_acc (can be adjusted)
  mutate(pm_count = sum(pm_acc < 30 ,na.rm=TRUE)) %>%
  ungroup() %>% 
  dplyr::select(participant, block_num, accessed_pm, pm_acc, pm_count) %>% 
  mutate(block_num = block_num -1)

data_clean <- data %>% 
  dplyr::select(participant, block_num, misc_key.started, space.rt, 
                misc_key.keys, section_id) %>% 
  left_join(., df_pm) %>% 
  filter(section_id == "M") %>%
  dplyr::rename(time = misc_key.started) %>%
  filter(!is.na(time)) %>% 
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
  group_by(participant, block_num) %>% 
  fill(cc_time, .direction = "down") %>% 
  mutate(cc_time = ifelse(is.na(cc_time), 0, cc_time)) %>% 
  mutate(
    cc_time = ifelse(clock_check == 1, lag(cc_time), cc_time)) %>% 
  ungroup() %>%  
  mutate(
    time_since_last_cc = time_since_start - cc_time
  ) %>% 
  mutate(known_t_to_target = block_duration - cc_time) %>% 
  # # Normalize for model
  # group_by(participant, block_num) %>% 
  # mutate(rel_known_t_to_target = known_t_to_target / block_duration) %>% 
  # mutate(rel_time_since_last_cc = time_since_last_cc / known_t_to_target)
  ungroup()

data_model <- data_clean %>% 
  filter(accessed_pm == 1) %>% 
  rename(
    "id" = "participant",
    "t_since_last_check" = time_since_last_cc,
    "check" = clock_check
  ) %>% 
  select(
    id, t_since_last_check, known_t_to_target, check
  ) %>% 
  na.omit() %>% 
  mutate(id = dense_rank(id))

write.csv(data_model, "cleaned_data.csv")
