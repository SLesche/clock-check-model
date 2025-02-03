library(tidyverse)
library(data.table)

# get data 
data <- fread("data_experiment2_uncleaned.csv")

# create df for pm per block 
df_pm <- data %>%
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
  # dplyr::select(participant, block_num, misc_key.started, space.rt, 
  #               misc_key.keys, section_id) %>% 
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
  mutate(known_t_to_target = 300 - cc_time) %>% 
  # # Normalize for model
  # group_by(participant, block_num) %>% 
  # mutate(rel_known_t_to_target = known_t_to_target / block_duration) %>% 
  # mutate(rel_time_since_last_cc = time_since_last_cc / known_t_to_target)
  ungroup() %>%
  mutate(subject_id = dense_rank(participant))

# write.csv(data_clean, "data_clean.csv")

data_weibull_model <- data_clean %>% 
  filter(clock_check == 1) %>% 
  select(
    participant, subject_id, block_num, known_t_to_target, cc_time, time_since_start, time_to_end, time_since_last_cc, block_duration
  ) %>% 
  mutate(cens = 0) %>% 
  filter(time_since_last_cc != 0) %>% 
  # filter(time_to_end != 0) %>% 
  na.omit() %>% 
  left_join(., df_pm)

last_cc_data <- data.frame(
  subject = rep(unique(data_weibull_model$subject_id), each = 8),
  block = rep(0:7, n = length(unique(data_weibull_model$subject_id)))
)

last_ccs <- vector(mode = "list", length = nrow(last_cc_data))

for (i in 1:nrow(last_cc_data)){
  # find last cc of this subject
  last_cc_info = data_weibull_model %>% filter(subject_id == last_cc_data[i, "subject"], block_num == last_cc_data[i, "block"]) %>% filter(cc_time == max(cc_time)) %>% head(1) 
  
  # Add a final (censored cc)
  if (nrow(last_cc_info) > 0){
    additional_cc = data.frame(
      participant = last_cc_info$participant,
      subject_id = last_cc_info$subject_id,
      block_num = last_cc_data[i, "block"],
      known_t_to_target = 300 - last_cc_info$time_since_start,
      cc_time = last_cc_info$time_since_start,
      time_since_start = last_cc_info$block_duration,
      time_to_end = 0,
      time_since_last_cc = last_cc_info$time_to_end,
      block_duration = last_cc_info$block_duration
    )
  } else {
    additional_cc = data.frame(
      participant = last_cc_data[i, "subject"],
      subject_id = last_cc_data[i, "subject"],
      block_num = last_cc_data[i, "block"],
      known_t_to_target = 300,
      cc_time = 0,
      time_since_start = data_clean %>% filter(subject_id == last_cc_data[i, "subject"], block_num == last_cc_data[i, "block"]) %>% pull(block_duration) %>% mean(na.rm = TRUE),
      time_to_end = 0,
      time_since_last_cc = data_clean %>% filter(subject_id == last_cc_data[i, "subject"], block_num == last_cc_data[i, "block"]) %>% pull(block_duration) %>% mean(na.rm = TRUE),
      block_duration = data_clean %>% filter(subject_id == last_cc_data[i, "subject"], block_num == last_cc_data[i, "block"]) %>% pull(block_duration) %>% mean(na.rm = TRUE)
    )
  }
  
  last_ccs[[i]] = additional_cc
}

last_cc_added <- data.table::rbindlist(last_ccs)
last_cc_added$cens = 1
last_cc_added <- last_cc_added %>% 
  left_join(., df_pm)

full_weibull <- data_weibull_model %>% 
  rbind(., last_cc_added) %>% 
  filter(time_since_last_cc != 0) %>%
  mutate(
    r_check = time_since_last_cc / block_duration,
    r_to_target = known_t_to_target / block_duration,
  ) %>% # filter(r > 1) %>% View()
  # mutate(
  #   r_check = ifelse(cens == 1, block_duration)
  # ) %>% 
  arrange(participant, block_num, cc_time)

write.csv(full_weibull, "weibull_data.csv")

data_diffusion_model <- data_clean %>% 
  filter(accessed_pm == 1, clock_check == 1) %>% 
  select(
    participant, block_num, known_t_to_target, cc_time, time_since_start, time_to_end, time_since_last_cc, block_duration
  ) %>% 
  filter(time_to_end != 0) %>% 
  na.omit() %>% 
  mutate(participant = dense_rank(participant))

write.csv(data_diffusion_model, "diffusion_data.csv")
