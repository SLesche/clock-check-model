# data cleaning pmp 
# strg shift c = #
# strg shift m = %>% 

library(tidyverse)
library(data.table)

# preparations to clean R before starting
rm(list = ls())
# which packages
search()

# get data
path <- getwd()

# get data 
df_full <- read.csv("data_experiment2_uncleaned.csv")

df_ot_acc <- df_full %>%
  # fill so that they are not deleted later 
  fill(pm_difficulty.response, .direction = "up") %>%
  fill(reasons_missing.response, .direction = "up") %>%
  # filters practice blocks
  filter(block_duration == 300) %>%
  dplyr::select(participant, block_num, section_id, correct, answer.rt, 
                misc_key.started, return.rt, space.rt, 
                misc_key.keys, accessed_pm, pm_difficulty.response, 
                reasons_missing.response) %>%
  # make stuff more handy 
  dplyr::rename(pm_diff = pm_difficulty.response) %>% 
  dplyr::rename(time = misc_key.started) %>% 
  dplyr::rename(why_miss = reasons_missing.response) %>% 
  # baseline block is block_num = 0 
  mutate(block_num = case_when(section_id == "M" ~ block_num + 1, 
                               TRUE ~ block_num)) %>% 
  mutate(RT_ms = answer.rt*1000) %>%
  # filter reaction times below 150 ms 
  filter(RT_ms>150) %>%
  # make all NAs (no response) 3000 ms (max)
  mutate(RT_ms = ifelse(is.na(RT_ms) & correct == 0, 3000, RT_ms)) %>%
  filter(!is.na(RT_ms)) %>% 
  # calculate total rts to cut later 
  group_by(participant) %>% 
  mutate(ot_rt_total = mean(RT_ms),
         SD_ot_rt_total = sd(RT_ms), CutLow = ot_rt_total - 2.5*SD_ot_rt_total,
         CutHigh = ot_rt_total + 2.5*SD_ot_rt_total)%>%
  # filter RTs +/- 2,5 SD
  filter(RT_ms > CutLow, RT_ms < CutHigh) %>%
  # calculate mean rt for participant
  mutate(ot_rt_mean = mean(RT_ms)) %>% 
  # calculate mean rt within blocks a z-standardise values 
  group_by(participant, block_num) %>% 
  mutate(ot_rt_block = mean(RT_ms),
         ot_rt_block_z = ot_rt_block - ot_rt_mean) %>% 
  # calculate accuracy total
  group_by(participant) %>%
  mutate(ot_acc_tot = mean(correct, na.rm = TRUE)) %>% 
  # calculate accuracy per block
  group_by(participant, block_num) %>% 
  mutate(acc_block = mean(correct, na.rm = TRUE)) %>% 
  # divide block into 4 timepoints 
  group_by(participant, block_num) %>% 
  mutate(
    start = first(time),
    end = last(time),
    block_duration = end - start,
    time_since_start = time - start,
    timepoint = case_when(
      time_since_start >= 0 & time_since_start < block_duration/4 ~ "T1",
      time_since_start >= block_duration/4 & time_since_start < 2 * block_duration/4 ~ "T2",
      time_since_start >= 2 * block_duration/4 & time_since_start < 3 * block_duration/4 ~ "T3",
      TRUE ~ "T4")) %>%
  # calculate mean RT for each timepoint 
  group_by(participant, block_num, timepoint) %>% 
  mutate(ot_rt_T = mean(RT_ms), 
         ot_acc_T = mean(correct)) %>% 
  # calculate rts in T4
  group_by(participant, block_num) %>% 
  # + value = slower towards end of block, - value = fast towards end 
  mutate(ot_rt_T4 = case_when(timepoint == "T4" ~ ot_rt_T / ot_rt_block),
         ot_acc_T4 = case_when(timepoint == "T4" ~ ot_acc_T / acc_block)) %>%  
  fill(ot_rt_T4, .direction = "up") %>%
  fill(ot_acc_T4, .direction = "up") %>% 
  dplyr::select(participant, block_num, ot_rt_mean, ot_rt_block, ot_rt_block_z, 
                timepoint, ot_rt_T, ot_rt_T4, ot_acc_tot, acc_block,
                ot_acc_T4, pm_diff, why_miss) %>% 
  unique()

df_ot_rt <- df_full %>%
  # filters practice blocks
  filter(block_duration == 300) %>%
  dplyr::select(participant, block_num, section_id, answer.rt, correct, 
                misc_key.started) %>%
  # make stuff more handy 
  dplyr::rename(time = misc_key.started) %>% 
  # baseline block is block_num = 0 
  mutate(block_num = case_when(section_id == "M" ~ block_num + 1, 
                               TRUE ~ block_num)) %>% 
  # divide block into 4 timepoints 
  group_by(participant, block_num) %>% 
  mutate(
    start = first(time),
    end = last(time),
    block_duration = end - start,
    time_since_start = time - start,
    timepoint = case_when(
      time_since_start >= 0 & time_since_start < block_duration/4 ~ "T1",
      time_since_start >= block_duration/4 & time_since_start < 2 * block_duration/4 ~ "T2",
      time_since_start >= 2 * block_duration/4 & time_since_start < 3 * block_duration/4 ~ "T3",
      TRUE ~ "T4")) %>%
  mutate(RT_ms = answer.rt*1000) %>%
  # filter reaction times below 150 ms 
  filter(RT_ms>150) %>%
  # make all NAs (no response) 3000 ms (max)
  mutate(RT_ms = ifelse(is.na(RT_ms) & correct == 0, 3000, RT_ms)) %>%
  filter(!is.na(RT_ms)) %>% 
  filter(correct == 1) %>% 
  # calculate total rts to cut later 
  group_by(participant) %>% 
  mutate(ot_rt_total = mean(RT_ms),
         SD_ot_rt_total = sd(RT_ms), CutLow = ot_rt_total - 2.5*SD_ot_rt_total,
         CutHigh = ot_rt_total + 2.5*SD_ot_rt_total)%>%
  # filter RTs +/- 2,5 SD
  filter(RT_ms > CutLow, RT_ms < CutHigh) %>%
  # calculate mean rt for participant
  mutate(ot_rt_cor_mean = mean(RT_ms)) %>% 
  # calculate mean rt for participant in blocks 
  group_by(participant, block_num) %>% 
  mutate(ot_rt_cor_block = mean(RT_ms)) %>%  
  # calculate mean RT for each timepoint 
  group_by(participant, block_num, timepoint) %>% 
  mutate(ot_rt_cor_T = mean(RT_ms)) %>% 
  # calculate rts in T4
  group_by(participant, block_num) %>% 
  # + value = slower towards end of block, - value = fast towards end 
  mutate(ot_rt_cor_T4 = case_when(timepoint == "T4" ~ ot_rt_cor_T / ot_rt_cor_block)) %>% # 
  fill(ot_rt_cor_T4, .direction = "up") %>%
  dplyr::select(participant, block_num, ot_rt_cor_mean, ot_rt_cor_block,  
                timepoint, ot_rt_cor_T, ot_rt_cor_T4) %>% 
  unique()

df_new_1 <- df_ot_acc %>% 
  dplyr::select(participant, block_num, 
                ot_rt_mean, ot_rt_block, ot_rt_block_z, ot_rt_T4, 
                ot_acc_tot, acc_block, ot_acc_T4, 
                pm_diff, why_miss) %>%  
  unique()

df_new_2 <- df_ot_rt %>% 
  dplyr::select(participant, block_num, 
                ot_rt_cor_mean, ot_rt_cor_block, ot_rt_cor_T4) %>%  
  unique()

# create df for clock checks per participant and block 
df_cc <- df_full %>%  
  dplyr::select(participant, block_num, misc_key.started, space.rt, 
         misc_key.keys, section_id) %>% 
  filter(section_id == "M") %>%
  # values have to be moved one space up, otherwise wrong block 
  mutate(block_num = case_when(section_id == "M" ~ block_num + 1, 
                               TRUE ~ block_num)) %>% 
  dplyr::rename(time = misc_key.started) %>% 
  # create 4 timepoints within each block
  group_by(participant, block_num) %>% 
  mutate(
    start = first(time),
    end = last(time),
    block_duration = end - start,
    time_since_start = time - start,
    timepoint = case_when(
      time_since_start >= 0 & time_since_start < block_duration/4 ~ "T1",
      time_since_start >= block_duration/4 & time_since_start < 2 * block_duration/4 ~ "T2",
      time_since_start >= 2 * block_duration/4 & time_since_start < 3 * block_duration/4 ~ "T3",
      TRUE ~ "T4")) %>%
  ungroup() %>% 
  separate_rows(space.rt, convert = TRUE) %>% # sometimes more space.rt in one row
  group_by(participant, block_num, timepoint) %>% 
  # calulate clock checks per timepoint 
  mutate(cc_timepoint = sum(space.rt > 0, na.rm = TRUE)) %>% # cc / tp 
  dplyr::select(participant, block_num, timepoint, cc_timepoint) %>%
  unique() %>% 
  ungroup() %>% 
  # get complete number of clock checks per block 
  group_by(participant, block_num) %>% 
  mutate(cc_count = sum(cc_timepoint)) %>% # cc / block total 
  # strategy of clockchecking dividing ccs in T4 by ccs in the block (per block)
  mutate(cc_strat = ifelse(timepoint == "T4", cc_timepoint / cc_count, NA) * 100) %>% 
  fill(cc_strat, .direction = "up") %>% 
  dplyr::select(participant, block_num, cc_count, cc_strat) %>% 
  unique() %>%
  ungroup() %>% 
  # if no clock check occurs in block strat = 0 
  mutate(cc_strat = case_when(cc_count == 0 ~ 0, TRUE ~ cc_strat)) %>% 
  group_by(participant) %>% 
  # measure of strategy for each participant 
  mutate(strat_all = sum(cc_strat, na.rm = TRUE)/ 8)

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
  dplyr::select(participant, block_num, accessed_pm, pm_acc, pm_count)

# create df for all relevant data for analyses
# 8 rows per participant to have all relevant data
df_pmp <- df_new_1 %>% 
  left_join(df_new_2) %>% 
  left_join(df_cc) %>% 
  left_join(df_pm) %>% 
  arrange(participant)

length(unique(df_pmp$participant))

df_slowing <- df_ot_rt %>% 
  filter(block_num != 0) %>% 
  left_join(df_pm) %>%  
  dplyr::select(participant, timepoint, ot_rt_cor_T, accessed_pm) %>% 
  group_by(timepoint, accessed_pm) %>% 
  mutate(T_cor_mean = mean(ot_rt_cor_T),
         T_cor_sd = sd(ot_rt_cor_T), 
         T_cor_se = sd(ot_rt_cor_T) / sqrt(n()))
  
# save in folder for analyses 
path_out = 'where you want to save it'
write.csv(df_pmp, file.path(path_out, "pm_task"), row.names=FALSE)
write.csv(df_slowing, file.path(path_out, "slowing_exp2"), row.names=FALSE)

