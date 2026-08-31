library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tidyverse)
library(readxl)


rosmap_demos <- read_excel("/Users/ga0034de/Documents/rosmap_analyses_feb26/ROSMAP_demos2026.xlsx")

rosmap_demos <- rosmap_demos %>%
  select(c("projid", "study", "msex", "educ", "age_bl", "dcfdx_bl"))
# filter only ROS and MAP studies

# 0 pad to 8 digits projid
rosmap_demos <- rosmap_demos %>%
  mutate(projid = str_pad(projid, width = 8, side = "left", pad = "0"))

df_2 <- read_csv("/Users/ga0034de/Documents/rosmap_analyses_feb26/mean_within_conn_demos.csv")
# count duplicates in first column
df_2 <- df_2 %>%
  group_by(sub) %>%
  mutate(ses_count = n()) %>%
  ungroup()

# new df with unique sub column
df_2_unique <- df_2 %>%
  group_by(sub) %>%
  slice(1) %>%
  ungroup() %>%
  select(c("sub", "ses_count"))

#plot
ggplot(df_2_unique, aes(x = ses_count)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Distribution of Session Counts per Subject",
       x = "Number of Sessions",
       y = "Frequency") +
  theme_minimal()
# make table
ses_count_table <- df_2_unique %>%
  group_by(ses_count) %>%
  summarise(count = n())
print(ses_count_table)

# add sub- to projid in rosmap_demos
rosmap_demos <- rosmap_demos %>%
  mutate(projid = paste0("sub-", projid))
# left merge rosmap_demos and df_2_unique by projid and sub
merged_df <- left_join(df_2, rosmap_demos, by = c("sub" = "projid"))

# merge df_2_unique and rosmap_demos by sub and projid
sexcount_df <- left_join(df_2_unique, rosmap_demos, by = c("sub" = "projid"))
# get count of study column
study_count <- sexcount_df %>%
  group_by(study) %>%
  summarise(count = n())
print(study_count)

# get sex count
sex <- sexcount_df %>%
  group_by(msex) %>%
  summarise(count = n())
print(sex)

sex_ses <- merged_df %>%
  group_by(msex.x) %>%
  summarise(count = n())
print(sex_ses)

site <- merged_df %>%
  group_by(site) %>%
  summarise(count = n())
print(site)
# filter out FD > 0.3
merged_fd <- merged_df %>%
  filter(mean_FD <= 0.25)

sex_fd <- merged_fd %>%
  group_by(msex.x) %>%
  summarise(count = n())
print(sex_fd)

site_fd <- merged_fd %>%
  group_by(site) %>%
  summarise(count = n())
print(site_fd)

s_count <- sexcount_df %>%
  select(c("sub", "study"))

# filter out those rows NOT in ROS or MAP study
s_count <- s_count %>%
  filter(!study %in% c("ROS", "MAP"))

#save
write_csv(s_count, "projid_study_count.csv")

df3 <- read_csv("/Volumes/GabrieleSSD/backup/Desktop/R_qc_metrics/derivatives_list.csv")


# make sub_ses column by pastiung sub- and ses-
df3 <- df3 %>%
  mutate(sub_ses = paste0(sub_id, "_", ses_id))
merged_df <- left_join(merged_df, df3, by = c("sub" = "sub_id", "ses" = "ses_id"))

#get count of eyes colunm
eyes <- merged_df %>%
  group_by(eyes) %>%
  summarise(count = n())
print(eyes)

filtered <- merged_df %>%
  filter(mean_FD <= 0.25)
eyes_fd <- filtered %>%
  group_by(eyes) %>%
  summarise(count = n())
print(eyes_fd)

df_3 <- read_csv("/Volumes/GabrieleSSD/backup/Desktop/R_qc_metrics/derivatives_list.csv")

# get count of visit 00
visit_count <- df_3 %>%
  group_by(visit) %>%
  summarise(count = n())
print(visit_count)
