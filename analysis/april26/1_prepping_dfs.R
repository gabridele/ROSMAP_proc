library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tidyverse)
library(readxl)
library(magrittr)
library(purrr)
library(ggpubr)
library(readxl)

# ================================
#### loading of DFs
# ================================

# get todays date in DDMM format
today_date <- format(Sys.Date(), "%d%m")

ders_withage <- read_csv("sheets/derivatives_list_with_age.csv")
# make sub_ses column to be able to merge

ders_withage <- ders_withage %>%
  mutate(sub_ses = paste0(ders_withage$sub_id, "_", ders_withage$ses_id))

ders_withage <- ders_withage %>%
  select(sub_id, ses_id, sub_ses, scanner, protocol, site, age_scandate, distortion_correction, eyes)

other_demos <- read_csv("sheets/OLD_mean_within_conn_demos.csv")
other_demos <- other_demos %>%
  select(c(1:7))

other_demos <- other_demos %>%
  select(-c(site,protocol, msex))

other_demos <- other_demos %>%
  mutate(sub_ses = paste0(other_demos$sub, "_", other_demos$ses))

other_demos <- other_demos %>%
  filter(sub_ses %in% ders_withage$sub_ses)

# ignore warning about '.' its a placeholder for undisclosed values
rosmap_demos <- read_excel("sheets/ROSMAP_demos2026.xlsx")

rosmap_demos <- rosmap_demos %>%
  select(c("projid", "study", "msex", "educ", "age_bl", "dcfdx_bl"))
# filter only ROS and MAP studies

# 0 pad to 8 digits projid
rosmap_demos <- rosmap_demos %>%
  mutate(projid = str_pad(projid, width = 8, side = "left", pad = "0"))

rosmap_demos <- rosmap_demos %>% 
  mutate(projid = paste0("sub-", projid))

rosmap_demos <- rosmap_demos %>% 
  filter(projid %in% ders_withage$sub_id)

# ================================
#### checking for duplicates and counting
# ================================

# count duplicates in first column
other_demos <- other_demos %>%
  group_by(sub) %>%
  mutate(ses_count = n()) %>%
  ungroup()

# new df with unique sub column
other_demos_unique <- other_demos %>%
  group_by(sub) %>%
  slice(1) %>%
  ungroup() %>%
  select(c("sub", "ses_count"))

#plot
ggplot(other_demos_unique, aes(x = ses_count)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Distribution of Session Counts per Subject",
       x = "Number of Sessions",
       y = "Frequency") +
  theme_minimal()

# make table
ses_count_table <- other_demos_unique %>%
  group_by(ses_count) %>%
  summarise(count = n())
print(ses_count_table)

# left merge rosmap_demos and other_demos_unique by projid and sub
merged_df <- left_join(ders_withage, rosmap_demos, by = c("sub_id" = "projid"))
merged_df <- left_join(merged_df, other_demos, by = "sub_ses")

summary(merged_df)


# merge other_demos_unique and rosmap_demos by sub and projid

sexcount_df <- left_join(other_demos_unique, rosmap_demos, by = c("sub" = "projid"))
# get count of study column
study_count <- sexcount_df %>%
  group_by(study) %>%
  summarise(count = n())
print(study_count)


# ================================
#### loading of connectivity data
# ================================

demos_0426 <- read_csv("sheets/v1.3/mean_connectivity_230626.csv")

demos_connectivity <- left_join(merged_df, demos_0426, by = c("sub_ses" = "timeseries"))

# add syn_bin to demos_connectivity from column distortion_correction
demos_connectivity <- demos_connectivity %>%
  mutate(syn_bin = ifelse(distortion_correction == "SyN", 1, 0))

variables = read_excel("sheets/variables_ses_specific_may26.xlsx") %>%
  select(sub_id, ses_id, dcfdx)

# make . entry in dcfdx column to be NA
variables <- variables %>%
  mutate(dcfdx = ifelse(dcfdx == ".", NA, dcfdx))

# transform dcfdx 1 into NCI, 2 into CI, 3 into MCI, 4 into AD
variables <- variables %>%
  mutate(dcfdx = case_when(
    dcfdx == "1" ~ "NCI",
    dcfdx == "2" ~ "MCI",
    dcfdx == "3" ~ "MCI",
    dcfdx == "4" ~ "AD",
    dcfdx == "5" ~ "AD",
    dcfdx == "6" ~ "other",
    TRUE ~ as.character(dcfdx)
  ))

# merge the two dfs by sub_id and ses_id 
demos_connectivity <- demos_connectivity %>%
  left_join(variables, by = c("sub_id", "ses_id"))
demos_connectivity <- demos_connectivity %>%
  select(-c(sub, ses))

# MAKE SURE ALL COLUMNS ARE THE CORRECT TYPE
demos_connectivity <- demos_connectivity %>%
  mutate(
    sub_id = factor(sub_id),
    ses_id = factor(ses_id),
    ses_count = factor(ses_count),
    mean_FD = as.numeric(mean_FD),
    msex = factor(
      msex,
      levels = c(0, 1),
      labels = c("female", "male")
    ),
    site = factor(site),
    age_scandate = as.numeric(age_scandate),
    age_bl = as.numeric(age_bl),
    distortion_correction = factor(distortion_correction),
    eyes = factor(eyes),
    dcfdx = factor(dcfdx, levels = c("NCI", "MCI", "AD", "other")),
    syn_bin = factor(
      syn_bin,
      levels = c(0, 1),
      labels = c("not SyN", "SyN")
    )
  )

# MAKE ALL CONNECTIVITY COLUMNS NUMERIC
colnames(demos_connectivity)
connectivity_cols <- c(17:39)
demos_connectivity <- demos_connectivity %>%
  mutate(across(all_of(connectivity_cols), as.numeric))

## ADD SCANDATE
demos_connectivity <- read_csv("sheets/age_atscan.csv") %>%
  separate(col = "scandate_visit_projID", into = c("scandate", "visit", "projID"), sep = "_") %>%
  select(c("ses_id", "sub_id", "scandate")) %>%
  right_join(demos_connectivity, by = c("sub_id", "ses_id"))

# MAKE SCANDATE A DATE OBJECT
demos_connectivity <- demos_connectivity %>%
  mutate(scandate = as.Date(as.character(scandate), format = "%Y%m%d"))

# GET NUMERIC SES ID FOR CALCULATING YEARS FROM BASELINE
demos_connectivity <- demos_connectivity %>%
  mutate(
    ses_num = as.numeric(str_extract(ses_id, "\\d+"))
  ) %>%
  mutate(sub_id = factor(sub_id))

# compute YEARS FROM BASELINE for each subject
demos_connectivity <- demos_connectivity %>%
  group_by(sub_id) %>%
  mutate(
    baseline_date = scandate[which.min(ses_num)],
    years_from_baseline = interval(baseline_date, scandate) / years(1)
  ) %>%
  ungroup()

# REARRANGE COLUMNS
demos_connectivity <- demos_connectivity %>%
  select(sub_id, ses_id, sub_ses, scandate, everything())

demos_connectivity <- demos_connectivity %>%
  mutate(dcfdx = factor(dcfdx, levels = c("NCI", "MCI", "AD")))

write_csv(demos_connectivity, paste0("sheets/v1.3/demos_conn_", today_date, ".csv"))
