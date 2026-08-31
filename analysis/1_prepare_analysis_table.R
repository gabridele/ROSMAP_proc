# Prepare the analysis-ready ROSMAP connectivity table.
#
# This script merges scan metadata, restricted ROSMAP demographics,
# session-specific diagnosis, scan dates, and network-connectivity summaries.
# Private inputs are resolved through paths.R

# Shared input/output path configuration. See README.md for environment variables.
.paths_file <- if (file.exists(file.path("analysis", "paths.R"))) {
  file.path("analysis", "paths.R")
} else {
  "paths.R"
}
source(.paths_file)
rm(.paths_file)

library(tidyverse)
library(readxl)
library(lubridate)

# ================================
# Load source tables
# ================================

ders_withage <- read_csv(require_file(data("derivatives_list_with_age.csv")))
# make sub_ses column to be able to merge

ders_withage <- ders_withage %>%
  mutate(sub_ses = paste0(sub_id, "_", ses_id))

ders_withage <- ders_withage %>%
  select(sub_id, ses_id, sub_ses, scanner, protocol, site, age_scandate, distortion_correction, eyes)

other_demos <- read_csv(require_file(data("OLD_mean_within_conn_demos.csv")))
other_demos <- other_demos %>%
  select(c(1:7))

other_demos <- other_demos %>%
  select(-c(site,protocol, msex))

other_demos <- other_demos %>%
  mutate(sub_ses = paste0(sub, "_", ses))

other_demos <- other_demos %>%
  filter(sub_ses %in% ders_withage$sub_ses)

# ignore warning about '.' its a placeholder for undisclosed values
rosmap_demos <- read_excel(require_file(data("ROSMAP_demos2026.xlsx")))

rosmap_demos <- rosmap_demos %>%
  select(c("projid", "study", "msex", "educ", "age_bl", "dcfdx_bl", "dcfdx_lv"))
# filter only ROS and MAP studies

# 0 pad to 8 digits projid
rosmap_demos <- rosmap_demos %>%
  mutate(projid = str_pad(projid, width = 8, side = "left", pad = "0"))

rosmap_demos <- rosmap_demos %>% 
  mutate(projid = paste0("sub-", projid))

rosmap_demos <- rosmap_demos %>% 
  filter(projid %in% ders_withage$sub_id)

# ================================
# Duplicate/session-count checks
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

# Plot the distribution of available sessions per participant.
p_session_counts <- ggplot(other_demos_unique, aes(x = ses_count)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Distribution of Session Counts per Subject",
       x = "Number of Sessions",
       y = "Frequency") +
  theme_minimal()

ggsave(
  output("prepared", "session_count_distribution.png"),
  p_session_counts,
  width = 6,
  height = 4,
  dpi = 300
)

# Tabulate session counts.
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
# Load connectivity summaries
# ================================

# the csv file is a summary of the mean within and between connectivity for each subject and session
demos <- read_csv(require_file(data("atlas_mean_connectivity456.csv")))

# Older FC exports called the scan identifier `timeseries`; the cleaned
# compute_fc_meanconn.py writes the clearer `sub_ses` name. Accept both so
# historical tables can still be reproduced without editing the script.
if ("timeseries" %in% names(demos) && !"sub_ses" %in% names(demos)) {
  demos <- demos %>% rename(sub_ses = timeseries)
}
if (!"sub_ses" %in% names(demos)) {
  stop("Connectivity summary must contain a 'sub_ses' identifier column.")
}

demos_connectivity <- left_join(merged_df, demos, by = "sub_ses")

# add syn_bin to demos_connectivity from column distortion_correction
demos_connectivity <- demos_connectivity %>%
  mutate(syn_bin = ifelse(distortion_correction == "SyN", 1, 0))

variables <- read_excel(require_file(data("variables_ses_specific_may26.xlsx"))) %>%
  select(sub_id, ses_id, dcfdx)

# make . entry in dcfdx column to be NA
variables <- variables %>%
  mutate(dcfdx = ifelse(dcfdx == ".", NA, dcfdx))

# Recode study-specific diagnosis codes into analysis categories.
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

demos_connectivity <- demos_connectivity %>%
  mutate(dcfdx_bl = case_when(
    dcfdx_bl == "1" ~ "NCI",
    dcfdx_bl == "2" ~ "MCI",
    dcfdx_bl == "3" ~ "MCI",
    dcfdx_bl == "4" ~ "AD",
    dcfdx_bl == "5" ~ "AD",
    dcfdx_bl == "6" ~ "other",
    TRUE ~ as.character(dcfdx_bl)
  )) %>%
  mutate(dcfdx_lv = case_when(
    dcfdx_lv == "1" ~ "NCI",
    dcfdx_lv == "2" ~ "MCI",
    dcfdx_lv == "3" ~ "MCI",
    dcfdx_lv == "4" ~ "AD",
    dcfdx_lv == "5" ~ "AD",
    dcfdx_lv == "6" ~ "other",
    TRUE ~ as.character(dcfdx_lv)
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

# Convert every connectivity measure supplied by the FC summary table to numeric.
# Selecting by the source table names is robust to changes in demographic column order.
connectivity_cols <- setdiff(names(demos), "sub_ses")
demos_connectivity <- demos_connectivity %>%
  mutate(across(all_of(connectivity_cols), as.numeric))

## ADD SCANDATE
demos_connectivity <- read_csv(require_file(data("age_atscan.csv"))) %>%
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

prepared_output <- Sys.getenv(
  "ROSMAP_effect_covs_PREPARED_CSV",
  unset = output("prepared", "demos_conn.csv")
)
dir.create(dirname(prepared_output), recursive = TRUE, showWarnings = FALSE)
write_csv(demos_connectivity, prepared_output)

# print number of sub_id 
demos_connectivity %>%
  summarise(n_subs = n_distinct(sub_id), n_ses = n_distinct(ses_id)) %>%
  print()

# Derive last-visit values after sorting visits numerically within participant.
demos_connectivity <- demos_connectivity %>%
  arrange(sub_id, ses_num) %>%
  group_by(sub_id) %>%
  mutate(
    age_lv = last(age_scandate),
    dcfdx_lv = last(dcfdx)
  ) %>%
  ungroup()

summary_tbl <- demos_connectivity %>%
  arrange(sub_id, ses_id) %>%
  group_by(sub_id) %>% 
  summarise(
    age_bl = first(age_bl),
    age_lv = last(age_scandate),
    sex = first(msex),
    education = first(educ),
    dcfdx_lv = last(dcfdx),
    .groups = "drop"
  ) %>%
  summarise(
    N = n(),
    `Age BL (SD)` = sprintf("%.1f (%.1f)", mean(age_bl, na.rm = TRUE), sd(age_bl, na.rm = TRUE)),
    `Age LV (SD)` = sprintf("%.1f (%.1f)", mean(age_lv, na.rm = TRUE), sd(age_lv, na.rm = TRUE)),
    `% Female` = sprintf("%.1f", mean(sex %in% c("Female", "F", "female"), na.rm = TRUE) * 100),
    `Education (SD)` = sprintf("%.1f (%.1f)", mean(education, na.rm = TRUE), sd(education, na.rm = TRUE)),
    `%MCI LV` = sprintf("%.1f", mean(dcfdx_lv == "MCI", na.rm = TRUE) * 100),
    `%AD LV` = sprintf("%.1f", mean(dcfdx_lv == "AD", na.rm = TRUE) * 100)
  )

print(summary_tbl)
write_csv(summary_tbl, output("prepared", "sample_summary_subject_level.csv"))

summary_tbl_all_sessions <- demos_connectivity %>%
  summarise(
    N = n(),
    `Age BL (SD)` = sprintf("%.1f (%.1f)", mean(age_bl, na.rm = TRUE), sd(age_bl, na.rm = TRUE)),
    `Age LV (SD)` = sprintf("%.1f (%.1f)", mean(age_lv, na.rm = TRUE), sd(age_lv, na.rm = TRUE)),
    `% Female` = sprintf("%.1f", mean(msex %in% c("Female", "F", "female"), na.rm = TRUE) * 100),
    `Education (SD)` = sprintf("%.1f (%.1f)", mean(educ, na.rm = TRUE), sd(educ, na.rm = TRUE)),
    `%MCI LV` = sprintf("%.1f", mean(dcfdx_lv == "MCI", na.rm = TRUE) * 100),
    `%AD LV` = sprintf("%.1f", mean(dcfdx_lv == "AD", na.rm = TRUE) * 100)
  )

print(summary_tbl_all_sessions)
write_csv(
  summary_tbl_all_sessions,
  output("prepared", "sample_summary_session_level.csv")
)
