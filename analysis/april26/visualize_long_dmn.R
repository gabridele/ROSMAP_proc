library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggeffects)
library(lme4)
library(lmerTest)
library(lubridate)

# ===============================================================
# plotting raw data with lines connecting each subject's sessions
# ===============================================================

demos_withinconn <- read.csv("demos_withinconn.csv") %>%
  mutate(
    ses_num = as.numeric(str_extract(ses, "\\d+")),
    sub = factor(sub),
    mean_FD = as.numeric(mean_FD),
    msex = factor(msex),
    site = factor(site),
    age_scandate = as.numeric(age_scandate),
    distortion_correction = factor(distortion_correction),
    eyes = factor(eyes)
  )

demos_withinconn <- demos_withinconn %>%
  mutate(
    ses_num = as.numeric(str_extract(ses, "\\d+"))
  ) %>%
  mutate(sub = factor(sub))

## add scandate
demos_withinconn <- read_csv("age_atscan.csv") %>%
  separate(col = "scandate_visit_projID", into = c("scandate", "visit", "projID"), sep = "_") %>%
  select(c("ses_id", "sub_id", "scandate")) %>%
  right_join(demos_withinconn, by = c("sub_id", "ses_id"))

# make scandate format yyyymmdd into a datee
demos_withinconn <- demos_withinconn %>%
  mutate(scandate = as.Date(as.character(scandate), format = "%Y%m%d"))

# compute years from baseline for each subject
demos_withinconn <- demos_withinconn %>%
  group_by(sub_id) %>%
  mutate(
    baseline_date = scandate[which.min(ses_num)],
    years_from_baseline = interval(baseline_date, scandate) / years(1)
  ) %>%
  ungroup()

make_long <- function(data) {
  data %>%
    pivot_longer(
      cols = all_of(target_cols),
      names_to = "network",
      values_to = "within_conn"
    ) %>%
    mutate(network = factor(network, levels = target_cols)) %>%
    filter(
      !is.na(ses_num),
      !is.na(mean_FD),
      !is.na(within_conn),
      !is.na(msex),
      !is.na(site),
      !is.na(age_scandate),
      !is.na(distortion_correction),
      !is.na(eyes)
    )
}
data_long <- make_long(demos_withinconn)

network_to_plot <- "Default"

df_one_net <- data_long %>%
  filter(network == network_to_plot)

# ============================================================

df_segments <- df_one_net %>%
  arrange(sub, years_from_baseline) %>%
  group_by(sub) %>%
  mutate(
    x_start = years_from_baseline,
    x_end = lead(years_from_baseline),
    y_start = within_conn,
    y_end = lead(within_conn),
    direction = case_when(
      y_end > y_start ~ "up",
      y_end < y_start ~ "down",
      TRUE ~ "flat"
    )
  ) %>%
  filter(!is.na(x_end)) %>%
  ungroup()

pre_fd <- ggplot() +
  geom_segment(
    data = df_segments,
    aes(
      x = x_start,
      xend = x_end,
      y = y_start,
      yend = y_end,
      color = direction
    ),
    alpha = 0.6,
    linewidth = 0.6
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = seq(
      0,
      ceiling(max(pred_dmn$years_from_baseline, na.rm = TRUE)),
      by = 1
    ),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = c(
    "up" = "#2C7BB6",
    "down" = "#D7191C"
  )) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste(network_to_plot, "connectivity across years"),
    x = "Years from baseline",
    y = "DMN",
    color = "Step direction"
  )
print(pre_fd)
#save
ggsave(
  "dmn_withinconn_longitudinal_preFD.png",
  plot = pre_fd,
  width = 12,
  height = 8,
  dpi = 300
)

# -------
## split into 4 groups for better visualization
set.seed(123)

sub_groups <- df_one_net %>%
  distinct(sub) %>%
  mutate(group = sample(rep(1:4, length.out = nrow(.))))
df_one_net <- df_one_net %>%
  left_join(sub_groups, by = "sub")

df_segments <- df_segments %>%
  left_join(sub_groups, by = "sub")

pre_fd_split <- ggplot() +
  geom_segment(
    data = df_segments,
    aes(
      x = x_start,
      xend = x_end,
      y = y_start,
      yend = y_end,
      color = direction
    ),
    alpha = 0.6,
    linewidth = 0.6
  ) +
  
  facet_wrap(~ group) +
  
  scale_color_manual(values = c(
    "up" = "#2C7BB6",
    "down" = "#D7191C"
  )) +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = seq(
      0,
      ceiling(max(pred_dmn$years_from_baseline, na.rm = TRUE)),
      by = 1
    ),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste(network_to_plot, "connectivity across years (split view)"),
    x = "Years from baseline",
    y = "DMN",
    color = "Step direction"
  )
print(pre_fd_split)
#save

ggsave(
  "dmn_withinconn_longitudinal_preFD_split.png",
  plot = pre_fd_split,
  width = 12,
  height = 8,
  dpi = 300)

####### =======================================================
# same but after filtering for mean_FD < 0.25
####### =======================================================



data_long_post <- data_long %>%
  filter(mean_FD < 0.25)

df_one_net_postFD <- data_long_post %>%
  filter(network == network_to_plot)

df_segments <- df_one_net_postFD %>%
  arrange(sub, years_from_baseline) %>%
  group_by(sub) %>%
  mutate(
    x_start = years_from_baseline,
    x_end = lead(years_from_baseline),
    y_start = within_conn,
    y_end = lead(within_conn),
    direction = case_when(
      y_end > y_start ~ "up",
      y_end < y_start ~ "down",
      TRUE ~ "flat"
    )
  ) %>%
  filter(!is.na(x_end)) %>%
  ungroup()

post_fd <- ggplot() +
  geom_segment(
    data = df_segments,
    aes(
      x = x_start,
      xend = x_end,
      y = y_start,
      yend = y_end,
      color = direction
    ),
    alpha = 0.6,
    linewidth = 0.6
  ) +
  scale_color_manual(values = c(
    "up" = "#2C7BB6",
    "down" = "#D7191C",
    "flat" = "grey50"
  )) +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = seq(
      0,
      ceiling(max(pred_dmn$years_from_baseline, na.rm = TRUE)),
      by = 1
    ),
    expand = c(0, 0)
  ) +  theme_minimal(base_size = 13) +
  labs(
    title = paste(network_to_plot, "connectivity across years"),
    x = "Years from baseline",
    y = "DMN",
    color = "Step direction"
  )
print(post_fd)
# save
ggsave(
  "dmn_withinconn_longitudinal_postFD.png",
  plot = post_fd,
  width = 12,
  height = 8,
  dpi = 300
)

# -------
## split into 4 groups for better visualization
set.seed(123)

sub_groups <- df_one_net_postFD %>%
  distinct(sub) %>%
  mutate(group = sample(rep(1:4, length.out = nrow(.))))
df_one_net_postFD <- df_one_net_postFD %>%
  left_join(sub_groups, by = "sub")

df_segments <- df_segments %>%
  left_join(sub_groups, by = "sub") %>%
  filter(sub %in% df_one_net_postFD$sub) # keep only segments for subjects that passed the FD filter

post_fd_split <- ggplot() +
  geom_segment(
    data = df_segments,
    aes(
      x = x_start,
      xend = x_end,
      y = y_start,
      yend = y_end,
      color = direction
    ),
    alpha = 0.6,
    linewidth = 0.6
  ) +
  
  facet_wrap(~ group) +
  
  scale_color_manual(values = c(
    "up" = "#2C7BB6",
    "down" = "#D7191C",
    "flat" = "grey50"
  )) +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = seq(
      0,
      ceiling(max(pred_dmn$years_from_baseline, na.rm = TRUE)),
      by = 1
    ),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste(network_to_plot, "connectivity across years (split view)"),
    x = "Years from baseline",
    y = "DMN",
    color = "Step direction"
  )
print(post_fd_split)
# save
ggsave(
  "dmn_withinconn_longitudinal_postFD_split.png",
  plot = post_fd_split,
  width = 12,
  height = 8,
  dpi = 300
)

########################################################
# new plotting with lmer fitted lines
########################################################
# ============================================================
# with random effects 

demos_withinconn <- read.csv("demos_withinconn.csv") %>%
  mutate(
    ses_num = as.numeric(str_extract(ses, "\\d+")),
    sub = factor(sub),
    mean_FD = as.numeric(mean_FD),
    msex = factor(msex),
    site = factor(site),
    age_scandate = as.numeric(age_scandate),
    distortion_correction = factor(distortion_correction),
    eyes = factor(eyes)
  )

demos_withinconn <- demos_withinconn %>%
  mutate(
    ses_num = as.numeric(str_extract(ses, "\\d+"))
  ) %>%
  mutate(sub = factor(sub))

## add scandate
demos_withinconn <- read_csv("age_atscan.csv") %>%
  separate(col = "scandate_visit_projID", into = c("scandate", "visit", "projID"), sep = "_") %>%
  select(c("ses_id", "sub_id", "scandate")) %>%
  right_join(demos_withinconn, by = c("sub_id", "ses_id"))

# make scandate format yyyymmdd into a datee
demos_withinconn <- demos_withinconn %>%
  mutate(scandate = as.Date(as.character(scandate), format = "%Y%m%d"))

# compute years from baseline for each subject
demos_withinconn <- demos_withinconn %>%
  group_by(sub_id) %>%
  mutate(
    baseline_date = scandate[which.min(ses_num)],
    years_from_baseline = interval(baseline_date, scandate) / years(1)
  ) %>%
  ungroup()

make_long <- function(data) {
  data %>%
    pivot_longer(
      cols = all_of(target_cols),
      names_to = "network",
      values_to = "within_conn"
    ) %>%
    mutate(network = factor(network, levels = target_cols)) %>%
    filter(
      !is.na(ses_num),
      !is.na(mean_FD),
      !is.na(within_conn),
      !is.na(msex),
      !is.na(site),
      !is.na(age_scandate),
      !is.na(distortion_correction),
      !is.na(eyes)
    )
}
data_long <- make_long(demos_withinconn)

network_to_plot <- "Default"

df_one_net <- data_long %>%
  filter(network == network_to_plot)

# ============================================================
# 3. Create prediction data for each subject
# ============================================================


model_dmn <- lmer(
  within_conn ~ years_from_baseline + mean_FD + msex + site +
    age_scandate + eyes +
    (1 + years_from_baseline | sub_id),
  data = df_one_net
)

summary(model_dmn)
model_data <- model.frame(model_dmn)

# ============================================================
# 3. Build prediction data for each subject
# ============================================================
# For each subject, create a smooth sequence of session values.
# Other covariates are held at that subject's typical/observed values.

pred_dmn <- model_data %>%
  group_by(sub_id) %>%
  summarise(
    years_from_baseline = list(seq(
      min(years_from_baseline, na.rm = TRUE),
      max(years_from_baseline, na.rm = TRUE),
      length.out = 50
    )),
    mean_FD = mean(mean_FD, na.rm = TRUE),
    msex = first(as.character(msex)),
    site = first(as.character(site)),
    age_scandate = mean(age_scandate, na.rm = TRUE),
    eyes = first(as.character(eyes)),
    .groups = "drop"
  ) %>%
  unnest(years_from_baseline) %>%
  mutate(
    sub_id = factor(sub_id, levels = levels(model_data$sub_id)),
    msex = factor(msex, levels = levels(model_data$msex)),
    site = factor(site, levels = levels(model_data$site)),
    eyes = factor(eyes, levels = levels(model_data$eyes))
  )

model_data <- model.frame(model_dmn) %>%
  as_tibble()

# Check factor levels in the fitted model data
sapply(model_data, function(x) {
  if (is.factor(x)) nlevels(x) else NA
})

pred_dmn$predicted_dmn <- predict(
  model_dmn,
  newdata = pred_dmn,
  re.form = NULL,
  allow.new.levels = FALSE
)

# ============================================================
# Add subject-level direction: up vs down
# ============================================================

pred_dmn <- pred_dmn %>%
  arrange(sub_id, years_from_baseline) %>%
  group_by(sub_id) %>%
  mutate(
    fitted_delta = last(predicted_dmn) - first(predicted_dmn),
    direction = case_when(
      fitted_delta > 0 ~ "up",
      fitted_delta < 0 ~ "down",
      TRUE ~ "flat"
    )
  ) %>%
  ungroup()

# ============================================================
# Plot one model-predicted straight line per subject
# ============================================================

ggplot(
  pred_dmn,
  aes(
    x = years_from_baseline,
    y = predicted_dmn,
    group = sub_id
  )
) +
  geom_line(
    aes(color = direction),
    alpha = 0.45,
    linewidth = 0.6
  ) +
  scale_color_manual(
    values = c(
      "up" = "#2C7BB6",
      "down" = "#D7191C"
    )
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = seq(
      0,
      ceiling(max(pred_dmn$years_from_baseline, na.rm = TRUE)),
      by = 1
    ),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  labs(
    title = "Model-predicted DMN trajectories across time",
    subtitle = "One straight subject-specific fitted line per subject",
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity",
    color = "Direction"
  )
