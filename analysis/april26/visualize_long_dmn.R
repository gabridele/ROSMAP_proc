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
# Load data and prepare long format
# ===============================================================

demos_withinconn <- read.csv("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/sheets/v1.3/demos_conn_2406.csv")

target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

fd_threshold <- 0.25

network_colors <- c(
  "Vis" = "#9B59B6",
  "SomMot" = "#6C8EBF",
  "Default" = "#D36B78",
  "Limbic" = "#C9D39A",
  "DorsAttn" = "#3C8D2F",
  "SalVentAttn" = "#C84CCF",
  "Cont" = "#E5B53A"
)

make_long <- function(data) {
  data %>%
    pivot_longer(
      cols = all_of(target_cols),
      names_to = "network",
      values_to = "within_conn"
    ) %>%
    mutate(
      network = factor(network, levels = target_cols)
    ) %>%
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

site_cols <- c(
  "BNK" = "#0072B2",
  "MG" = "#E69F00",
  "RIRC" = "#009E73",
  "UC" = "#d50700"
)

# ===============================================================
# LMER fitted lines: pre-FD
# ===============================================================

data_long <- make_long(demos_withinconn)

network_to_plot <- "Default"

df_one_net <- data_long %>%
  filter(network == network_to_plot)

model_dmn <- lmer(
  within_conn ~ years_from_baseline + mean_FD + msex + site +
    age_scandate + eyes + dcfdx + syn_bin +
    (1 + years_from_baseline | sub_id),
  data = df_one_net,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 100000)
  )
)

summary(model_dmn)

model_data <- model.frame(model_dmn)

# ===============================================================
# Build prediction data for each subject: pre-FD
# ===============================================================

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
    dcfdx = first(as.character(dcfdx)),
    syn_bin = first(as.character(syn_bin)),
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

# ===============================================================
# Add subject-level direction: pre-FD
# Flat = fitted delta within +/- 1 SD around zero
# ===============================================================

pred_dmn <- pred_dmn %>%
  arrange(sub_id, years_from_baseline) %>%
  group_by(sub_id) %>%
  mutate(
    fitted_delta = last(predicted_dmn) - first(predicted_dmn)
  ) %>%
  ungroup()

flat_threshold <- sd(
  pred_dmn %>%
    distinct(sub_id, fitted_delta) %>%
    pull(fitted_delta),
  na.rm = TRUE
)

pred_dmn <- pred_dmn %>%
  mutate(
    direction = case_when(
      fitted_delta > flat_threshold ~ "up",
      fitted_delta < -flat_threshold ~ "down",
      TRUE ~ "flat"
    )
  )

# ===============================================================
# LMER fitted lines: post-FD, mean_FD < 0.25
# ===============================================================

demos_withinconn_postFD <- demos_withinconn %>%
  filter(mean_FD < 0.25)

data_long_postFD <- make_long(demos_withinconn_postFD)

network_to_plot <- "Default"

df_one_net_postFD <- data_long_postFD %>%
  filter(network == network_to_plot)

df_one_net_postFD_2plus <- df_one_net_postFD %>%
  group_by(sub_id) %>%
  filter(n() >= 2) %>%
  ungroup()

df_one_net_postFD_2plus %>%
  count(sub_id) %>%
  count(n)

df_one_net_postFD_2plus %>%
  count(sub_id) %>%
  summarise(
    n_subjects = n(),
    min_obs = min(n),
    mean_obs = mean(n),
    median_obs = median(n),
    max_obs = max(n),
    n_with_1_obs = sum(n == 1),
    n_with_2plus_obs = sum(n >= 2),
    n_with_3plus_obs = sum(n >= 3)
  )

model_dmn_postFD <- lmer(
  within_conn ~ years_from_baseline + mean_FD + msex + site +
    age_scandate + eyes + dcfdx + syn_bin +
    (1 + years_from_baseline | sub_id),
  data = df_one_net_postFD_2plus,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 100000)
  )
)

summary(model_dmn_postFD)

model_data_postFD <- model.frame(model_dmn_postFD) %>%
  as_tibble()

# ===============================================================
# Build prediction data for each subject: post-FD
# ===============================================================

pred_dmn_postFD <- model_data_postFD %>%
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
    dcfdx = first(as.character(dcfdx)),
    syn_bin = first(as.character(syn_bin)),
    .groups = "drop"
  ) %>%
  unnest(years_from_baseline) %>%
  mutate(
    sub_id = factor(sub_id, levels = levels(model_data_postFD$sub_id)),
    msex = factor(msex, levels = levels(model_data_postFD$msex)),
    site = factor(site, levels = levels(model_data_postFD$site)),
    eyes = factor(eyes, levels = levels(model_data_postFD$eyes)),
    dcfdx = factor(dcfdx, levels = levels(model_data_postFD$dcfdx)),
    syn_bin = factor(syn_bin, levels = levels(model_data_postFD$syn_bin))
  )

pred_dmn_postFD$predicted_dmn <- predict(
  model_dmn_postFD,
  newdata = pred_dmn_postFD,
  re.form = NULL,
  allow.new.levels = FALSE
)

y_limits_all <- range(
  c(pred_dmn$predicted_dmn, pred_dmn_postFD$predicted_dmn),
  na.rm = TRUE
)
# ===============================================================
# Add subject-level direction: post-FD
# ===============================================================

pred_dmn_postFD <- pred_dmn_postFD %>%
  arrange(sub_id, years_from_baseline) %>%
  group_by(sub_id) %>%
  mutate(
    fitted_delta = last(predicted_dmn) - first(predicted_dmn)
  ) %>%
  ungroup()

flat_threshold_postFD <- sd(
  pred_dmn_postFD %>%
    distinct(sub_id, fitted_delta) %>%
    pull(fitted_delta),
  na.rm = TRUE
)

pred_dmn_postFD <- pred_dmn_postFD %>%
  mutate(
    direction = case_when(
      fitted_delta > flat_threshold_postFD ~ "up",
      fitted_delta < -flat_threshold_postFD ~ "down",
      TRUE ~ "flat"
    )
  )

# ===============================================================
# Plot model-predicted lines: pre-FD, colored by direction
# ===============================================================

prefd_direction <- ggplot(
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
      "down" = "#D7191C",
      "flat" = "grey60"
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
  scale_y_continuous(limits = y_limits_all) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  labs(
    title = "Model-predicted DMN trajectories across time",
    subtitle = paste0(
      "Flat = subject-specific fitted change between -",
      round(flat_threshold, 4),
      " and +",
      round(flat_threshold, 4),
      "; larger positive/negative changes labeled up/down"
    ),
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity",
    color = "Direction"
  )

ggsave(
  filename = "analysis/april26/plots/predicted_dmn_prefd_direction.pdf",
  plot = prefd_direction,
  width = 18,
  height = 9,
  dpi = 300
)
# ===============================================================
# Plot model-predicted lines: post-FD, colored by direction
# ===============================================================

postFD_fitted_plot <- ggplot(
  pred_dmn_postFD,
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
      "down" = "#D7191C",
      "flat" = "grey60"
    )
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = seq(
      0,
      ceiling(max(pred_dmn_postFD$years_from_baseline, na.rm = TRUE)),
      by = 1
    ),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = y_limits_all) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  labs(
    title = "Post-FD-filtered model-predicted DMN trajectories across time",
    subtitle = paste0(
      "Flat = subject-specific fitted change between −",
      round(flat_threshold_postFD, 4),
      " and +",
      round(flat_threshold_postFD, 4),
      "; larger positive/negative changes labeled up/down"
    ),
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity",
    color = "Direction"
  )

print(postFD_fitted_plot)

# ===============================================================
# Pre-FD fitted lines: color by site at first visit
# ===============================================================

site_first_visit <- model_data %>%
  arrange(sub_id, years_from_baseline) %>%
  group_by(sub_id) %>%
  summarise(
    site_first = first(site),
    .groups = "drop"
  )

pred_dmn <- pred_dmn %>%
  left_join(site_first_visit, by = "sub_id")

ggplot(
  pred_dmn,
  aes(
    x = years_from_baseline,
    y = predicted_dmn,
    group = sub_id
  )
) +
  geom_line(
    aes(color = site_first),
    alpha = 0.45,
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
  scale_color_manual(
    values = site_cols,
    na.value = "grey80"
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
    subtitle = "Lines colored by site at first visit",
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity",
    color = "Site at first visit"
  )

# ===============================================================
# Pre-FD fitted lines: facet by first-visit site,
# color by whether subject changed site
# ===============================================================

site_change_info <- model_data %>%
  arrange(sub_id, years_from_baseline) %>%
  group_by(sub_id) %>%
  summarise(
    site_first = first(site),
    n_sites = n_distinct(site),
    changed_site = if_else(n_sites > 1, "Changed site", "Same site"),
    .groups = "drop"
  )

pred_dmn <- pred_dmn %>%
  select(-any_of(c(
    "site_first",
    "site_first.x",
    "site_first.y",
    "n_sites",
    "changed_site",
    "changed_site.x",
    "changed_site.y"
  ))) %>%
  left_join(site_change_info, by = "sub_id")

site_facet_change_plot <- ggplot(
  pred_dmn,
  aes(
    x = years_from_baseline,
    y = predicted_dmn,
    group = sub_id
  )
) +
  geom_line(
    aes(color = changed_site),
    alpha = 0.45,
    linewidth = 0.6
  ) +
  facet_wrap(~ site_first) +
  scale_color_manual(
    values = c(
      "Same site" = "#ad0909",
      "Changed site" = "#5281be"
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
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Model-predicted DMN trajectories across time by first-visit site",
    subtitle = "Panels show site at first visit; line color shows whether subject changed site across visits",
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity",
    color = "Site status"
  )

print(site_facet_change_plot)

# ===============================================================
# Pre-FD fitted lines: color by diagnosis at last visit
# ===============================================================

dx_last_visit <- model_data %>%
  arrange(sub_id, years_from_baseline) %>%
  group_by(sub_id) %>%
  summarise(
    dx_last = last(dcfdx),
    .groups = "drop"
  )

pred_dmn <- pred_dmn %>%
  select(-any_of(c("dx_last", "dx_last.x", "dx_last.y"))) %>%
  left_join(dx_last_visit, by = "sub_id")

dx_last_plot <- ggplot(
  pred_dmn,
  aes(
    x = years_from_baseline,
    y = predicted_dmn,
    group = sub_id
  )
) +
  geom_line(
    aes(color = dx_last),
    alpha = 0.55,
    linewidth = 0.65
  ) +
  scale_color_brewer(
    palette = "Dark2",
    na.value = "grey70"
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
    axis.ticks = element_line(color = "black"),
    legend.position = "right"
  ) +
  labs(
    title = "Model-predicted DMN trajectories across time",
    subtitle = "One fitted subject-specific line per subject; color indicates diagnosis at last visit",
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity",
    color = "Diagnosis at last visit"
  )

print(dx_last_plot)

# ===============================================================
# Pre-FD fitted lines: facet by diagnosis at last visit
# ===============================================================

dx_last_facet_plot <- ggplot(
  pred_dmn,
  aes(
    x = years_from_baseline,
    y = predicted_dmn,
    group = sub_id
  )
) +
  geom_line(
    color = "#2C7BB6",
    alpha = 0.45,
    linewidth = 0.6
  ) +
  facet_wrap(~ dx_last) +
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
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Model-predicted DMN trajectories across time by diagnosis at last visit",
    subtitle = "Each panel shows subjects grouped by diagnosis at their last observed visit",
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity"
  )

print(dx_last_facet_plot)
