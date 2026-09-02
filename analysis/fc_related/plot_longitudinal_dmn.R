# Visualize longitudinal Default Mode Network (DMN) trajectories.
#
# The mixed-effects model is fitted to the Default within-network connectivity
# measure. Subject-specific (conditional) fitted trajectories are shown.

# -----------------------------------------------------------------------------
# Shared paths
# -----------------------------------------------------------------------------
.script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.script_arg)) {
  dirname(normalizePath(sub("^--file=", "", .script_arg[[1]]), mustWork = FALSE))
} else {
  getwd()
}
.paths_candidates <- unique(c(
  file.path("analysis", "paths.R"),
  file.path(.script_dir, "paths.R"),
  file.path(.script_dir, "..", "paths.R"),
  "paths.R",
  file.path("..", "paths.R")
))
.paths_file <- .paths_candidates[file.exists(.paths_candidates)][1]
if (is.na(.paths_file)) stop("Could not locate analysis/paths.R")
source(.paths_file)
rm(.script_arg, .script_dir, .paths_candidates, .paths_file)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(readr)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# Analysis conventions
# -----------------------------------------------------------------------------
NETWORK_TO_PLOT <- "Default"
MIN_VISITS <- 2L
PREDICTION_POINTS <- 50L

SEX_LEVELS <- c("female", "male")
SITE_LEVELS <- c("BNK", "UC", "MG", "RIRC")
EYES_LEVELS <- c("closed", "open")
DX_LEVELS <- c("NCI", "MCI", "AD", "other")
SYN_LEVELS <- c("not SyN", "SyN")

SITE_COLORS <- c(
  "BNK" = "#0072B2",
  "MG" = "#E69F00",
  "RIRC" = "#009E73",
  "UC" = "#D55E00"
)

DIRECTION_COLORS <- c(
  "up" = "#2C7BB6",
  "down" = "#D7191C",
  "flat" = "grey60"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
validate_columns <- function(data, columns) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "Prepared analysis table is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

factor_with_levels <- function(x, levels, variable_name) {
  observed <- unique(stats::na.omit(as.character(x)))
  unexpected <- setdiff(observed, levels)
  if (length(unexpected) > 0) {
    stop(
      "Unexpected level(s) in ", variable_name, ": ",
      paste(unexpected, collapse = ", ")
    )
  }
  factor(x, levels = levels)
}

prepare_dmn_data <- function(data, apply_fd_filter = FALSE) {
  required <- c(
    "sub_id",
    "years_from_baseline",
    "mean_FD",
    "msex",
    "site",
    "age_bl",
    "eyes",
    "dcfdx",
    "syn_bin",
    NETWORK_TO_PLOT
  )
  validate_columns(data, required)

  out <- data %>%
    transmute(
      sub_id = factor(sub_id),
      years_from_baseline = as.numeric(years_from_baseline),
      mean_FD = as.numeric(mean_FD),
      msex = factor_with_levels(msex, SEX_LEVELS, "msex"),
      site = factor_with_levels(site, SITE_LEVELS, "site"),
      age_bl = as.numeric(age_bl),
      eyes = factor_with_levels(eyes, EYES_LEVELS, "eyes"),
      dcfdx = factor_with_levels(dcfdx, DX_LEVELS, "dcfdx"),
      syn_bin = factor_with_levels(syn_bin, SYN_LEVELS, "syn_bin"),
      within_conn = as.numeric(.data[[NETWORK_TO_PLOT]])
    ) %>%
    filter(
      complete.cases(
        sub_id,
        years_from_baseline,
        mean_FD,
        msex,
        site,
        age_bl,
        eyes,
        dcfdx,
        syn_bin,
        within_conn
      )
    )

  if (apply_fd_filter) {
    out <- out %>%
      filter(mean_FD < FD_THRESHOLD)
  }

  out %>%
    arrange(sub_id, years_from_baseline) %>%
    group_by(sub_id) %>%
    filter(n() >= MIN_VISITS) %>%
    ungroup() %>%
    droplevels()
}

fit_dmn_model <- function(data, label) {
  model <- lmerTest::lmer(
    within_conn ~ years_from_baseline + mean_FD + msex + site +
      age_bl + eyes + dcfdx + syn_bin +
      (1 + years_from_baseline | sub_id),
    data = data,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000)
    )
  )

  print(summary(model))

  if (lme4::isSingular(model, tol = 1e-4)) {
    warning(label, " DMN model is singular; inspect the random-effects structure.")
  }

  model
}

build_subject_predictions <- function(model) {
  model_data <- model.frame(model) %>%
    as_tibble() %>%
    arrange(sub_id, years_from_baseline)

  pred <- model_data %>%
    group_by(sub_id) %>%
    summarise(
      years_from_baseline = list(seq(
        min(years_from_baseline, na.rm = TRUE),
        max(years_from_baseline, na.rm = TRUE),
        length.out = PREDICTION_POINTS
      )),
      mean_FD = mean(mean_FD, na.rm = TRUE),
      age_bl = mean(age_bl, na.rm = TRUE),
      msex = first(as.character(msex)),
      site = first(as.character(site)),
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
      eyes = factor(eyes, levels = levels(model_data$eyes)),
      dcfdx = factor(dcfdx, levels = levels(model_data$dcfdx)),
      syn_bin = factor(syn_bin, levels = levels(model_data$syn_bin))
    )

  # Conditional prediction: include participant random intercepts and slopes.
  pred$predicted_dmn <- predict(
    model,
    newdata = pred,
    re.form = NULL,
    allow.new.levels = FALSE
  )

  list(model_data = model_data, predictions = pred)
}

add_direction_classification <- function(predictions) {
  subject_delta <- predictions %>%
    arrange(sub_id, years_from_baseline) %>%
    group_by(sub_id) %>%
    summarise(
      fitted_delta = last(predicted_dmn) - first(predicted_dmn),
      .groups = "drop"
    )

  flat_threshold <- sd(subject_delta$fitted_delta, na.rm = TRUE)
  if (!is.finite(flat_threshold)) flat_threshold <- 0

  subject_delta <- subject_delta %>%
    mutate(
      direction = case_when(
        fitted_delta > flat_threshold ~ "up",
        fitted_delta < -flat_threshold ~ "down",
        TRUE ~ "flat"
      )
    )

  list(
    predictions = predictions %>%
      left_join(subject_delta, by = "sub_id"),
    threshold = flat_threshold,
    subject_delta = subject_delta
  )
}

subject_metadata <- function(model_data) {
  model_data %>%
    arrange(sub_id, years_from_baseline) %>%
    group_by(sub_id) %>%
    summarise(
      n_visits = n(),
      site_first = as.character(first(site)),
      n_sites = n_distinct(site),
      changed_site = if_else(n_sites > 1, "Changed site", "Same site"),
      dx_last = as.character(last(dcfdx)),
      .groups = "drop"
    )
}

save_model_coefficients <- function(model, condition) {
  coef_table <- as.data.frame(summary(model)$coefficients) %>%
    tibble::rownames_to_column("term") %>%
    mutate(condition = condition, .before = 1)

  write_csv(
    coef_table,
    output(
      "fc_related",
      paste0("dmn_model_coefficients_", condition, ".csv")
    )
  )
}

make_direction_plot <- function(predictions, threshold, title, y_limits) {
  ggplot(
    predictions,
    aes(
      x = years_from_baseline,
      y = predicted_dmn,
      group = sub_id,
      color = direction
    )
  ) +
    geom_line(alpha = 0.45, linewidth = 0.6) +
    scale_color_manual(values = DIRECTION_COLORS) +
    scale_x_continuous(
      limits = c(0, NA),
      breaks = seq(
        0,
        ceiling(max(predictions$years_from_baseline, na.rm = TRUE)),
        by = 1
      ),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = y_limits) +
    theme_classic(base_size = 13) +
    labs(
      title = title,
      subtitle = paste0(
        "Conditional subject-specific predictions. Flat = fitted change within +/- ",
        round(threshold, 4),
        " (1 SD of subject fitted changes)."
      ),
      x = "Years from baseline",
      y = "Predicted DMN within-network connectivity",
      color = "Direction"
    )
}

# -----------------------------------------------------------------------------
# Load and prepare data
# -----------------------------------------------------------------------------
demos_withinconn <- read_csv(
  require_file(demos("demos_conn.csv")),
  show_col_types = FALSE
)

pre_data <- prepare_dmn_data(
  demos_withinconn,
  apply_fd_filter = FALSE
)

post_data <- prepare_dmn_data(
  demos_withinconn,
  apply_fd_filter = TRUE
)

sample_summary <- bind_rows(
  pre_data %>%
    summarise(
      condition = "pre",
      scans = n(),
      participants = n_distinct(sub_id)
    ),
  post_data %>%
    summarise(
      condition = "post",
      scans = n(),
      participants = n_distinct(sub_id)
    )
)
print(sample_summary)
write_csv(
  sample_summary,
  output("fc_related", "dmn_longitudinal_sample_summary.csv")
)

# -----------------------------------------------------------------------------
# Fit pre- and post-FD models
# -----------------------------------------------------------------------------
model_pre <- fit_dmn_model(pre_data, "Pre-FD")
model_post <- fit_dmn_model(post_data, "Post-FD")

save_model_coefficients(model_pre, "pre_fd")
save_model_coefficients(model_post, "post_fd")

pre_fit <- build_subject_predictions(model_pre)
post_fit <- build_subject_predictions(model_post)

pre_direction <- add_direction_classification(pre_fit$predictions)
post_direction <- add_direction_classification(post_fit$predictions)

pred_pre <- pre_direction$predictions
pred_post <- post_direction$predictions

meta_pre <- subject_metadata(pre_fit$model_data)
meta_post <- subject_metadata(post_fit$model_data)

pred_pre <- pred_pre %>% left_join(meta_pre, by = "sub_id")
pred_post <- pred_post %>% left_join(meta_post, by = "sub_id")

# Save one row per subject/condition
# up/flat/down classification.
trajectory_summary <- bind_rows(
  pre_direction$subject_delta %>%
    left_join(meta_pre, by = "sub_id") %>%
    mutate(
      condition = "pre",
      flat_threshold = pre_direction$threshold,
      .before = 1
    ),
  post_direction$subject_delta %>%
    left_join(meta_post, by = "sub_id") %>%
    mutate(
      condition = "post",
      flat_threshold = post_direction$threshold,
      .before = 1
    )
)

write_csv(
  trajectory_summary,
  output("fc_related", "dmn_subject_trajectory_summary.csv")
)

# -----------------------------------------------------------------------------
# Main pre/post direction figures
# -----------------------------------------------------------------------------
y_limits_all <- range(
  c(pred_pre$predicted_dmn, pred_post$predicted_dmn),
  na.rm = TRUE
)

pre_direction_plot <- make_direction_plot(
  pred_pre,
  threshold = pre_direction$threshold,
  title = "Pre-FD-filter DMN trajectories",
  y_limits = y_limits_all
)

post_direction_plot <- make_direction_plot(
  pred_post,
  threshold = post_direction$threshold,
  title = paste0("Post-FD-filter DMN trajectories (mean FD < ", FD_THRESHOLD, ")"),
  y_limits = y_limits_all
)

print(pre_direction_plot)
print(post_direction_plot)

ggsave(
  filename = output("fc_related", "predicted_dmn_pre_fd_direction.pdf"),
  plot = pre_direction_plot,
  width = 18,
  height = 9
)

ggsave(
  filename = output("fc_related", "predicted_dmn_post_fd_direction.pdf"),
  plot = post_direction_plot,
  width = 18,
  height = 9
)

# -----------------------------------------------------------------------------
# Descriptive pre-FD trajectory plots
# -----------------------------------------------------------------------------
site_plot <- ggplot(
  pred_pre,
  aes(
    x = years_from_baseline,
    y = predicted_dmn,
    group = sub_id,
    color = site_first
  )
) +
  geom_line(alpha = 0.45, linewidth = 0.6) +
  scale_color_manual(values = SITE_COLORS, na.value = "grey80") +
  scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  theme_classic(base_size = 13) +
  labs(
    title = "Pre-FD conditional DMN trajectories by first-visit site",
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity",
    color = "First-visit site"
  )

site_change_plot <- ggplot(
  pred_pre,
  aes(
    x = years_from_baseline,
    y = predicted_dmn,
    group = sub_id,
    color = changed_site
  )
) +
  geom_line(alpha = 0.45, linewidth = 0.6) +
  facet_wrap(~site_first) +
  scale_color_manual(
    values = c(
      "Same site" = "#AD0909",
      "Changed site" = "#5281BE"
    )
  ) +
  scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  theme_classic(base_size = 13) +
  labs(
    title = "Pre-FD conditional DMN trajectories by first-visit site",
    subtitle = "Line color indicates whether acquisition site changed across fitted visits",
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity",
    color = "Site status"
  )

dx_last_plot <- ggplot(
  pred_pre,
  aes(
    x = years_from_baseline,
    y = predicted_dmn,
    group = sub_id,
    color = dx_last
  )
) +
  geom_line(alpha = 0.55, linewidth = 0.65) +
  scale_color_brewer(palette = "Dark2", na.value = "grey70") +
  scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  theme_classic(base_size = 13) +
  labs(
    title = "Pre-FD conditional DMN trajectories by last-visit diagnosis",
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity",
    color = "Last-visit diagnosis"
  )

dx_last_facet_plot <- ggplot(
  pred_pre,
  aes(
    x = years_from_baseline,
    y = predicted_dmn,
    group = sub_id
  )
) +
  geom_line(color = "#2C7BB6", alpha = 0.45, linewidth = 0.6) +
  facet_wrap(~dx_last) +
  scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  theme_classic(base_size = 13) +
  labs(
    title = "Pre-FD conditional DMN trajectories by last-visit diagnosis",
    x = "Years from baseline",
    y = "Predicted DMN within-network connectivity"
  )

exploratory_plots <- list(
  dmn_pre_fd_by_first_site = site_plot,
  dmn_pre_fd_by_site_change = site_change_plot,
  dmn_pre_fd_by_last_diagnosis = dx_last_plot,
  dmn_pre_fd_by_last_diagnosis_facet = dx_last_facet_plot
)

for (plot_name in names(exploratory_plots)) {
  print(exploratory_plots[[plot_name]])
  ggsave(
    filename = output("fc_related", paste0(plot_name, ".pdf")),
    plot = exploratory_plots[[plot_name]],
    width = 12,
    height = 8
  )
}