# Longitudinal covariate-effects analysis for between-network connectivity.
# Fits participant-level mixed models, estimated marginal means, and optional
# partial-residual plots before/after the FD exclusion threshold.

# Shared input/output path configuration. See analysis/README.md.
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

source(file.path(ANALYSIS_DIR, "effect_covs", "covs_utils.R"))
rm(.script_arg, .script_dir, .paths_candidates, .paths_file)


# ============================================================
# 1. Load and prepare data
# ============================================================

demos_betweenconn <- read.csv(require_file(demos("demos_conn.csv")))

# Missingness check before analysis-specific filtering.
print(colSums(is.na(demos_betweenconn)))

fd_threshold <- FD_THRESHOLD

# Normalize model-variable types.
demos_betweenconn <- demos_betweenconn %>%
  ec_prepare_model_variables(longitudinal = TRUE) %>%
  filter(dcfdx != "other") %>%
  droplevels()

# ============================================================
# 2. Colors for between-network combos
# ============================================================

# Make sure color vector is in the same order as target_combos
between_network_colors <- between_network_colors[target_combos]

# ============================================================
# 3. Model formula
# ============================================================

model_formula <- between_conn ~ mean_FD + msex + site + age_scandate +
  eyes + dcfdx + syn_bin + (1 | sub_id)

# ============================================================
# 2. Helper functions
# ============================================================

make_long <- function(data) {
  ec_make_long(
    data = data,
    measure_cols = target_combos,
    measure_name = "network_combo",
    value_name = "between_conn",
    required_vars = c(
      "sub_id", "mean_FD", "msex", "site", "age_scandate",
      "eyes", "dcfdx", "syn_bin"
    )
  )
}


# ============================================================
# 3. Fit adjusted models once per dataset
# ============================================================

fit_network_models <- function(data_long) {
  ec_fit_models_by_measure(
    data_long = data_long,
    measure_levels = target_combos,
    measure_col = "network_combo",
    model_formula = model_formula,
    mixed = TRUE
  )
}

# ============================================================
# 4. Add fitted values from the model list
# ============================================================

get_predicted_data <- function(data_long, models) {
  ec_get_fitted_data(
    data_long = data_long,
    models = models,
    measure_levels = target_combos,
    measure_col = "network_combo",
    fixed_only = TRUE
  )
}

# ============================================================
# 4. Model-based pairwise comparisons with emmeans
# ============================================================

fit_model_pairwise <- function(models, covariate) {
  ec_pairwise_from_models(
    models = models,
    covariate = covariate,
    measure_col = "network_combo",
    measure_levels = target_combos
  )
}

# ============================================================
# Extract partial residuals
# ============================================================

get_partial_residual_data <- function(models, covariate) {
  ec_partial_residuals_from_models(
    models = models,
    covariate = covariate,
    measure_col = "network_combo",
    measure_levels = target_combos,
    value_col = "partial_residual"
  )
}


# ============================================================
# Plot partial residual distributions
# ============================================================

plot_partial_residual_covariate <- function(
    models,
    covariate,
    pairwise_results,
    title
) {
  partial_data <- get_partial_residual_data(models, covariate)

  ec_plot_factor_distribution(
    plot_data = partial_data,
    covariate = covariate,
    pairwise_results = pairwise_results,
    title = title,
    subtitle = paste0(
      "Violin/box/jitter show partial residuals; ",
      "stars indicate BH-FDR-corrected emmeans contrasts across connectivity outcomes"
    ),
    y_label = "Partial residual",
    measure_col = "network_combo",
    value_col = "partial_residual",
    group_col = "group",
    palette = between_network_colors
  )
}

# ============================================================
# 6. Plot predicted distributions
# ============================================================

plot_factor_covariate <- function(
    predicted_data,
    covariate,
    pairwise_results,
    title,
    y_limits = NULL,
    y_breaks = NULL
) {
  ec_plot_factor_distribution(
    plot_data = predicted_data,
    covariate = covariate,
    pairwise_results = pairwise_results,
    title = title,
    subtitle = paste0(
      "Violin/box/jitter show fixed-effect predictions; ",
      "stars indicate BH-FDR-corrected emmeans contrasts across connectivity outcomes"
    ),
    y_label = "Predicted between-network connectivity",
    measure_col = "network_combo",
    value_col = "predicted_conn",
    group_col = covariate,
    palette = between_network_colors,
    y_limits = y_limits,
    y_breaks = y_breaks
  )
}

# ============================================================
# 7. Print pairwise tables
# ============================================================

print_pairwise_table <- ec_print_pairwise_table

# ============================================================
# 8. Wrapper
# ============================================================

run_factor_analysis <- function(
  predicted_data,
  models,
  covariate,
  analysis_label,
  file_suffix,
  y_limits = NULL,
  y_breaks = NULL,
  make_partial_residuals = FALSE
) {
  pairwise_results <- fit_model_pairwise(models, covariate)

  print_pairwise_table(
    pairwise_results,
    paste0(analysis_label, ": model-based pairwise comparisons for ", covariate)
  )

  p <- plot_factor_covariate(
    predicted_data = predicted_data,
    covariate = covariate,
    pairwise_results = pairwise_results,
    title = paste0(analysis_label, ": model-predicted distribution by ", covariate),
    y_limits = y_limits,
    y_breaks = y_breaks
  )

  print(p)
  ggsave(
    filename = output(
      "covariates",
      paste0("betweenconn_predicted_", covariate, "_", file_suffix, ".png")
    ),
    plot = p,
    width = 13,
    height = 9,
    dpi = 300
  )

  p_partial <- NULL
  if (make_partial_residuals) {
    p_partial <- plot_partial_residual_covariate(
      models = models,
      covariate = covariate,
      pairwise_results = pairwise_results,
      title = paste0("Partial residual distribution by ", covariate)
    )
    print(p_partial)
    ggsave(
      filename = output(
        "covariates",
        paste0("betweenconn_partial_residuals_", covariate, "_", file_suffix, ".png")
      ),
      plot = p_partial,
      width = 13,
      height = 9,
      dpi = 300
    )
  }

  list(pairwise = pairwise_results, plot = p, partial_plot = p_partial)
}

fixed_y_limits <- c(-0.2, 0.25)
fixed_y_breaks <- seq(-0.2, 0.25, by = 0.25)

# ============================================================
# 9. Pre-filtering analysis
# ============================================================

data_long_pre <- make_long(demos_betweenconn)
models_pre <- fit_network_models(data_long_pre)
predicted_pre <- get_predicted_data(data_long_pre, models_pre)

pre_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    predicted_data = predicted_pre,
    models = models_pre,
    covariate = .x,
    analysis_label = "Pre-filtering",
    file_suffix = "preFDfilter",
    y_limits = fixed_y_limits,
    y_breaks = fixed_y_breaks,
    make_partial_residuals = TRUE
  )
)

names(pre_results) <- covariates_to_run

# ============================================================
# 10. Post-filtering analysis: mean_FD < 0.25
# ============================================================

data_long_post <- data_long_pre %>%
  ec_apply_fd_filter(threshold = fd_threshold)
models_post <- fit_network_models(data_long_post)
predicted_post <- get_predicted_data(data_long_post, models_post)

post_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    predicted_data = predicted_post,
    models = models_post,
    covariate = .x,
    analysis_label = paste0("Post-filtering, mean_FD < ", fd_threshold),
    file_suffix = "postFDfilter",
    y_limits = fixed_y_limits,
    y_breaks = fixed_y_breaks
  )
)

names(post_results) <- covariates_to_run

# ============================================================
# 11. Combine all pairwise results
# ============================================================

all_pairwise_results <- bind_rows(
  map_dfr(names(pre_results), function(cov) {
    pre_results[[cov]]$pairwise %>%
      mutate(filter_status = "pre")
  }),
  map_dfr(names(post_results), function(cov) {
    post_results[[cov]]$pairwise %>%
      mutate(filter_status = "post")
  })
)

all_pairwise_results_display <- all_pairwise_results %>%
  mutate(
    estimate = round(
      estimate,
      4
    ),

    se = round(
      se,
      4
    ),

    statistic = round(
      as.numeric(statistic),
      3
    ),

    p_raw = signif(
      p_raw,
      3
    ),

    p_tukey = signif(
      p_tukey,
      3
    ),

    q_across = signif(
      q_across,
      3
    )
  ) %>%
  as_tibble()

print(all_pairwise_results_display, n = Inf)

# Preserve full numerical precision in machine-readable results.
write_csv(all_pairwise_results, output("covariates", "all_pairwise_results_betweenconn_covs_long.csv"))