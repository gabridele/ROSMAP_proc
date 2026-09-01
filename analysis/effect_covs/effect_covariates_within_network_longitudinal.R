# Longitudinal covariate-effects analysis for within-network connectivity.
# Fits participant-level mixed models and produces model-predicted covariate
# plots and pairwise comparison tables before/after FD filtering.

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

demos_withinconn <- read.csv(require_file(demos("demos_conn.csv")))

demos_withinconn <- demos_withinconn %>%
  ec_prepare_model_variables(longitudinal = TRUE) %>%
  mutate(
    dcfdx = factor(
      dcfdx,
      levels = c("NCI", "MCI", "AD", "other")
    )
  ) %>%
  filter(dcfdx != "other") %>%
  droplevels()

# Longitudinal random-intercept model
model_formula <- within_conn ~ mean_FD + msex + site + age_scandate +
  eyes + dcfdx + syn_bin + (1 | sub_id)

# ============================================================
# 2. Helper functions
# ============================================================

make_long <- function(data) {
  ec_make_long(
    data = data,
    measure_cols = target_cols,
    measure_name = "network",
    value_name = "within_conn",
    required_vars = c(
      "sub_id", "mean_FD", "msex", "site", "age_scandate",
      "eyes", "dcfdx", "syn_bin"
    )
  )
}



# ============================================================
# Output settings
# ============================================================

out_dir <- output_dir("covariates_within_longitudinal")

fixed_y_limits <- c(0, 0.6)
fixed_y_breaks <- seq(0, 0.6, by = 0.25)


# ============================================================
# Fit network models ONCE per dataset
# ============================================================

fit_network_models <- function(data_long) {
  ec_fit_models_by_measure(
    data_long = data_long,
    measure_levels = target_cols,
    measure_col = "network",
    model_formula = model_formula,
    mixed = TRUE
  )
}


# ============================================================
# Get model-predicted values
# ============================================================

get_predicted_data <- function(data_long, models) {
  ec_get_fitted_data(
    data_long = data_long,
    models = models,
    measure_levels = target_cols,
    measure_col = "network",
    fixed_only = TRUE
  )
}


# ============================================================
# Pairwise comparisons
# ============================================================

fit_model_pairwise <- function(models, covariate) {
  ec_pairwise_from_models(
    models = models,
    covariate = covariate,
    measure_col = "network",
    measure_levels = target_cols
  )
}


# ============================================================
# Extract partial residuals
# ============================================================

get_partial_residual_data <- function(models, covariate) {
  ec_partial_residuals_from_models(
    models = models,
    covariate = covariate,
    measure_col = "network",
    measure_levels = target_cols,
    value_col = "value"
  )
}


# ============================================================
# Common plot function
# ============================================================

plot_factor_distribution <- function(
    plot_data,
    covariate,
    pairwise_results,
    title,
    subtitle,
    y_lab,
    y_limits = NULL,
    y_breaks = NULL
) {
  ec_plot_factor_distribution(
    plot_data = plot_data,
    covariate = covariate,
    pairwise_results = pairwise_results,
    title = title,
    subtitle = subtitle,
    y_label = y_lab,
    measure_col = "network",
    value_col = "value",
    group_col = "group",
    palette = network_colors,
    y_limits = y_limits,
    y_breaks = y_breaks,
    facet_ncol = 3,
    facet_axes_all = TRUE,
    x_text_angle = 0,
    x_text_size = 20,
    y_text_size = 20
  )
}

# ============================================================
# Save PDF + editable SVG
# ============================================================

save_plot <- function(plot, filename, width = 13, height = 7) {
  ec_save_pdf_svg(
    plot = plot,
    out_dir = out_dir,
    filename = filename,
    width = width,
    height = height
  )
}

# ============================================================
# Run one covariate
# ============================================================

run_factor_analysis <- function(
  predicted_data,
  models,
  covariate,
  analysis_label,
  file_suffix
) {

  # ----------------------------------------------------------
  # Pairwise statistics
  # ----------------------------------------------------------

  pairwise_results <- fit_model_pairwise(
    models,
    covariate
  )

  print(
    pairwise_results
  )


  # ----------------------------------------------------------
  # Predicted-value data
  # ----------------------------------------------------------

  pred_plot_data <- predicted_data %>%
    transmute(
      network,
      group = .data[[covariate]],
      value = predicted_conn
    )


  # ----------------------------------------------------------
  # Predicted-value plot
  # ----------------------------------------------------------

  p_predicted <- plot_factor_distribution(
    plot_data = pred_plot_data,
    covariate = covariate,
    pairwise_results = pairwise_results,

    title = paste0(
      analysis_label,
      ": model-predicted distribution by ",
      covariate
    ),

    subtitle =
      "Violin/box/jitter show fixed-effect predictions; stars indicate BH-FDR-corrected emmeans contrasts across connectivity outcomes",

    y_lab =
      "Predicted within-network connectivity",

    y_limits = fixed_y_limits,
    y_breaks = fixed_y_breaks
  )


  # ----------------------------------------------------------
  # Partial residual data
  # ----------------------------------------------------------

  partial_plot_data <- get_partial_residual_data(
    models,
    covariate
  )


  # ----------------------------------------------------------
  # Partial residual plot
  # ----------------------------------------------------------

  p_partial <- plot_factor_distribution(
    plot_data = partial_plot_data,
    covariate = covariate,
    pairwise_results = pairwise_results,

    title = paste0(
      analysis_label,
      ": partial residual distribution by ",
      covariate
    ),

    subtitle =
      "Violin/box/jitter show partial residuals; stars indicate BH-FDR-corrected emmeans contrasts across connectivity outcomes",

    y_lab =
      "Within-network connectivity (partial residual)",

    # DON'T force partial residuals to 0–0.6
    y_limits = NULL,
    y_breaks = NULL
  )


  # ----------------------------------------------------------
  # Save PDF + SVG
  # ----------------------------------------------------------

  save_plot(
    p_predicted,
    paste0(
      "withinconn_predicted_",
      covariate,
      "_",
      file_suffix
    )
  )

  save_plot(
    p_partial,
    paste0(
      "withinconn_partial_residuals_",
      covariate,
      "_",
      file_suffix
    )
  )


  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

  list(
    pairwise = pairwise_results,
    predicted_plot = p_predicted,
    partial_plot = p_partial
  )
}

# ============================================================
# PRE-FILTER
# ============================================================

data_long_pre <- make_long(
  demos_withinconn
)

models_pre <- fit_network_models(
  data_long_pre
)

predicted_pre <- get_predicted_data(
  data_long_pre,
  models_pre
)

pre_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    predicted_data = predicted_pre,
    models = models_pre,
    covariate = .x,
    analysis_label = "Pre-filtering",
    file_suffix = "preFDfilter"
  )
) %>%
  set_names(covariates_to_run)

# ============================================================
# POST-FILTER
# ============================================================

data_long_post <- data_long_pre %>%
  ec_apply_fd_filter(threshold = fd_threshold)

models_post <- fit_network_models(
  data_long_post
)

predicted_post <- get_predicted_data(
  data_long_post,
  models_post
)

post_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    predicted_data = predicted_post,
    models = models_post,
    covariate = .x,
    analysis_label = paste0(
      "Post-filtering, mean_FD < ",
      fd_threshold
    ),
    file_suffix = "postFDfilter"
  )
) %>%
  set_names(covariates_to_run)


# ============================================================
# Combine pairwise results
# ============================================================

all_pairwise_results <- bind_rows(

  map_dfr(
    pre_results,
    ~ .x$pairwise %>%
      mutate(filter_status = "pre")
  ),

  map_dfr(
    post_results,
    ~ .x$pairwise %>%
      mutate(filter_status = "post")
  )
) %>%
  as_tibble()

print(
  all_pairwise_results,
  n = Inf
)

write_csv(
  all_pairwise_results,
  file.path(
    out_dir,
    "all_pairwise_results_withinconn_covs_long.csv"
  )
)