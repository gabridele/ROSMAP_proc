# Baseline motion-sensitivity analysis for within-network connectivity.
# Keeps the earliest session per participant and compares adjusted FD effects
# before and after the prespecified mean-FD exclusion threshold.

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

# Missingness sanity check before analysis-specific filtering.
print(colSums(is.na(demos_withinconn)))

# Keep the earliest numeric session per participant and normalize model types.
demos_min_ses <- demos_withinconn %>%
  ec_select_baseline() %>%
  ec_prepare_model_variables()

fd_threshold <- FD_THRESHOLD

# Adjusted baseline FD model.
model_formula <- within_conn ~ mean_FD + msex + site + age_scandate +
  syn_bin + eyes + dcfdx

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
      "mean_FD", "msex", "site", "age_scandate",
      "syn_bin", "eyes", "dcfdx"
    )
  )
}

fit_fd_models <- function(data_long) {
  models <- ec_fit_models_by_measure(
    data_long = data_long,
    measure_levels = target_cols,
    measure_col = "network",
    model_formula = model_formula,
    mixed = FALSE,
    optimizer = "bobyqa"
  )

  list(
    models = models,
    results = ec_fd_results_from_models(
      models = models,
      measure_col = "network",
      measure_levels = target_cols
    )
  )
}

get_adjusted_predictions <- function(models) {
  ec_fd_predictions_from_models(
    models = models,
    measure_col = "network",
    measure_levels = target_cols,
    fixed_only = FALSE
  )
}

plot_fd_effects <- function(data_long, pred_adjusted, model_results, title) {
  ec_plot_fd_effects(
    data_long = data_long,
    predictions = pred_adjusted,
    model_results = model_results,
    title = title,
    measure_col = "network",
    value_col = "within_conn",
    palette = network_colors,
    y_label = "Within-network connectivity",
    label_box = FALSE,
    label_size = 3,
    rotate_x = FALSE
  )
}

print_model_table <- function(model_results, title) {
  ec_print_fd_model_table(model_results, title)
}

# ============================================================
# 3. Pre-filtering analysis
# ============================================================

data_long_pre <- make_long(demos_min_ses)

fit_pre <- fit_fd_models(data_long_pre)
model_results_pre <- fit_pre$results
pred_adjusted_pre <- get_adjusted_predictions(fit_pre$models)

print_model_table(
  model_results_pre,
  "Pre-filtering adjusted FD model results"
)
# save table
write_csv(model_results_pre, output("motion", "fd_effects_withinconn_bl_preFDfilter_modelresults.csv"))

p_pre <- plot_fd_effects(
  data_long_pre,
  pred_adjusted_pre,
  model_results_pre,
  "Adjusted FD effect by network: baseline only"
)

print(p_pre)

# save plot
ggsave(
  output("motion", "fd_effects_withinconn_bl.png"),
  plot = p_pre,
  width = 12,
  height = 8,
  dpi = 300
)

# ============================================================
# 4. Post-filtering analysis: mean_FD < 0.25
# ============================================================

data_long_post <- data_long_pre %>%
  ec_apply_fd_filter(threshold = fd_threshold)

fit_post <- fit_fd_models(data_long_post)
model_results_post <- fit_post$results
pred_adjusted_post <- get_adjusted_predictions(fit_post$models)

print_model_table(
  model_results_post,
  paste0("Post-filtering adjusted FD model results: mean_FD < ", fd_threshold)
)

# save table
write_csv(model_results_post, output("motion", "fd_effects_withinconn_bl_postFDfilter_modelresults.csv"))

p_post <- plot_fd_effects(
  data_long_post,
  pred_adjusted_post,
  model_results_post,
  paste0("Adjusted FD effect by network: baseline only, mean_FD < ", fd_threshold)
)

print(p_post)

# save plot
ggsave(
  output("motion", "fd_effects_withinconn_bl_postFDfilter.png"),
  plot = p_post,
  width = 12,
  height = 8,
  dpi = 300
)
