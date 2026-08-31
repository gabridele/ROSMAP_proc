# Longitudinal motion-sensitivity analysis for within-network connectivity.
# Fits mixed-effects models before and after FD exclusion and generates the
# model tables and diagnostic figures used to assess motion bias.

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
rm(.script_arg, .script_dir, .paths_candidates, .paths_file)

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggeffects)
library(lme4)
library(lmerTest)
library(visreg)
library(patchwork)
library(svglite)


# ============================================================
# 1. Load and prepare data
# ============================================================

demos_withinconn <- read.csv(
  require_file(demos("demos_conn.csv"))
)

demos_withinconn <- demos_withinconn %>%
  mutate(
    ses_num = as.numeric(str_extract(ses_id, "\\d+")),
    sub_id = factor(sub_id)
  )


# Networks to analyse
target_cols <- c(
  "Vis",
  "SomMot",
  "DorsAttn",
  "SalVentAttn",
  "Limbic",
  "Cont",
  "Default"
)


# FD exclusion threshold
fd_threshold <- FD_THRESHOLD

# Network colours
network_colors <- c(
  "Vis" = "#9B59B6",
  "SomMot" = "#6C8EBF",
  "Default" = "#D36B78",
  "Limbic" = "#C9D39A",
  "DorsAttn" = "#3C8D2F",
  "SalVentAttn" = "#C84CCF",
  "Cont" = "#E5B53A"
)


# ============================================================
# 2. Convert network data to long format
# ============================================================

make_long <- function(data) {

  data %>%
    pivot_longer(
      cols = all_of(target_cols),
      names_to = "network",
      values_to = "within_conn"
    ) %>%
    mutate(
      network = factor(
        network,
        levels = target_cols
      )
    ) %>%
    filter(
      !is.na(sub_id),
      !is.na(ses_num),
      !is.na(mean_FD),
      !is.na(within_conn),
      !is.na(msex),
      !is.na(site),
      !is.na(age_scandate),
      !is.na(eyes),
      !is.na(syn_bin),
      !is.na(dcfdx)
    )
}


# ============================================================
# 3. Create ALL and FD-FILTERED datasets
# ============================================================

data_long_all <- make_long(
  demos_withinconn
)


# Keep observations with mean FD < threshold
data_long_fd_filtered <- data_long_all %>%
  filter(
    mean_FD < fd_threshold
  )


# ============================================================
# 4. Print sample-size information
# ============================================================

cat(
  "\n============================================================\n",
  "SAMPLE SIZE BEFORE AND AFTER FD FILTERING\n",
  "============================================================\n",
  sep = ""
)


# Because data are in long format, count unique subject-session
# combinations rather than network rows.

sample_summary <- bind_rows(

  data_long_all %>%
    distinct(sub_id, ses_num) %>%
    summarise(
      dataset = "All observations",
      n_subjects = n_distinct(sub_id),
      n_scans = n()
    ),

  data_long_fd_filtered %>%
    distinct(sub_id, ses_num) %>%
    summarise(
      dataset = paste0(
        "FD < ",
        fd_threshold
      ),
      n_subjects = n_distinct(sub_id),
      n_scans = n()
    )
)

print(sample_summary)


# Number of unique scans excluded
n_scans_all <- data_long_all %>%
  distinct(sub_id, ses_num) %>%
  nrow()

n_scans_filtered <- data_long_fd_filtered %>%
  distinct(sub_id, ses_num) %>%
  nrow()

cat(
  "\nScans excluded by FD threshold:",
  n_scans_all - n_scans_filtered,
  "\n"
)


# ============================================================
# 5. Global plotting / output settings
# ============================================================

out_dir <- output_dir("fd_within_longitudinal")


# Same model is fitted in both datasets.
#
# ALL:
#   Uses the full available FD range.
#
# FD FILTERED:
#   Uses only observations with mean_FD < 0.25.
#
# mean_FD remains the predictor in both models.

fd_formula <- within_conn ~
  mean_FD +
  ses_num +
  msex +
  site +
  age_scandate +
  eyes +
  syn_bin +
  dcfdx +
  (1 | sub_id)


# Keep same y-axis across figures for direct comparison
fixed_y_limits <- c(-0.2, 0.7)

fixed_y_breaks <- seq(
  -0.2,
  0.8,
  by = 0.2
)


# ============================================================
# 6. Fit network models + extract statistics
# ============================================================

fit_fd_models <- function(data_long) {

  models <- map(
    target_cols,
    function(net) {

      df_net <- data_long %>%
        filter(
          network == net
        )

      lmer(
        fd_formula,
        data = df_net,
        control = lmerControl(
          optimizer = "bobyqa"
        )
      )
    }
  ) %>%
    set_names(
      target_cols
    )


  # ----------------------------------------------------------
  # Extract FD coefficient from each network model
  # ----------------------------------------------------------

  results <- imap_dfr(
    models,
    function(model, net) {

      coefs <- summary(model)$coefficients

      tibble(
        network = net,

        beta_adjusted = coefs[
          "mean_FD",
          "Estimate"
        ],

        t_val_adjusted = coefs[
          "mean_FD",
          "t value"
        ],

        p_adjusted = coefs[
          "mean_FD",
          "Pr(>|t|)"
        ]
      )
    }
  ) %>%

    mutate(

      # BH/FDR correction across the 7 networks
      q_adjusted = p.adjust(
        p_adjusted,
        method = "BH"
      ),

      sig_adjusted = case_when(
        q_adjusted < 0.001 ~ "***",
        q_adjusted < 0.01  ~ "**",
        q_adjusted < 0.05  ~ "*",
        TRUE ~ ""
      ),

      label = sprintf(
        "β = %.3f %s\nq = %.3g",
        beta_adjusted,
        sig_adjusted,
        q_adjusted
      ),

      network = factor(
        network,
        levels = target_cols
      )
    )


  list(
    models = models,
    results = results
  )
}


# ============================================================
# 7. Fixed-effect predictions
# ============================================================

get_adjusted_predictions <- function(models) {

  imap_dfr(
    models,
    function(model, net) {

      predict_response(
        model,
        terms = "mean_FD [all]",
        type = "fixed"
      ) %>%
        as.data.frame() %>%
        mutate(
          network = net
        )
    }
  ) %>%

    mutate(
      network = factor(
        network,
        levels = target_cols
      )
    )
}


# ============================================================
# 8. Common figure theme
#    Readable at manuscript-figure scale
# ============================================================

network_plot_theme <- function() {

  theme_minimal(
    base_size = 17
  ) +

    theme(

      # No legend
      legend.position = "none",


      # Network names above facets
      strip.text = element_text(
        size = 16,
        face = "bold"
      ),


      # Main plot title
      plot.title = element_text(
        size = 19,
        face = "bold",
        margin = margin(
          b = 10
        )
      ),


      # Axis titles
      axis.title.x = element_text(
        size = 17,
        face = "bold",
        margin = margin(
          t = 8
        )
      ),

      axis.title.y = element_text(
        size = 17,
        face = "bold",
        margin = margin(
          r = 8
        )
      ),


      # Tick labels
      # Original script used size = 9.
      # Keep axis text legible in exported figures.
      axis.text.x = element_text(
        size = 14
      ),

      axis.text.y = element_text(
        size = 14
      ),


      # Remove minor grid
      panel.grid.minor = element_blank()
    )
}


# ============================================================
# 9. Fitted-value plot
# ============================================================

plot_fd_effects <- function(
  data_long,
  predictions,
  model_results,
  title,
  y_limits = fixed_y_limits,
  y_breaks = fixed_y_breaks
) {

  labels <- model_results %>%
    mutate(
      x = Inf,
      y = Inf
    )


  ggplot(
    data_long,
    aes(
      x = mean_FD,
      y = within_conn,
      color = network
    )
  ) +

    # Raw observations
    geom_point(
      alpha = 0.3,
      size = 1
    ) +


    # Confidence interval
    geom_ribbon(
      data = predictions,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = network
      ),
      inherit.aes = FALSE,
      alpha = 0.25
    ) +


    # Adjusted fixed-effect line
    geom_line(
      data = predictions,
      aes(
        x = x,
        y = predicted,
        color = network
      ),
      inherit.aes = FALSE,
      linewidth = 1.2
    ) +


    # Model statistics
    geom_label(
      data = labels,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 5,
      color = "black",
      fill = "white",
      alpha = 0.8,
      linewidth = 0
    ) +


    # Axes + tick labels on every panel
    facet_wrap(
      ~ network,
      ncol = 3,
      scales = "fixed",
      axes = "all",
      axis.labels = "all"
    ) +


    scale_color_manual(
      values = network_colors
    ) +

    scale_fill_manual(
      values = network_colors
    ) +


    scale_y_continuous(
      breaks = y_breaks
    ) +


    coord_cartesian(
      ylim = y_limits
    ) +


    network_plot_theme() +


    labs(
      title = title,
      x = "Mean FD",
      y = "Within-network connectivity"
    )
}


# ============================================================
# 10. Get visreg partial residual data
# ============================================================

get_partial_residual_data <- function(
  models,
  predictor = "mean_FD"
) {

  # Run visreg only once per model
  visreg_objects <- imap(
    models,
    function(model, net) {

      visreg(
        model,
        predictor,
        plot = FALSE,
        predict = list(
          re.form = NA
        )
      )
    }
  )


  # ----------------------------------------------------------
  # Partial residual points
  # ----------------------------------------------------------

  points <- imap_dfr(
    visreg_objects,
    function(v, net) {

      residual_col <- intersect(
        c(
          "visreg_res",
          "visregRes"
        ),
        names(v$res)
      )


      if (length(residual_col) == 0) {

        stop(
          paste(
            "Could not identify visreg residual column for",
            net
          )
        )
      }


      tibble(
        x = v$res[[predictor]],

        partial_residual =
          v$res[[residual_col[1]]],

        network = net
      )
    }
  )


  # ----------------------------------------------------------
  # Fitted fixed-effect lines
  # ----------------------------------------------------------

  lines <- imap_dfr(
    visreg_objects,
    function(v, net) {

      fit_col <- intersect(
        c(
          "visreg_fit",
          "visregFit"
        ),
        names(v$fit)
      )


      if (length(fit_col) == 0) {

        stop(
          paste(
            "Could not identify visreg fitted column for",
            net
          )
        )
      }


      tibble(
        x = v$fit[[predictor]],
        fitted = v$fit[[fit_col[1]]],
        network = net
      )
    }
  )


  list(

    points = points %>%
      mutate(
        network = factor(
          network,
          levels = target_cols
        )
      ),

    lines = lines %>%
      mutate(
        network = factor(
          network,
          levels = target_cols
        )
      )
  )
}


# ============================================================
# 11. Partial residual plot
# ============================================================

plot_fd_partial_residuals <- function(
  models,
  model_results,
  predictor = "mean_FD",
  title,
  y_limits = fixed_y_limits,
  y_breaks = fixed_y_breaks
) {

  vr <- get_partial_residual_data(
    models,
    predictor
  )


  labels <- model_results %>%
    mutate(
      x = Inf,
      y = Inf
    )


  ggplot(
    vr$points,
    aes(
      x = x,
      y = partial_residual,
      color = network
    )
  ) +

    # Partial residual observations
    geom_point(
      alpha = 0.3,
      size = 1
    ) +


    # Adjusted fitted line
    geom_line(
      data = vr$lines,
      aes(
        x = x,
        y = fitted,
        color = network
      ),
      inherit.aes = FALSE,
      linewidth = 1.2
    ) +


    # Model statistics
    geom_label(
      data = labels,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 5,
      color = "black",
      fill = "white",
      alpha = 0.8,
      linewidth = 0
    ) +


    facet_wrap(
      ~ network,
      ncol = 3,
      scales = "fixed",
      axes = "all",
      axis.labels = "all"
    ) +


    scale_color_manual(
      values = network_colors
    ) +


    scale_y_continuous(
      breaks = y_breaks
    ) +


    coord_cartesian(
      ylim = y_limits
    ) +


    network_plot_theme() +


    labs(
      title = title,
      x = "Mean framewise displacement",
      y = "Within-network connectivity (partial residual)"
    )
}


# ============================================================
# 12. Save figure as PDF + editable SVG
# ============================================================

save_figure <- function(
  plot,
  filename,
  width = 13,
  height = 7
) {

  # PDF — Cairo handles Unicode characters such as β
  ggsave(
    filename = file.path(
      out_dir,
      paste0(filename, ".pdf")
    ),
    plot = plot,
    device = cairo_pdf,
    width = width,
    height = height
  )

  # SVG for downstream vector editing
  ggsave(
    filename = file.path(
      out_dir,
      paste0(filename, ".svg")
    ),
    plot = plot,
    device = svglite::svglite,
    width = width,
    height = height
  )
}


# ============================================================
# 13. Print model table
# ============================================================

print_model_table <- function(
  results,
  title
) {

  cat(
    "\n============================================================\n",
    title,
    "\n============================================================\n",
    sep = ""
  )


  results %>%

    select(
      network,
      beta_adjusted,
      t_val_adjusted,
      p_adjusted,
      q_adjusted,
      sig_adjusted
    ) %>%

    mutate(

      beta_adjusted = round(
        beta_adjusted,
        4
      ),

      t_val_adjusted = round(
        t_val_adjusted,
        3
      ),

      p_adjusted = signif(
        p_adjusted,
        3
      ),

      q_adjusted = signif(
        q_adjusted,
        3
      )
    ) %>%

    as_tibble() %>%

    print(
      n = Inf
    )
}


# ============================================================
# 14. Run one complete FD analysis
# ============================================================

run_fd_analysis <- function(
  data_long,
  analysis_label,
  file_suffix
) {

  cat(
    "\n\n############################################################\n",
    analysis_label,
    "\n############################################################\n",
    sep = ""
  )


  # ----------------------------------------------------------
  # Models
  # ----------------------------------------------------------

  fit <- fit_fd_models(
    data_long
  )

  models <- fit$models

  results <- fit$results


  # ----------------------------------------------------------
  # Predictions
  # ----------------------------------------------------------

  predictions <- get_adjusted_predictions(
    models
  )


  # ----------------------------------------------------------
  # Print model results
  # ----------------------------------------------------------

  print_model_table(
    results,
    paste0(
      analysis_label,
      " adjusted FD model results"
    )
  )


  # ----------------------------------------------------------
  # Save model results CSV
  # ----------------------------------------------------------

  write_csv(
    results,
    file.path(
      out_dir,
      paste0(
        "fd_effects_withinconn_",
        file_suffix,
        "_modelresults.csv"
      )
    )
  )


  # ----------------------------------------------------------
  # Fitted-value plot
  # ----------------------------------------------------------

  p_fitted <- plot_fd_effects(
    data_long = data_long,
    predictions = predictions,
    model_results = results,
    title = paste0(
      analysis_label,
      ": adjusted FD effect by network"
    )
  )


  # ----------------------------------------------------------
  # Partial-residual plot
  # ----------------------------------------------------------

  p_partial <- plot_fd_partial_residuals(
    models = models,
    model_results = results,
    predictor = "mean_FD",
    title = paste0(
      analysis_label,
      ": partial residual FD effect by network"
    )
  )


  # ----------------------------------------------------------
  # Print plots
  # ----------------------------------------------------------

  print(
    p_fitted
  )

  print(
    p_partial
  )


  # ----------------------------------------------------------
  # Save PDF + SVG
  # ----------------------------------------------------------

  save_figure(
    p_fitted,
    paste0(
      "fd_effects_withinconn_",
      file_suffix,
      "_fitted"
    )
  )


  save_figure(
    p_partial,
    paste0(
      "fd_effects_withinconn_",
      file_suffix,
      "_partial_residuals"
    )
  )


  # ----------------------------------------------------------
  # Return everything
  # ----------------------------------------------------------

  list(
    models = models,
    results = results,
    predictions = predictions,
    fitted_plot = p_fitted,
    partial_plot = p_partial
  )
}


# ============================================================
# 15. RUN MODEL 1:
#     ALL OBSERVATIONS — NO FD THRESHOLD EXCLUSION
# ============================================================

results_all <- run_fd_analysis(

  data_long = data_long_all,

  analysis_label =
    "All observations",

  file_suffix =
    "all"
)


# ============================================================
# 16. RUN MODEL 2:
#     FD-FILTERED OBSERVATIONS
#     Keep mean FD < 0.25
# ============================================================

results_fd_filtered <- run_fd_analysis(

  data_long = data_long_fd_filtered,

  analysis_label = paste0(
    "FD filtered (mean FD < ",
    fd_threshold,
    ")"
  ),

  file_suffix =
    "fd_filtered"
)


# ============================================================
# 17. OPTIONAL: direct comparison of model coefficients
# ============================================================

comparison_results <- results_all$results %>%

  select(
    network,
    beta_all = beta_adjusted,
    t_all = t_val_adjusted,
    p_all = p_adjusted,
    q_all = q_adjusted
  ) %>%

  left_join(

    results_fd_filtered$results %>%

      select(
        network,
        beta_fd_filtered = beta_adjusted,
        t_fd_filtered = t_val_adjusted,
        p_fd_filtered = p_adjusted,
        q_fd_filtered = q_adjusted
      ),

    by = "network"
  )


cat(
  "\n============================================================\n",
  "COMPARISON: ALL vs FD-FILTERED MODELS\n",
  "============================================================\n",
  sep = ""
)

print(
  comparison_results,
  n = Inf
)


# Save comparison table
write_csv(
  comparison_results,
  file.path(
    out_dir,
    "fd_effects_withinconn_all_vs_fd_filtered.csv"
  )
)

# ============================================================
# 18. SELECTED-NETWORK COMPARISON FIGURE:
#     Partial residuals before vs after FD filtering
#
#     Top row:    Default | Limbic | SomMot — ALL
#     Bottom row: Default | Limbic | SomMot — FD < 0.25
# ============================================================


# Networks to show
selected_networks <- c(
  "Default",
  "Limbic",
  "SomMot"
)


# ------------------------------------------------------------
# Extract partial residual data for selected networks
# ------------------------------------------------------------

get_selected_residuals <- function(
  models,
  model_results,
  condition_label
) {

  # Only use the three selected networks
  selected_models <- models[
    selected_networks
  ]


  # Get partial residual points + fitted lines
  vr <- get_partial_residual_data(
    selected_models,
    predictor = "mean_FD"
  )


  # Add condition label
  points <- vr$points %>%
    filter(
      network %in% selected_networks
    ) %>%
    mutate(
      network = factor(
        network,
        levels = selected_networks
      ),
      condition = condition_label
    )


  lines <- vr$lines %>%
    filter(
      network %in% selected_networks
    ) %>%
    mutate(
      network = factor(
        network,
        levels = selected_networks
      ),
      condition = condition_label
    )


  # Model statistics for labels
  labels <- model_results %>%
    filter(
      network %in% selected_networks
    ) %>%
    mutate(
      network = factor(
        network,
        levels = selected_networks
      ),
      condition = condition_label,
      x = Inf,
      y = Inf,

      # ASCII "beta" avoids PDF encoding problems
      plot_label = sprintf(
        "beta = %.3f %s\nq = %.3g",
        beta_adjusted,
        sig_adjusted,
        q_adjusted
      )
    )


  list(
    points = points,
    lines = lines,
    labels = labels
  )
}


# ------------------------------------------------------------
# ALL observations
# ------------------------------------------------------------

pres_all <- get_selected_residuals(
  models = results_all$models,
  model_results = results_all$results,
  condition_label = "Before FD exclusion"
)


# ------------------------------------------------------------
# AFTER FD filtering
# ------------------------------------------------------------

pres_filtered <- get_selected_residuals(
  models = results_fd_filtered$models,
  model_results = results_fd_filtered$results,
  condition_label = "After FD exclusion"
)


# ============================================================
# 19. Make 2 x 3 selected-network plot
# ============================================================

# ============================================================
# Selected-network residual figure
# Separate plots so BEFORE and AFTER can have different x axes
# ============================================================


# ------------------------------------------------------------
# BEFORE FD exclusion
# ------------------------------------------------------------

p_before <- ggplot(
  pres_all$points,
  aes(
    x = x,
    y = partial_residual,
    color = network
  )
) +

  geom_point(
    alpha = 0.35,
    size = 0.8
  ) +

  geom_line(
    data = pres_all$lines,
    aes(
      x = x,
      y = fitted,
      color = network
    ),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +

  geom_label(
    data = pres_all$labels,
    aes(
      x = x,
      y = y,
      label = plot_label
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.1,
    size = 3.2,
    color = "black",
    fill = "white",
    alpha = 0.85,
    linewidth = 0
  ) +

  facet_grid(
    "Before FD exclusion" ~ network
  ) +

  scale_color_manual(
    values = network_colors
  ) +

  scale_y_continuous(
    breaks = fixed_y_breaks
  ) +

  coord_cartesian(
    ylim = fixed_y_limits
  ) +

  labs(
    x = "Mean framewise displacement",
    y = "Within-network connectivity\n(partial residual)"
  ) +

  theme_minimal(base_size = 10) +

  theme(
    legend.position = "none",

    strip.text.x = element_text(
      size = 11,
      face = "bold"
    ),

    # "Before FD exclusion" on side
    strip.text.y = element_text(
      size = 10,
      face = "bold"
    ),

    axis.title.x = element_text(
      size = 11,
      face = "bold"
    ),

    axis.title.y = element_text(
      size = 11,
      face = "bold"
    ),

    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),

    panel.grid.minor = element_blank(),

    panel.spacing.x = unit(
      0.35,
      "cm"
    )
  )


# ------------------------------------------------------------
# AFTER FD exclusion
# ------------------------------------------------------------

p_after <- ggplot(
  pres_filtered$points,
  aes(
    x = x,
    y = partial_residual,
    color = network
  )
) +

  geom_point(
    alpha = 0.35,
    size = 0.8
  ) +

  geom_line(
    data = pres_filtered$lines,
    aes(
      x = x,
      y = fitted,
      color = network
    ),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +

  geom_label(
    data = pres_filtered$labels,
    aes(
      x = x,
      y = y,
      label = plot_label
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.1,
    size = 3.2,
    color = "black",
    fill = "white",
    alpha = 0.85,
    linewidth = 0
  ) +
  facet_grid(
    "After FD exclusion" ~ network
  ) +

  scale_color_manual(
    values = network_colors
  ) +

  scale_y_continuous(
    breaks = fixed_y_breaks
  ) +

  coord_cartesian(
    ylim = fixed_y_limits
  ) +

  labs(
    title = paste0("FD < ", fd_threshold),
    x = "Mean framewise displacement",
    y = "Within-network connectivity\n(partial residual)"
  ) +

  theme_minimal(
    base_size = 10
  ) +

  theme(
    legend.position = "none",

    plot.title = element_text(
      size = 11,
      face = "bold"
    ),

    strip.text = element_text(
      size = 11,
      face = "bold"
    ),

    axis.title.x = element_text(
      size = 11,
      face = "bold"
    ),

    axis.title.y = element_text(
      size = 11,
      face = "bold"
    ),

    axis.text.x = element_text(
      size = 9
    ),

    axis.text.y = element_text(
      size = 9
    ),

    panel.grid.minor = element_blank(),

    panel.spacing.x = unit(
      0.35,
      "cm"
    )
  )


# ============================================================
# Stack the two rows
# ============================================================

p_residuals_selected <-
  p_before /
  p_after


print(p_residuals_selected)


# ============================================================
# Save selected-network comparison
# ============================================================

ggsave(
  file.path(
    out_dir,
    "fd_partial_residuals_selected_networks_before_after.pdf"
  ),
  p_residuals_selected,
  device = cairo_pdf,
  width = 18,
  height = 9,
  units = "cm"
)

ggsave(
  file.path(
    out_dir,
    "fd_partial_residuals_selected_networks_before_after.png"
  ),
  p_residuals_selected,
  width = 18,
  height = 9,
  units = "cm"
)
