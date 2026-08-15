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
  "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/sheets/v1.3/demos_conn_2807.csv"
)

demos_withinconn <- demos_withinconn %>%
  mutate(
    ses_num = as.numeric(str_extract(ses_id, "\\d+")),
    sub_id = factor(sub_id)
  )

target_cols <- c(
  "Vis",
  "SomMot",
  "DorsAttn",
  "SalVentAttn",
  "Limbic",
  "Cont",
  "Default"
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
      network = factor(network, levels = target_cols)
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
# 3. Global plotting / output settings
# ============================================================

out_dir <- "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26"

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

fixed_y_limits <- c(-0.2, 0.7)
fixed_y_breaks <- seq(-0.2, 0.8, by = 0.2)


# ============================================================
# 4. Fit network models + extract statistics
# ============================================================

fit_fd_models <- function(data_long) {

  models <- map(target_cols, function(net) {

    df_net <- data_long %>%
      filter(network == net)

    lmer(
      fd_formula,
      data = df_net,
      control = lmerControl(
        optimizer = "bobyqa"
      )
    )
  }) %>%
    set_names(target_cols)


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
      # BH/FDR correction across 7 networks
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
# 5. Fixed-effect predictions
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
# 6. Common figure theme
# ============================================================

network_plot_theme <- function() {

  theme_minimal(base_size = 13) +

    theme(
      legend.position = "none",

      strip.text = element_text(
        face = "bold"
      ),

      panel.grid.minor = element_blank(),

      axis.text.x = element_text(
        size = 9
      ),

      axis.text.y = element_text(
        size = 9
      )
    )
}


# ============================================================
# 7. Fitted-value plot
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

    geom_point(
      alpha = 0.3,
      size = 1
    ) +

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

    geom_line(
      data = predictions,
      aes(
        x = x,
        y = predicted,
        color = network
      ),
      inherit.aes = FALSE,
      linewidth = 1
    ) +

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
      size = 3,
      color = "black",
      fill = "white",
      alpha = 0.8,
      linewidth = 0
    ) +

    # Axes + tick labels on EVERY panel
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
      x = "Mean framewise displacement",
      y = "Within-network connectivity"
    )
}


# ============================================================
# 8. Get visreg partial residual data
# ============================================================

get_partial_residual_data <- function(
  models,
  predictor = "mean_FD"
) {

  # Run visreg only ONCE per model
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


  # -------------------------
  # Partial residual points
  # -------------------------

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


  # -------------------------
  # Fitted fixed-effect lines
  # -------------------------

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
# 9. Partial residual plot
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

    geom_point(
      alpha = 0.3,
      size = 1
    ) +

    geom_line(
      data = vr$lines,
      aes(
        x = x,
        y = fitted,
        color = network
      ),
      inherit.aes = FALSE,
      linewidth = 1
    ) +

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
      size = 3,
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
# 10. Save figure as PDF + editable SVG
# ============================================================

save_figure <- function(
  plot,
  filename,
  width = 13,
  height = 7
) {

  # PDF
  ggsave(
    filename = file.path(
      out_dir,
      paste0(filename, ".pdf")
    ),
    plot = plot,
    width = width,
    height = height
  )

  # SVG — editable in Illustrator
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
# 11. Print model table
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
    print(n = Inf)
}


# ============================================================
# 12. Run entire FD analysis
# ============================================================

run_fd_analysis <- function(
  data_long,
  analysis_label,
  file_suffix
) {

  # -------------------------
  # Models
  # -------------------------

  fit <- fit_fd_models(
    data_long
  )

  models <- fit$models
  results <- fit$results


  # -------------------------
  # Predictions
  # -------------------------

  predictions <- get_adjusted_predictions(
    models
  )


  # -------------------------
  # Results
  # -------------------------

  print_model_table(
    results,
    paste0(
      analysis_label,
      " adjusted FD model results"
    )
  )

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


  # -------------------------
  # Fitted-value plot
  # -------------------------

  p_fitted <- plot_fd_effects(
    data_long = data_long,
    predictions = predictions,
    model_results = results,
    title = paste0(
      analysis_label,
      ": adjusted FD effect by network"
    )
  )


  # -------------------------
  # Partial-residual plot
  # -------------------------

  p_partial <- plot_fd_partial_residuals(
    models = models,
    model_results = results,
    predictor = "mean_FD",
    title = paste0(
      analysis_label,
      ": partial residual FD effect by network"
    )
  )


  # Print
  print(p_fitted)
  print(p_partial)


  # Save PDF + SVG
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


  # Return everything
  list(
    models = models,
    results = results,
    predictions = predictions,
    fitted_plot = p_fitted,
    partial_plot = p_partial
  )
}
