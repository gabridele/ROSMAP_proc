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
# 3. Fit the seven network models
# ============================================================

fit_fd_models <- function(data_long) {

  # Fit one model per network
  models <- map(
    target_cols,
    function(net) {

      df_net <- data_long %>%
        filter(network == net)

      lmer(
        within_conn ~
          mean_FD +
          ses_num +
          msex +
          site +
          age_scandate +
          eyes +
          syn_bin +
          dcfdx +
          (1 | sub_id),
        data = df_net
      )
    }
  )

  names(models) <- target_cols


  # Extract model results
  results <- imap_dfr(
    models,
    function(model, net) {

      coef_table <- summary(model)$coefficients

      tibble(
        network = net,
        beta_adjusted =
          coef_table["mean_FD", "Estimate"],

        t_val_adjusted =
          coef_table["mean_FD", "t value"],

        p_adjusted =
          coef_table["mean_FD", "Pr(>|t|)"]
      )
    }
  ) %>%

    # FDR correction across the 7 networks
    mutate(
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

  # Return BOTH models and results
  list(
    models = models,
    results = results
  )
}


# ============================================================
# 4. Get adjusted model predictions
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
        mutate(network = net)
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
# 5. Plot ordinary adjusted predictions
# ============================================================

plot_fd_effects <- function(
  data_long,
  pred_adjusted,
  model_results,
  title
) {

  label_pos <- model_results %>%
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
      data = pred_adjusted,
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
      data = pred_adjusted,
      aes(
        x = x,
        y = predicted,
        color = network
      ),
      inherit.aes = FALSE,
      linewidth = 1
    ) +

    geom_label(
      data = label_pos,
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
      scales = "fixed"
    ) +

    scale_y_continuous(
      breaks = seq(-0.2, 0.8, by = 0.2)
    ) +

    coord_cartesian(
      ylim = c(-0.2, 0.7)
    ) +

    scale_color_manual(
      values = network_colors
    ) +

    scale_fill_manual(
      values = network_colors
    ) +

    theme_minimal(
      base_size = 13
    ) +

    theme(
      legend.position = "none",
      strip.text = element_text(
        face = "bold"
      ),
      panel.grid.minor = element_blank()
    ) +

    labs(
      title = title,
      x = "Mean framewise displacement",
      y = "Within-network connectivity"
    )
}


# ============================================================
# 6. Partial residual plots
# ============================================================

plot_network_partial_residuals <- function(
  models,
  model_results,
  predictor = "mean_FD",
  network_colors,
  x_lab = "Mean framewise displacement",
  y_lab = "Partial residual",
  y_limits = c(-0.2, 0.7)
) {

  network_levels <- names(models)

  # ----------------------------------------------------------
  # Extract partial residuals from each model
  # ----------------------------------------------------------

  partial_points <- imap_dfr(
    models,
    function(model, net) {

      v <- visreg(
        model,
        predictor,
        plot = FALSE,
        predict = list(re.form = NA)
      )

      tibble(
        x = v$res[[predictor]],
        partial_residual = v$res$visreg_res,
        network = net
      )
    }
  ) %>%
    mutate(
      network = factor(
        network,
        levels = network_levels
      )
    )


  # ----------------------------------------------------------
  # Extract fitted fixed-effect line
  # ----------------------------------------------------------

  partial_lines <- imap_dfr(
    models,
    function(model, net) {

      v <- visreg(
        model,
        predictor,
        plot = FALSE,
        predict = list(re.form = NA)
      )

      tibble(
        x = v$fit[[predictor]],
        fitted = v$fit$visreg_fit,
        network = net
      )
    }
  ) %>%
    mutate(
      network = factor(
        network,
        levels = network_levels
      )
    )


  # ----------------------------------------------------------
  # Beta + q-value labels
  # ----------------------------------------------------------

  label_pos <- model_results %>%
    mutate(
      network = factor(
        network,
        levels = network_levels
      ),
      label = sprintf(
        "β = %.3f %s\nq = %.3g",
        beta_adjusted,
        sig_adjusted,
        q_adjusted
      ),
      x = Inf,
      y = Inf
    )


  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  ggplot(
    partial_points,
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
      data = partial_lines,
      aes(
        x = x,
        y = fitted,
        color = network
      ),
      inherit.aes = FALSE,
      linewidth = 1
    ) +

    geom_label(
      data = label_pos,
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
      scales = "fixed"
    ) +

    scale_y_continuous(
      breaks = seq(-0.2, 0.8, by = 0.2)
    ) +

    coord_cartesian(
      ylim = y_limits
    ) +

    scale_color_manual(
      values = network_colors
    ) +

    theme_minimal(
      base_size = 13
    ) +

    theme(
      legend.position = "none",
      strip.text = element_text(
        face = "bold"
      ),
      panel.grid.minor = element_blank()
    ) +

    labs(
      x = x_lab,
      y = y_lab
    )
}


# ============================================================
# 7. Print results table
# ============================================================

print_model_table <- function(
  model_results,
  title
) {

  cat(
    "\n============================================================\n"
  )

  cat(
    title,
    "\n"
  )

  cat(
    "============================================================\n"
  )

  model_results %>%
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

    print()
}


# ============================================================
# 8. PRE-FILTERING ANALYSIS
# ============================================================

data_long_pre <- make_long(
  demos_withinconn
)

pre <- fit_fd_models(
  data_long_pre
)

models_pre <- pre$models
model_results_pre <- pre$results

pred_adjusted_pre <- get_adjusted_predictions(
  models_pre
)


# Print results
print_model_table(
  model_results_pre,
  "Pre-filtering adjusted FD model results"
)


# Save results
write_csv(
  model_results_pre,
  "fd_effects_withinconn_longitudinal_preFDfilter_modelresults.csv"
)


# ------------------------------------------------------------
# Adjusted fitted-value plot
# ------------------------------------------------------------

p_pre_fitted <- plot_fd_effects(
  data_long_pre,
  pred_adjusted_pre,
  model_results_pre,
  "Adjusted FD effect by network: pre-filtering"
)

print(
  p_pre_fitted
)

ggsave(
  "fd_effects_withinconn_longitudinal_preFDfilter_fitted.png",
  plot = p_pre_fitted,
  width = 12,
  height = 8,
  dpi = 600
)


# ------------------------------------------------------------
# Partial residual plot
# ------------------------------------------------------------

p_pre_partial <- plot_network_partial_residuals(
  models = models_pre,
  model_results = model_results_pre,
  predictor = "mean_FD",
  network_colors = network_colors,
  x_lab = "Mean framewise displacement",
  y_lab = "Within-network connectivity (partial residual)",
  y_limits = c(-0.2, 0.7)
)

print(
  p_pre_partial
)

ggsave(
  "fd_effects_withinconn_longitudinal_preFDfilter_partial_residuals.png",
  plot = p_pre_partial,
  width = 12,
  height = 8,
  dpi = 600
)


# ============================================================
# 9. POST-FILTERING ANALYSIS
# ============================================================

data_long_post <- data_long_pre %>%
  filter(
    mean_FD < fd_threshold
  )

post <- fit_fd_models(
  data_long_post
)

models_post <- post$models
model_results_post <- post$results

pred_adjusted_post <- get_adjusted_predictions(
  models_post
)


# Print results
print_model_table(
  model_results_post,
  paste0(
    "Post-filtering adjusted FD model results: mean_FD < ",
    fd_threshold
  )
)


# Save results
write_csv(
  model_results_post,
  "fd_effects_withinconn_longitudinal_postFDfilter_modelresults.csv"
)


# ------------------------------------------------------------
# Adjusted fitted-value plot
# ------------------------------------------------------------

p_post_fitted <- plot_fd_effects(
  data_long_post,
  pred_adjusted_post,
  model_results_post,
  paste0(
    "Adjusted FD effect by network, mean_FD < ",
    fd_threshold
  )
)

print(
  p_post_fitted
)

ggsave(
  "fd_effects_withinconn_longitudinal_postFDfilter_fitted.png",
  plot = p_post_fitted,
  width = 12,
  height = 8,
  dpi = 600
)


# ------------------------------------------------------------
# Partial residual plot
# ------------------------------------------------------------

p_post_partial <- plot_network_partial_residuals(
  models = models_post,
  model_results = model_results_post,
  predictor = "mean_FD",
  network_colors = network_colors,
  x_lab = "Mean framewise displacement",
  y_lab = "Within-network connectivity (partial residual)",
  y_limits = c(-0.2, 0.7)
)

print(
  p_post_partial
)

ggsave(
  "fd_effects_withinconn_longitudinal_postFDfilter_partial_residuals.png",
  plot = p_post_partial,
  width = 12,
  height = 8,
  dpi = 600
)