library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(emmeans)
library(ggpubr)
library(readxl)
library(lme4)
library(lmerTest)
library(visreg)

emm_options(pbkrtest.limit = 10000)
emm_options(lmerTest.limit = 50000)

# ============================================================
# 1. Load and prepare data
# ============================================================

demos_withinconn <- read.csv("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/sheets/v1.3/demos_conn_2807.csv")

demos_withinconn <- demos_withinconn %>%
  mutate(
    dcfdx = factor(
      dcfdx,
      levels = c("NCI", "MCI", "AD", "other"),
      labels = c("NCI", "MCI", "AD", "other")
    ))

# drop rows that have other as dfcdx
demos_withinconn <- demos_withinconn %>%
  filter(dcfdx != "other") %>%
  droplevels()

# Network columns
target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

fd_threshold <- 0.25

# Type conversion

network_colors <- c(
  "Vis" = "#9B59B6",
  "SomMot" = "#6C8EBF",
  "Default" = "#D36B78",
  "Limbic" = "#C9D39A",
  "DorsAttn" = "#3C8D2F",
  "SalVentAttn" = "#C84CCF",
  "Cont" = "#E5B53A"
)

covariates_to_run <- c(
  "msex",
  "site",
  "eyes",
  "syn_bin",
  "dcfdx"
)

# Longitudinal random-intercept model
model_formula <- within_conn ~ mean_FD + msex + site + age_scandate +
  eyes + dcfdx + syn_bin + (1 | sub_id)

# ============================================================
# 2. Helper functions
# ============================================================

make_long <- function(data) {
  data %>%
    pivot_longer(
      cols = all_of(target_cols),
      names_to = "network",
      values_to = "within_conn"
    ) %>%
    mutate(
      network = factor(network, levels = target_cols),
      site = factor(
        site,
        levels = c("BNK", "UC", "MG", "RIRC")
      )
    ) %>%
    filter(
      !is.na(sub_id),
      !is.na(mean_FD),
      !is.na(within_conn),
      !is.na(msex),
      !is.na(site),
      !is.na(age_scandate),
      !is.na(eyes),
      !is.na(dcfdx),
      !is.na(syn_bin)
    ) %>%
    droplevels()
}

sig_from_p <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

# ============================================================
# 3. Fit network models ONCE
# ============================================================

fit_network_models <- function(data_long) {

  models <- map(target_cols, function(net) {

    df_net <- data_long %>%
      filter(network == net) %>%
      droplevels()

    lmer(
      model_formula,
      data = df_net,
      control = lmerControl(optimizer = "bobyqa")
    )
  })

  names(models) <- target_cols

  models
}


# ============================================================
# 4. Add predicted values
# ============================================================

get_predicted_data <- function(data_long, models) {

  map_dfr(target_cols, function(net) {

    df_net <- data_long %>%
      filter(network == net) %>%
      droplevels()

    df_net %>%
      mutate(
        predicted_conn = fitted(models[[net]])
      )

  }) %>%
    mutate(
      network = factor(network, levels = target_cols)
    )
}


# ============================================================
# 5. Pairwise comparisons with emmeans
# ============================================================

fit_model_pairwise <- function(models, covariate) {

  imap_dfr(models, function(model, net) {

    emm <- emmeans(
      model,
      specs = as.formula(paste("pairwise ~", covariate)),
      adjust = "tukey",
      lmer.df = "satterthwaite"
    )

    pairwise_df <- as.data.frame(emm$contrasts)

    stat_col <- intersect(
      c("t.ratio", "z.ratio"),
      names(pairwise_df)
    )[1]

    pairwise_df %>%
      transmute(
        network = net,
        covariate = covariate,
        contrast = contrast,
        estimate = estimate,
        se = SE,
        df = if ("df" %in% names(pairwise_df)) df else NA_real_,
        statistic = .data[[stat_col]],
        p_adj = p.value
      )

  }) %>%
    mutate(
      sig = sig_from_p(p_adj),
      network = factor(network, levels = target_cols)
    )
}


# ============================================================
# 6. General bracket annotations
#    Works for predicted OR partial residual data
# ============================================================

make_pairwise_brackets <- function(
  plot_data,
  pairwise_results,
  covariate,
  y_var
) {

  x_levels <- levels(
    factor(plot_data[[covariate]])
  )

  y_positions <- plot_data %>%
    group_by(network) %>%
    summarise(
      y_max = max(.data[[y_var]], na.rm = TRUE),
      y_min = min(.data[[y_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = y_max - y_min,
      y_range = ifelse(
        y_range == 0,
        1,
        y_range
      )
    )

  pairwise_results %>%
    filter(p_adj < 0.05) %>%

    separate(
      contrast,
      into = c("group1", "group2"),
      sep = " - ",
      remove = FALSE
    ) %>%

    mutate(
      group1_clean = str_remove(
        group1,
        paste0("^", covariate)
      ),

      group2_clean = str_remove(
        group2,
        paste0("^", covariate)
      ),

      group1 = ifelse(
        group1 %in% x_levels,
        group1,
        group1_clean
      ),

      group2 = ifelse(
        group2 %in% x_levels,
        group2,
        group2_clean
      )
    ) %>%

    select(
      -group1_clean,
      -group2_clean
    ) %>%

    left_join(
      y_positions,
      by = "network"
    ) %>%

    group_by(network) %>%
    arrange(p_adj, .by_group = TRUE) %>%

    mutate(
      bracket_number = row_number(),

      bracket_side = ifelse(
        bracket_number %% 2 == 1,
        "above",
        "below"
      ),

      bracket_rank = ceiling(
        bracket_number / 2
      ),

      y.position = ifelse(
        bracket_side == "above",

        y_max +
          bracket_rank *
          0.08 *
          y_range,

        y_min -
          bracket_rank *
          0.08 *
          y_range
      ),

      label = sig
    ) %>%

    ungroup()
}


# ============================================================
# 7. Predicted distributions
# ============================================================

plot_factor_covariate <- function(
  data_long,
  models,
  covariate,
  pairwise_results,
  title,
  y_limits = NULL,
  y_breaks = NULL
) {

  predicted_data <- get_predicted_data(
    data_long,
    models
  )

  pairwise_annot <- make_pairwise_brackets(
    plot_data = predicted_data,
    pairwise_results = pairwise_results,
    covariate = covariate,
    y_var = "predicted_conn"
  )

  p <- ggplot(
    predicted_data,
    aes(
      x = .data[[covariate]],
      y = predicted_conn,
      color = network,
      fill = network
    )
  ) +

    geom_violin(
      alpha = 0.18,
      linewidth = 0.3,
      trim = FALSE
    ) +

    geom_boxplot(
      width = 0.30,
      alpha = 0.75,
      outlier.shape = NA,
      linewidth = 0.35
    ) +

    geom_jitter(
      width = 0.12,
      alpha = 0.09,
      size = 0.4
    ) +

    facet_wrap(
      ~ network,
      ncol = 3,
      scales = "fixed"
    ) +

    scale_x_discrete(
      drop = FALSE
    ) +

    scale_color_manual(
      values = network_colors
    ) +

    scale_fill_manual(
      values = network_colors
    ) +

    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(
        mult = c(0.18, 0.18)
      )
    ) +

    coord_cartesian(
      ylim = y_limits
    ) +

    theme_minimal(
      base_size = 13
    ) +

    theme(
      legend.position = "none",
      strip.text = element_text(
        face = "bold"
      ),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    ) +

    labs(
      title = title,
      subtitle = paste0(
        "Violin/box/jitter show model-predicted values; ",
        "stars show Tukey-adjusted emmeans comparisons"
      ),
      x = covariate,
      y = "Predicted within-network connectivity"
    )

  if (nrow(pairwise_annot) > 0) {

    p <- p +
      ggpubr::stat_pvalue_manual(
        pairwise_annot,
        label = "label",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0,
        size = 4,
        bracket.size = 0.35,
        hide.ns = TRUE
      )
  }

  p
}


# ============================================================
# 8. Extract partial residuals for categorical covariate
# ============================================================

get_partial_residual_data <- function(
  models,
  covariate
) {

  imap_dfr(
    models,
    function(model, net) {

      v <- visreg(
        model,
        covariate,
        plot = FALSE,
        predict = list(
          re.form = NA
        )
      )

      # visreg >= 3.0
      if ("visreg_res" %in% names(v$res)) {
        residual_col <- "visreg_res"

      # fallback for older visreg
      } else if ("visregRes" %in% names(v$res)) {
        residual_col <- "visregRes"

      } else {
        stop(
          "Could not find the partial residual column in visreg output."
        )
      }

      out <- tibble(
        network = net,
        partial_residual = v$res[[residual_col]]
      )

      out[[covariate]] <- v$res[[covariate]]

      out %>%
        select(
          network,
          all_of(covariate),
          partial_residual
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
# 9. Plot partial residual distributions
# ============================================================

plot_partial_factor_covariate <- function(
  models,
  covariate,
  pairwise_results,
  title,
  y_limits = NULL,
  y_breaks = NULL
) {

  partial_data <- get_partial_residual_data(
    models = models,
    covariate = covariate
  )

  pairwise_annot <- make_pairwise_brackets(
    plot_data = partial_data,
    pairwise_results = pairwise_results,
    covariate = covariate,
    y_var = "partial_residual"
  )

  p <- ggplot(
    partial_data,
    aes(
      x = .data[[covariate]],
      y = partial_residual,
      color = network,
      fill = network
    )
  ) +

    # Same aesthetics as predicted plot
    geom_violin(
      alpha = 0.18,
      linewidth = 0.3,
      trim = FALSE
    ) +

    geom_boxplot(
      width = 0.30,
      alpha = 0.75,
      outlier.shape = NA,
      linewidth = 0.35
    ) +

    geom_jitter(
      width = 0.12,
      alpha = 0.09,
      size = 0.4
    ) +

    facet_wrap(
      ~ network,
      ncol = 3,
      scales = "fixed"
    ) +

    scale_x_discrete(
      drop = FALSE
    ) +

    scale_color_manual(
      values = network_colors
    ) +

    scale_fill_manual(
      values = network_colors
    ) +

    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(
        mult = c(0.18, 0.18)
      )
    ) +

    coord_cartesian(
      ylim = y_limits
    ) +

    theme_minimal(
      base_size = 13
    ) +

    theme(
      legend.position = "none",

      strip.text = element_text(
        face = "bold"
      ),

      panel.grid.minor = element_blank(),

      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    ) +

    labs(
      title = title,

      subtitle = paste0(
        "Violin/box/jitter show partial residuals; ",
        "stars show Tukey-adjusted emmeans comparisons"
      ),

      x = covariate,
      y = "Partial residual"
    )


  # Add same significance brackets
  if (nrow(pairwise_annot) > 0) {

    p <- p +
      ggpubr::stat_pvalue_manual(
        pairwise_annot,
        label = "label",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0,
        size = 4,
        bracket.size = 0.35,
        hide.ns = TRUE
      )
  }

  p
}

fixed_y_limits <- c(0, 0.6)
fixed_y_breaks <- seq(0, 0.6, by = 0.25)

# ============================================================
# 10. Print pairwise results
# ============================================================

print_pairwise_table <- function(
  pairwise_results,
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

  pairwise_results %>%
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
      p_adj = signif(
        p_adj,
        3
      )
    ) %>%
    as_tibble() %>%
    print(n = Inf)
}


# ============================================================
# 11. Wrapper
# ============================================================

run_factor_analysis <- function(
  data_long,
  models,
  covariate,
  analysis_label,
  file_suffix,
  y_limits = NULL,
  y_breaks = NULL
) {

  # ----------------------------------------------------------
  # Statistics
  # ----------------------------------------------------------

  pairwise_results <- fit_model_pairwise(
    models = models,
    covariate = covariate
  )

  print_pairwise_table(
    pairwise_results,
    paste0(
      analysis_label,
      ": model-based pairwise comparisons for ",
      covariate
    )
  )


  # ----------------------------------------------------------
  # Predicted-value plot
  # ----------------------------------------------------------

  p_predicted <- plot_factor_covariate(
    data_long = data_long,
    models = models,
    covariate = covariate,
    pairwise_results = pairwise_results,
    title = paste0(
      analysis_label,
      ": model-predicted distribution by ",
      covariate
    ),
    y_limits = y_limits,
    y_breaks = y_breaks
  )

  print(
    p_predicted
  )

  ggsave(
    filename = paste0(
      "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/",
      "withinconn_predicted_",
      covariate,
      "_",
      file_suffix,
      ".pdf"
    ),
    plot = p_predicted,
    width = 15,
    height = 9,
    dpi = 300
  )


  # ----------------------------------------------------------
  # Partial-residual plot
  # ----------------------------------------------------------

  p_partial <- plot_partial_factor_covariate(
    models = models,
    covariate = covariate,
    pairwise_results = pairwise_results,
    title = paste0(
      analysis_label,
      ": partial residual distribution by ",
      covariate
    ),
    y_limits = y_limits,
    y_breaks = y_breaks
  )

  print(
    p_partial
  )

  ggsave(
    filename = paste0(
      "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/",
      "withinconn_partial_residuals_",
      covariate,
      "_",
      file_suffix,
      ".pdf"
    ),
    plot = p_partial,
    width = 15,
    height = 9,
    dpi = 300
  )


  # ----------------------------------------------------------
  # Return everything
  # ----------------------------------------------------------

  list(
    pairwise = pairwise_results,
    predicted_plot = p_predicted,
    partial_plot = p_partial
  )
}

# ============================================================
# 12. Pre-filtering analysis
# ============================================================

data_long_pre <- make_long(
  demos_withinconn
)

# Fit 7 network models ONCE
models_pre <- fit_network_models(
  data_long_pre
)

pre_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    data_long = data_long_pre,
    models = models_pre,
    covariate = .x,
    analysis_label = "Pre-filtering",
    file_suffix = "preFDfilter",
    y_limits = fixed_y_limits,
    y_breaks = fixed_y_breaks
  )
)

names(pre_results) <- covariates_to_run

# ============================================================
# 13. Post-filtering analysis: mean_FD < 0.25
# ============================================================

data_long_post <- data_long_pre %>%
  filter(
    mean_FD < fd_threshold
  ) %>%
  droplevels()

# Fit 7 post-filter network models ONCE
models_post <- fit_network_models(
  data_long_post
)

post_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    data_long = data_long_post,
    models = models_post,
    covariate = .x,
    analysis_label = paste0(
      "Post-filtering, mean_FD < ",
      fd_threshold
    ),
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

all_pairwise_results_table <- all_pairwise_results %>%
  mutate(
    estimate = round(estimate, 4),
    se = round(se, 4),
    statistic = round(as.numeric(statistic), 3),
    p_adj = signif(p_adj, 3)
  ) %>%
  as_tibble()

print(all_pairwise_results_table, n = Inf)

write.csv(all_pairwise_results_table, "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/all_pairwise_results_withinconn_covs_long.csv", row.names = FALSE)