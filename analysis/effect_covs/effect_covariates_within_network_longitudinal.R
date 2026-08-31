# Longitudinal covariate-effects analysis for within-network connectivity.
# Fits participant-level mixed models and produces model-predicted covariate
# plots and pairwise comparison tables before/after FD filtering.

# Shared input/output path configuration. See README.md for environment variables.
.paths_file <- if (file.exists(file.path("analysis", "april26", "paths.R"))) {
  file.path("analysis", "april26", "paths.R")
} else {
  "paths.R"
}
source(.paths_file)
rm(.paths_file)

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
library(svglite)

emm_options(pbkrtest.limit = 10000)
emm_options(lmerTest.limit = 50000)

# ============================================================
# 1. Load and prepare data
# ============================================================

demos_withinconn <- read.csv(require_file(demos("demos_conn.csv")))

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
# Output settings
# ============================================================

out_dir <- output_dir("covariates_within_longitudinal")

fixed_y_limits <- c(0, 0.6)
fixed_y_breaks <- seq(0, 0.6, by = 0.25)


# ============================================================
# Fit network models ONCE per dataset
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

  set_names(models, target_cols)
}


# ============================================================
# Get model-predicted values
# ============================================================

get_predicted_data <- function(data_long, models) {

  map_dfr(target_cols, function(net) {

    data_long %>%
      filter(network == net) %>%
      droplevels() %>%
      mutate(
        predicted_conn = fitted(models[[net]])
      )
  }) %>%
    mutate(
      network = factor(network, levels = target_cols)
    )
}


# ============================================================
# Pairwise comparisons
# ============================================================

fit_model_pairwise <- function(models, covariate) {

  imap_dfr(models, function(model, net) {

    emm <- emmeans(
      model,
      specs = as.formula(paste("pairwise ~", covariate)),
      adjust = "tukey",
      lmer.df = "satterthwaite"
    )

    out <- as.data.frame(emm$contrasts)

    stat_col <- intersect(
      c("t.ratio", "z.ratio"),
      names(out)
    )[1]

    out %>%
      transmute(
        network = net,
        covariate = covariate,
        contrast,
        estimate,
        se = SE,
        df = if ("df" %in% names(out)) df else NA_real_,
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
# Extract partial residuals
# ============================================================

get_partial_residual_data <- function(models, covariate) {

  imap_dfr(models, function(model, net) {

    v <- visreg(
      model,
      covariate,
      plot = FALSE,
      predict = list(re.form = NA)
    )

    # visreg >= 3.0 uses visreg_res
    residual_col <- intersect(
      c("visreg_res", "visregRes"),
      names(v$res)
    )[1]

    if (is.na(residual_col)) {
      stop(
        paste(
          "Could not find partial residuals for",
          covariate,
          "in",
          net
        )
      )
    }

    tibble(
      network = net,
      group = v$res[[covariate]],
      value = v$res[[residual_col]]
    )
  }) %>%
    mutate(
      network = factor(network, levels = target_cols)
    )
}


# ============================================================
# Significance brackets
# ============================================================

make_pairwise_brackets <- function(
  plot_data,
  pairwise_results,
  covariate
) {

  x_levels <- if (is.factor(plot_data$group)) {
    levels(plot_data$group)
  } else {
    unique(as.character(plot_data$group))
  }

  y_positions <- plot_data %>%
    group_by(network) %>%
    summarise(
      y_min = min(value, na.rm = TRUE),
      y_max = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = y_max - y_min,
      y_range = ifelse(y_range == 0, 1, y_range)
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
      # Handles contrast labels with or without variable prefix
      group1 = ifelse(
        group1 %in% x_levels,
        group1,
        str_remove(group1, paste0("^", covariate))
      ),
      group2 = ifelse(
        group2 %in% x_levels,
        group2,
        str_remove(group2, paste0("^", covariate))
      )
    ) %>%
    left_join(
      y_positions,
      by = "network"
    ) %>%
    group_by(network) %>%
    arrange(p_adj, .by_group = TRUE) %>%
    mutate(
      bracket_number = row_number(),
      bracket_rank = ceiling(bracket_number / 2),

      y.position = ifelse(
        bracket_number %% 2 == 1,
        y_max + bracket_rank * 0.08 * y_range,
        y_min - bracket_rank * 0.08 * y_range
      ),

      label = sig
    ) %>%
    ungroup()
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

  pairwise_annot <- make_pairwise_brackets(
    plot_data = plot_data,
    pairwise_results = pairwise_results,
    covariate = covariate
  )

  p <- ggplot(
    plot_data,
    aes(
      x = group,
      y = value,
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

    # Axis ticks + labels on EVERY network panel
    facet_wrap(
      ~ network,
      ncol = 3,
      scales = "fixed",
      axes = "all",
      axis.labels = "all"
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
      breaks = if (is.null(y_breaks)) waiver() else y_breaks,
      expand = expansion(mult = c(0.18, 0.18))
    ) +

    coord_cartesian(
      ylim = y_limits
    ) +

    theme_minimal(base_size = 13) +

    theme(
      legend.position = "none",

      strip.text = element_text(
        face = "bold"
      ),

      panel.grid.minor = element_blank(),

      axis.text.x = element_text(
        size = 20
      ),

      axis.text.y = element_text(
        size = 20
      )
    ) +

    labs(
      title = title,
      subtitle = subtitle,
      x = covariate,
      y = y_lab
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
# Save PDF + editable SVG
# ============================================================

save_plot <- function(
  plot,
  filename,
  width = 13,
  height = 7
) {

  ggsave(
    file.path(
      out_dir,
      paste0(filename, ".pdf")
    ),
    plot = plot,
    width = width,
    height = height
  )

  ggsave(
    file.path(
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
# Run one covariate
# ============================================================

run_factor_analysis <- function(
  predicted_data,
  models,
  covariate,
  analysis_label,
  file_suffix,
  print_default_partial = FALSE
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
      "Violin/box/jitter show model-predicted values; stars show Tukey-adjusted emmeans comparisons",

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
      "Violin/box/jitter show partial residuals; stars show Tukey-adjusted emmeans comparisons",

    y_lab =
      "Within-network connectivity (partial residual)",

    # DON'T force partial residuals to 0–0.6
    y_limits = NULL,
    y_breaks = NULL
  )

  # ----------------------------------------------------------
  # Print
  # ----------------------------------------------------------

  if (print_default_partial) {

    default_partial_data <- partial_plot_data %>%
      filter(network == "Default") %>%
      droplevels()

    default_pairwise <- pairwise_results %>%
      filter(network == "Default")

    p_default_partial <- plot_factor_distribution(
      plot_data = default_partial_data,
      covariate = covariate,
      pairwise_results = default_pairwise,

      title = paste0(
        "Partial residual distribution by ",
        covariate
      ),

      subtitle =
        "Violin/box/jitter show partial residuals; stars show Tukey-adjusted emmeans comparisons",

      y_lab = "Within-network connectivity (partial residual)",

      y_limits = NULL,
      y_breaks = NULL
    )

    print(p_default_partial)
  }


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
    file_suffix = "preFDfilter",
    print_default_partial = TRUE
  )
) %>%
  set_names(covariates_to_run)

# ============================================================
# POST-FILTER
# ============================================================

data_long_post <- data_long_pre %>%
  filter(mean_FD < fd_threshold) %>%
  droplevels()

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