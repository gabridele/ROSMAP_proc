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
source(require_file("analysis/utils.R"))
rm(.script_arg, .script_dir, .paths_candidates, .paths_file)


# ============================================================
# 1. Load and prepare data
# ============================================================

demos_betweenconn <- read.csv(require_file(demos("demos_conn.csv")))

# Add numeric session
demos_betweenconn <- demos_betweenconn %>%
  mutate(
    ses_num = as.numeric(str_extract(ses_id, "\\d+"))
  )

# The prepared analysis table already contains diagnosis and SyN variables.

# Missingness check
print(colSums(is.na(demos_betweenconn)))

fd_threshold <- FD_THRESHOLD

# Convert model variables explicitly so results do not depend on CSV type inference.
demos_betweenconn <- demos_betweenconn %>%
  mutate(
    sub_id = factor(sub_id),
    mean_FD = as.numeric(mean_FD),
    msex = factor(msex),
    site = factor(site),
    age_scandate = as.numeric(age_scandate),
    eyes = factor(eyes),
    dcfdx = factor(dcfdx),
    syn_bin = factor(syn_bin)
  )

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
  data %>%
    pivot_longer(
      cols = all_of(target_combos),
      names_to = "network_combo",
      values_to = "between_conn"
    ) %>%
    mutate(
      network_combo = factor(network_combo, levels = target_combos)
    ) %>%
    filter(
      !is.na(sub_id),
      !is.na(mean_FD),
      !is.na(between_conn),
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
# 3. Fit adjusted model and add predicted values
# ============================================================

get_predicted_data <- function(data_long) {
  
  map_dfr(target_combos, function(net) {
    
    df_net <- data_long %>%
      filter(network_combo == net) %>%
      droplevels()
    
    model <- lmer(
      model_formula,
      data = df_net,
      control = lmerControl(optimizer = "bobyqa")
    )
    
    df_net %>%
      mutate(
        predicted_conn = fitted(model)
      )
  }) %>%
    mutate(
      network_combo = factor(network_combo, levels = target_combos)
    )
}

# ============================================================
# 4. Model-based pairwise comparisons with emmeans
# ============================================================

fit_model_pairwise <- function(data_long, covariate) {
  
  map_dfr(target_combos, function(net) {
    
    df_net <- data_long %>%
      filter(network_combo == net) %>%
      droplevels()
    
    model <- lmer(
      model_formula,
      data = df_net,
      control = lmerControl(optimizer = "bobyqa")
    )
    
    emm <- emmeans(
      model,
      specs = as.formula(paste("pairwise ~", covariate)),
      adjust = "tukey",
      lmer.df = "satterthwaite"
    )
    
    pairwise_df <- as.data.frame(emm$contrasts)
    
    stat_col <- intersect(c("t.ratio", "z.ratio"), names(pairwise_df))[1]
    
    pairwise_df %>%
      transmute(
        network_combo = net,
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
      network_combo = factor(network_combo, levels = target_combos)
    )
}

# ============================================================
# Extract partial residuals
# ============================================================

get_partial_residual_data <- function(data_long, covariate) {

  map_dfr(target_combos, function(net) {

    df_net <- data_long %>%
      filter(network_combo == net) %>%
      droplevels()

    model <- lmer(
      model_formula,
      data = df_net,
      control = lmerControl(optimizer = "bobyqa")
    )

    v <- visreg(
      model,
      covariate,
      plot = FALSE,
      predict = list(re.form = NA)
    )

    # visreg naming differs between versions
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
      network_combo = net,
      group = v$res[[covariate]],
      partial_residual = v$res[[residual_col]]
    )
  }) %>%
    mutate(
      network_combo = factor(
        network_combo,
        levels = target_combos
      )
    )
}


# ============================================================
# Brackets for partial-residual plots
# ============================================================

make_partial_residual_brackets <- function(
  partial_data,
  pairwise_results,
  covariate
) {

  x_levels <- if (is.factor(partial_data$group)) {
    levels(partial_data$group)
  } else {
    unique(as.character(partial_data$group))
  }

  y_positions <- partial_data %>%
    group_by(network_combo) %>%
    summarise(
      y_min = min(partial_residual, na.rm = TRUE),
      y_max = max(partial_residual, na.rm = TRUE),
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
      by = "network_combo"
    ) %>%
    group_by(network_combo) %>%
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
# Plot partial residual distributions
# ============================================================

plot_partial_residual_covariate <- function(
  data_long,
  covariate,
  pairwise_results,
  title
) {

  partial_data <- get_partial_residual_data(
    data_long = data_long,
    covariate = covariate
  )

  pairwise_annot <- make_partial_residual_brackets(
    partial_data = partial_data,
    pairwise_results = pairwise_results,
    covariate = covariate
  )

  p <- ggplot(
    partial_data,
    aes(
      x = group,
      y = partial_residual,
      color = network_combo,
      fill = network_combo
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
      ~ network_combo,
      scales = "fixed"
    ) +

    scale_x_discrete(
      drop = FALSE
    ) +

    scale_color_manual(
      values = between_network_colors
    ) +

    scale_fill_manual(
      values = between_network_colors
    ) +

    scale_y_continuous(
      expand = expansion(mult = c(0.18, 0.18))
    ) +

    theme_minimal(base_size = 13) +

    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    ) +

    labs(
      title = title,
      subtitle =
        "Violin/box/jitter show partial residuals; stars show Tukey-adjusted emmeans comparisons",
      x = covariate,
      y = "Partial residual"
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
# 5. Bracket annotations
# ============================================================

make_pairwise_brackets <- function(predicted_data, pairwise_results, covariate) {
  
  x_levels <- levels(factor(predicted_data[[covariate]]))
  
  y_positions <- predicted_data %>%
    group_by(network_combo) %>%
    summarise(
      y_max = max(predicted_conn, na.rm = TRUE),
      y_min = min(predicted_conn, na.rm = TRUE),
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
      group1_clean = str_remove(group1, paste0("^", covariate)),
      group2_clean = str_remove(group2, paste0("^", covariate)),
      group1 = ifelse(group1 %in% x_levels, group1, group1_clean),
      group2 = ifelse(group2 %in% x_levels, group2, group2_clean)
    ) %>%
    select(-group1_clean, -group2_clean) %>%
    left_join(y_positions, by = "network_combo") %>%
    group_by(network_combo) %>%
    arrange(p_adj, .by_group = TRUE) %>%
    mutate(
      bracket_number = row_number(),
      bracket_side = ifelse(bracket_number %% 2 == 1, "above", "below"),
      bracket_rank = ceiling(bracket_number / 2),
      y.position = ifelse(
        bracket_side == "above",
        y_max + bracket_rank * 0.08 * y_range,
        y_min - bracket_rank * 0.08 * y_range
      ),
      label = sig
    ) %>%
    ungroup()
}

# ============================================================
# 6. Plot predicted distributions
# ============================================================

plot_factor_covariate <- function(data_long, covariate, pairwise_results, title,
                                  y_limits = NULL,
                                  y_breaks = NULL) {
  
  predicted_data <- get_predicted_data(data_long)
  
  pairwise_annot <- make_pairwise_brackets(
    predicted_data = predicted_data,
    pairwise_results = pairwise_results,
    covariate = covariate
  )
  
  p <- ggplot(
    predicted_data,
    aes(
      x = .data[[covariate]],
      y = predicted_conn,
      color = network_combo,
      fill = network_combo
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
    facet_wrap(~ network_combo, scales = "fixed") +
    scale_x_discrete(drop = FALSE) +
    scale_color_manual(values = between_network_colors) +
    scale_fill_manual(values = between_network_colors) +
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0.18, 0.18))
    ) +
    coord_cartesian(
      ylim = y_limits
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = title,
      subtitle = "Violin/box/jitter show model-predicted values; stars show Tukey-adjusted emmeans comparisons",
      x = covariate,
      y = "Predicted between-network connectivity"
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
  
  return(p)
}

# ============================================================
# 7. Print pairwise tables
# ============================================================

print_pairwise_table <- function(pairwise_results, title) {
  
  cat("\n============================================================\n")
  cat(title, "\n")
  cat("============================================================\n")
  
  pairwise_results %>%
    mutate(
      estimate = round(estimate, 4),
      se = round(se, 4),
      statistic = round(as.numeric(statistic), 3),
      p_adj = signif(p_adj, 3)
    ) %>%
    as_tibble() %>%
    print(n = Inf)
}

# ============================================================
# 8. Wrapper
# ============================================================

run_factor_analysis <- function(
  data_long,
  covariate,
  analysis_label,
  file_suffix,
  y_limits = NULL,
  y_breaks = NULL,
  make_partial_residuals = FALSE
) {

  pairwise_results <- fit_model_pairwise(
    data_long = data_long,
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

  p <- plot_factor_covariate(
    data_long = data_long,
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

  print(p)

  ggsave(
    filename = output(
      "covariates",
      paste0(
        "betweenconn_predicted_",
        covariate,
        "_",
        file_suffix,
        ".png"
      )
    ),
    plot = p,
    width = 13,
    height = 9,
    dpi = 300
  )


  # ----------------------------------------------------------
  # Partial residual plot
  # Only produced when explicitly requested
  # ----------------------------------------------------------

  p_partial <- NULL

  if (make_partial_residuals) {

    p_partial <- plot_partial_residual_covariate(
      data_long = data_long,
      covariate = covariate,
      pairwise_results = pairwise_results,

      # no need to say pre-filtering in the title
      title = paste0(
        "Partial residual distribution by ",
        covariate
      )
    )

    print(p_partial)

    ggsave(
      filename = output(
        "covariates",
        paste0(
          "betweenconn_partial_residuals_",
          covariate,
          "_",
          file_suffix,
          ".png"
        )
      ),
      plot = p_partial,
      width = 13,
      height = 9,
      dpi = 300
    )
  }


  list(
    pairwise = pairwise_results,
    plot = p,
    partial_plot = p_partial
  )
}

fixed_y_limits <- c(-0.2, 0.25)
fixed_y_breaks <- seq(-0.2, 0.25, by = 0.25)

# ============================================================
# 9. Pre-filtering analysis
# ============================================================

data_long_pre <- make_long(demos_betweenconn)

pre_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    data_long = data_long_pre,
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
  filter(mean_FD < fd_threshold) %>%
  droplevels()

post_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    data_long = data_long_post,
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

all_pairwise_results_table <- all_pairwise_results %>%
  mutate(
    estimate = round(estimate, 4),
    se = round(se, 4),
    statistic = round(as.numeric(statistic), 3),
    p_adj = signif(p_adj, 3)
  ) %>%
  as_tibble()

print(all_pairwise_results_table, n = Inf)
# save final table
write_csv(all_pairwise_results_table, output("covariates", "all_pairwise_results_betweenconn_covs_long.csv"))