library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(emmeans)
library(ggpubr)
library(readxl)

# ============================================================
# 1. Load and prepare data
# ============================================================

demos_withinconn <- read.csv("sheets/v1.3/demos_conn_2306.csv")

demos_withinconn <- demos_withinconn %>%
  mutate(
    dcfdx = factor(
      dcfdx,
      levels = c("NCI", "MCI", "AD", "other"),
      labels = c("NCI", "MCI", "AD", "other")
    ))

# Missingness check
print(colSums(is.na(demos_withinconn)))

# Network columns
target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

fd_threshold <- 0.25

# drop rows that have other as dfcdx
demos_withinconn <- demos_withinconn %>%
  filter(dcfdx != "other") %>%
  droplevels()

# Network colors
network_colors <- c(
  "Vis" = "#9B59B6",
  "SomMot" = "#6C8EBF",
  "Default" = "#D36B78",
  "Limbic" = "#C9D39A",
  "DorsAttn" = "#3C8D2F",
  "SalVentAttn" = "#C84CCF",
  "Cont" = "#E5B53A"
)

# Categorical covariates to test/plot
covariates_to_run <- c(
  "msex",
  "site",
  "eyes",
  "syn_bin",
  "dcfdx"
)

# Adjusted model formula
model_formula <- within_conn ~ mean_FD + msex + site + age_scandate +
  eyes + dcfdx + syn_bin

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
      network = factor(network, levels = target_cols)
    ) %>%
    filter(
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
# 3. Fit adjusted model and add predicted values
# ============================================================

get_predicted_data <- function(data_long) {
  
  map_dfr(target_cols, function(net) {
    
    df_net <- data_long %>%
      filter(network == net) %>%
      droplevels()
    
    model <- lm(
      model_formula,
      data = df_net
    )
    
    df_net %>%
      mutate(
        predicted_conn = fitted(model)
      )
  }) %>%
    mutate(
      network = factor(network, levels = target_cols)
    )
}

# ============================================================
# 4. Model-based pairwise comparisons with emmeans
# ============================================================

fit_model_pairwise <- function(data_long, covariate) {
  
  map_dfr(target_cols, function(net) {
    
    df_net <- data_long %>%
      filter(network == net) %>%
      droplevels()
    
    model <- lm(
      model_formula,
      data = df_net
    )
    
    emm <- emmeans(
      model,
      specs = as.formula(paste("pairwise ~", covariate)),
      adjust = "tukey"
    )
    
    pairwise_df <- as.data.frame(emm$contrasts)
    
    stat_col <- intersect(c("t.ratio", "z.ratio"), names(pairwise_df))[1]
    
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
# 5. Bracket annotations
# ============================================================

make_pairwise_brackets <- function(predicted_data, pairwise_results, covariate) {
  
  x_levels <- levels(factor(predicted_data[[covariate]]))
  
  y_positions <- predicted_data %>%
    group_by(network) %>%
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
    left_join(y_positions, by = "network") %>%
    group_by(network) %>%
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
    facet_wrap(~ network, scales = "fixed") +
    scale_x_discrete(drop = FALSE) +
    scale_color_manual(values = network_colors) +
    scale_fill_manual(values = network_colors) +
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

run_factor_analysis <- function(data_long, covariate, analysis_label, file_suffix,
                                y_limits = NULL,
                                y_breaks = NULL) {
  
  pairwise_results <- fit_model_pairwise(
    data_long = data_long,
    covariate = covariate
  )
  
  print_pairwise_table(
    pairwise_results,
    paste0(analysis_label, ": model-based pairwise comparisons for ", covariate)
  )
  
  p <- plot_factor_covariate(
    data_long = data_long,
    covariate = covariate,
    pairwise_results = pairwise_results,
    title = paste0(analysis_label, ": model-predicted distribution by ", covariate),
    y_limits = y_limits,
    y_breaks = y_breaks
  )
  
  print(p)
  
  ggsave(
    filename = paste0("withinconn_predicted_", covariate, "_", file_suffix, ".png"),
    plot = p,
    width = 13,
    height = 9,
    dpi = 300
  )
  
  list(
    pairwise = pairwise_results,
    plot = p
  )
}

fixed_y_limits <- c(0, 0.6)
fixed_y_breaks <- seq(0, 0.6, by = 0.25)

# ============================================================
# 9. Pre-filtering analysis
# ============================================================

data_long_pre <- make_long(demos_withinconn)

pre_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    data_long = data_long_pre,
    covariate = .x,
    analysis_label = "Pre-filtering",
    file_suffix = "preFDfilter",
    y_limits = fixed_y_limits,
    y_breaks = fixed_y_breaks
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
write_csv(all_pairwise_results_table, "all_pairwise_results_withinconn_covs_bl.csv")
