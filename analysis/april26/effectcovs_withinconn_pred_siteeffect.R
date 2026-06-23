library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(tibble)

# ============================================================
# 1. Load data
# ============================================================

demos_withinconn <- read.csv("sheets/v1.3/demos_conn_2306.csv")

print(colSums(is.na(demos_withinconn)))

# Drop small "other" diagnosis group
demos_withinconn <- demos_withinconn %>%
  filter(dcfdx != "other")

# ============================================================
# 2. Settings
# ============================================================

target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

covariates_to_run <- c(
  "msex",
  "eyes",
  "syn_bin",
  "dcfdx"
)

model_formula <- within_conn ~ mean_FD + msex + age_scandate +
  eyes + dcfdx + syn_bin + (1 | site)

fixed_y_limits <- c(0, 0.6)
fixed_y_breaks <- seq(0, 0.6, by = 0.25)

emm_options(
  lmer.df = "satterthwaite",
  lmerTest.limit = 10000,
  pbkrtest.limit = 10000
)

# ============================================================
# 3. Helper functions
# ============================================================

sig_from_p <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

safe_filename <- function(x) {
  str_replace_all(as.character(x), "[^A-Za-z0-9_-]", "_")
}

make_long <- function(data) {
  data %>%
    pivot_longer(
      cols = all_of(target_cols),
      names_to = "network",
      values_to = "within_conn"
    ) %>%
    mutate(
      network = factor(network, levels = target_cols),
      site = factor(site),
      msex = factor(msex),
      eyes = factor(eyes),
      syn_bin = factor(syn_bin),
      dcfdx = factor(dcfdx)
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

data_long_pre <- make_long(demos_withinconn)

# ============================================================
# 4. Fit LMER and get conditional predicted values
# ============================================================

get_conditional_predicted_data <- function(data_long) {
  
  map_dfr(target_cols, function(net) {
    
    df_net <- data_long %>%
      filter(network == net) %>%
      droplevels() %>%
      as.data.frame()
    
    model <- lmer(
      model_formula,
      data = df_net,
      REML = FALSE
    )
    
    df_net %>%
      mutate(
        predicted_conditional = predict(model, re.form = NULL)
      )
  }) %>%
    mutate(
      network = factor(network, levels = target_cols)
    )
}

# ============================================================
# 5. emmeans pairwise comparisons
# ============================================================

fit_emmeans_pairwise <- function(data_long, covariate) {
  
  map_dfr(target_cols, function(net) {
    
    df_net <- data_long %>%
      filter(network == net) %>%
      droplevels() %>%
      as.data.frame()
    
    # Create a temporary global data object so emmeans/lmerTest
    # can find the model data during Satterthwaite df calculation.
    temp_data_name <- paste0(
      ".tmp_emm_data_",
      safe_filename(net),
      "_",
      safe_filename(covariate)
    )
    
    assign(temp_data_name, df_net, envir = .GlobalEnv)
    
    on.exit(
      rm(list = temp_data_name, envir = .GlobalEnv),
      add = TRUE
    )
    
    model <- eval(
      substitute(
        lmerTest::lmer(
          model_formula,
          data = TEMP_DATA,
          REML = FALSE
        ),
        list(TEMP_DATA = as.name(temp_data_name))
      ),
      envir = .GlobalEnv
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
# 6. Bracket annotations
# ============================================================

make_pairwise_brackets <- function(predicted_data, pairwise_results, covariate) {
  
  x_levels <- levels(factor(predicted_data[[covariate]]))
  
  y_positions <- predicted_data %>%
    group_by(network) %>%
    summarise(
      y_max = max(predicted_conditional, na.rm = TRUE),
      y_min = min(predicted_conditional, na.rm = TRUE),
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
      bracket_rank = bracket_number,
      y.position = y_max + bracket_rank * 0.08 * y_range,
      label = sig
    ) %>%
    ungroup()
}

# ============================================================
# 7. Print pairwise table
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
# 8. Plot conditional predictions with site-colored jitter + emmeans stars
# ============================================================

plot_conditional_by_covariate <- function(data_long, covariate, pairwise_results,
                                          file_suffix = "preFDfilter",
                                          y_limits = NULL,
                                          y_breaks = NULL) {
  
  predicted_data <- get_conditional_predicted_data(data_long)
  
  pairwise_annot <- make_pairwise_brackets(
    predicted_data = predicted_data,
    pairwise_results = pairwise_results,
    covariate = covariate
  )
  
  p <- ggplot(
    predicted_data,
    aes(
      x = .data[[covariate]],
      y = predicted_conditional
    )
  ) +
    geom_violin(
      aes(group = .data[[covariate]]),
      fill = "grey80",
      color = "grey40",
      alpha = 0.25,
      linewidth = 0.3,
      trim = FALSE
    ) +
    geom_boxplot(
      aes(group = .data[[covariate]]),
      width = 0.25,
      alpha = 0.45,
      outlier.shape = NA,
      linewidth = 0.35,
      color = "black"
    ) +
    geom_jitter(
      aes(color = site),
      width = 0.11,
      alpha = 0.2,
      size = 0.4
    ) +
    facet_wrap(~ network, scales = "fixed") +
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0.18, 0.18))
    ) +
    coord_cartesian(ylim = y_limits) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste0("LMER-predicted values by ", covariate),
      subtitle = "Model includes site as random intercept: (1 | site). Scatter plot is color coded for site. Stars are Tukey-adjusted emmeans comparisons.",
      x = covariate,
      y = "Predicted within-network connectivity",
      color = "Site"
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
  
  print(p)
  
  ggsave(
    filename = paste0(
      "withinconn_lmer_conditional_predicted_",
      covariate,
      "_colored_by_site_emmeans_",
      file_suffix,
      ".png"
    ),
    plot = p,
    width = 14,
    height = 9,
    dpi = 300
  )
  
  return(p)
}

# ============================================================
# 9. Run one covariate
# ============================================================

run_one_covariate <- function(data_long, covariate, file_suffix = "preFDfilter") {
  
  pairwise_results <- fit_emmeans_pairwise(
    data_long = data_long,
    covariate = covariate
  )
  
  print_pairwise_table(
    pairwise_results,
    paste0("LMER emmeans pairwise comparisons for ", covariate)
  )
  
  p <- plot_conditional_by_covariate(
    data_long = data_long,
    covariate = covariate,
    pairwise_results = pairwise_results,
    file_suffix = file_suffix,
    y_limits = fixed_y_limits,
    y_breaks = fixed_y_breaks
  )
  
  list(
    pairwise = pairwise_results,
    plot = p
  )
}

# ============================================================
# 10. Run all covariates
# ============================================================

results <- map(
  covariates_to_run,
  ~ run_one_covariate(
    data_long = data_long_pre,
    covariate = .x,
    file_suffix = "preFDfilter"
  )
)

names(results) <- covariates_to_run

# ============================================================
# 11. Save combined pairwise table
# ============================================================

all_pairwise_results <- map_dfr(names(results), function(cov) {
  results[[cov]]$pairwise
})

all_pairwise_results_table <- all_pairwise_results %>%
  mutate(
    estimate = round(estimate, 4),
    se = round(se, 4),
    statistic = round(as.numeric(statistic), 3),
    p_adj = signif(p_adj, 3)
  ) %>%
  as_tibble()

print(all_pairwise_results_table, n = Inf)

write_csv(
  all_pairwise_results_table,
  "all_pairwise_results_withinconn_lmer_siteRE_conditional_predictions.csv"
)


# ============================================================

# new test just lm

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(readxl)

# ============================================================
# 1. Load data
# ============================================================

demos_withinconn <- read.csv("sheets/v1.3/demos_conn_2306.csv")

# Drop other diagnosis group if needed
demos_withinconn <- demos_withinconn %>%
  filter(!is.na(dcfdx), dcfdx != "other")

print(colSums(is.na(demos_withinconn)))

# ============================================================
# 2. Settings
# ============================================================

target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

covariates_to_plot <- c(
  "msex",
  "eyes",
  "syn_bin",
  "dcfdx"
)

model_formula <- within_conn ~ mean_FD + msex + site + age_scandate +
  eyes + dcfdx + syn_bin

fixed_y_limits <- c(0, 0.6)
fixed_y_breaks <- seq(0, 0.6, by = 0.25)

# ============================================================
# 3. Long-format data
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
      site = factor(site),
      msex = factor(msex),
      eyes = factor(eyes),
      syn_bin = factor(syn_bin),
      dcfdx = factor(dcfdx)
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

data_long <- make_long(demos_withinconn)

# ============================================================
# 4. Fit lm per network and get predicted values
# ============================================================

get_predicted_data_lm <- function(data_long) {
  
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

predicted_data <- get_predicted_data_lm(data_long)

plot_predicted_by_site_hue <- function(predicted_data, hue_var,
                                        y_limits = NULL,
                                        y_breaks = NULL) {
  
  dodge_width <- 0.75
  
  p <- ggplot(
    predicted_data,
    aes(
      x = site,
      y = predicted_conn,
      fill = .data[[hue_var]],
      color = .data[[hue_var]]
    )
  ) +
    geom_violin(
      position = position_dodge(width = dodge_width),
      alpha = 0.20,
      linewidth = 0.3,
      trim = FALSE,
      scale = "width"
    ) +
    geom_boxplot(
      position = position_dodge(width = dodge_width),
      width = 0.18,
      alpha = 0.65,
      outlier.shape = NA,
      linewidth = 0.35
    ) +
    geom_point(
      position = position_jitterdodge(
        jitter.width = 0.12,
        dodge.width = dodge_width
      ),
      alpha = 0.25,
      size = 0.55
    ) +
    facet_wrap(~ network, scales = "fixed") +
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0.12, 0.12))
    ) +
    coord_cartesian(ylim = y_limits) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste0("lm-predicted within-network connectivity by site, effect of ", hue_var),
      x = "Site",
      y = "Predicted within-network connectivity",
      fill = hue_var,
      color = hue_var
    )
  
  print(p)
  
  ggsave(
    filename = paste0(
      "withinconn_lm_predicted_by_site_hue_",
      hue_var,
      ".png"
    ),
    plot = p,
    width = 15,
    height = 9,
    dpi = 300
  )
  
  return(p)
}

p_sex <- plot_predicted_by_site_hue(
  predicted_data = predicted_data,
  hue_var = "msex",
  y_limits = fixed_y_limits,
  y_breaks = fixed_y_breaks
)

# ============================================================
#### site interaction

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(readxl)
library(broom)
library(emmeans)

# ============================================================
# 1. Load data
# ============================================================

demos_withinconn <- read.csv("demos_conn_2505.csv")

# Drop other diagnosis group if needed
demos_withinconn <- demos_withinconn %>%
  filter(!is.na(dcfdx), dcfdx != "other")

print(colSums(is.na(demos_withinconn)))

# ============================================================
# 2. Settings
# ============================================================

target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

covariates_to_plot <- c(
  "msex",
  "eyes",
  "syn_bin",
  "dcfdx"
)

# Main model:
# msex * site tests whether the msex effect differs across sites
model_formula_interaction <- within_conn ~ mean_FD + msex * site + age_scandate +
  eyes + dcfdx + syn_bin

# Reduced model:
# same model, but without the msex x site interaction
model_formula_no_interaction <- within_conn ~ mean_FD + msex + site + age_scandate +
  eyes + dcfdx + syn_bin

fixed_y_limits <- c(0, 0.6)
fixed_y_breaks <- seq(0, 0.6, by = 0.25)

alpha_level <- 0.05

# ============================================================
# 3. Long-format data
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
      site = factor(site),
      msex = factor(msex),
      eyes = factor(eyes),
      syn_bin = factor(syn_bin),
      dcfdx = factor(dcfdx, levels = c("NCI", "MCI", "AD"))
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

data_long <- make_long(demos_withinconn)

# Check sample size by site and msex
site_msex_counts <- data_long %>%
  count(network, site, msex) %>%
  arrange(network, site, msex)

print(site_msex_counts)

write.csv(
  site_msex_counts,
  "withinconn_site_msex_counts.csv",
  row.names = FALSE
)

# ============================================================
# 4. Fit models per network
# ============================================================

fit_models_by_network <- function(data_long) {
  
  map(target_cols, function(net) {
    
    df_net <- data_long %>%
      filter(network == net) %>%
      droplevels()
    
    model_no_interaction <- lm(
      model_formula_no_interaction,
      data = df_net
    )
    
    model_interaction <- lm(
      model_formula_interaction,
      data = df_net
    )
    
    list(
      network = net,
      data = df_net,
      model_no_interaction = model_no_interaction,
      model_interaction = model_interaction
    )
  }) %>%
    set_names(target_cols)
}

models_by_network <- fit_models_by_network(data_long)

# ============================================================
# 5. Omnibus interaction tests: does msex effect differ by site?
# ============================================================

interaction_tests <- map_dfr(models_by_network, function(x) {
  
  test <- anova(
    x$model_no_interaction,
    x$model_interaction
  )
  
  tibble(
    network = x$network,
    df_added = test$Df[2],
    rss_no_interaction = test$RSS[1],
    rss_interaction = test$RSS[2],
    f_value = test$F[2],
    p_value = test$`Pr(>F)`[2]
  )
}) %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "fdr"),
    interaction_significant_raw = p_value < alpha_level,
    interaction_significant_fdr = p_fdr < alpha_level
  ) %>%
  arrange(p_value)

print(interaction_tests)

write.csv(
  interaction_tests,
  "withinconn_msex_site_interaction_tests.csv",
  row.names = FALSE
)

# Networks to follow up
# I recommend using FDR-corrected significance.
# If you want raw p-value significance instead, replace with interaction_significant_raw.

significant_networks <- interaction_tests %>%
  filter(interaction_significant_fdr) %>%
  pull(network)

print(significant_networks)

# ============================================================
# 6. Add fitted values from interaction model
# ============================================================

get_predicted_data_lm <- function(models_by_network) {
  
  map_dfr(models_by_network, function(x) {
    
    x$data %>%
      mutate(
        predicted_conn = fitted(x$model_interaction)
      )
  }) %>%
    mutate(
      network = factor(network, levels = target_cols)
    )
}

predicted_data <- get_predicted_data_lm(models_by_network)

# ============================================================
# 7. Post hoc tests
# ============================================================
# These only run for networks where the msex x site interaction
# is significant.
#
# There are two useful post hoc questions:
#
# A) Within each site, is there an msex difference?
#    emmeans(model, pairwise ~ msex | site)
#
# B) Does the msex difference differ between pairs of sites?
#    emmeans(model, ~ msex | site), then contrast interaction = "pairwise"
#
# B is the direct follow-up to the interaction.

run_posthoc_tests <- function(models_by_network, significant_networks) {
  
  if (length(significant_networks) == 0) {
    message("No significant msex x site interactions. Skipping post hoc tests.")
    return(
      list(
        msex_within_site = tibble(),
        site_differences_in_msex_effect = tibble()
      )
    )
  }
  
  # A) msex effect within each site
  msex_within_site <- map_dfr(significant_networks, function(net) {
    
    model <- models_by_network[[net]]$model_interaction
    
    emm <- emmeans(model, ~ msex | site)
    
    pairs(emm, adjust = "fdr") %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(network = net) %>%
      relocate(network)
  })
  
  # B) compare the msex effect across sites
  site_differences_in_msex_effect <- map_dfr(significant_networks, function(net) {
    
    model <- models_by_network[[net]]$model_interaction
    
    emm <- emmeans(model, ~ msex | site)
    
    contrast(
      emm,
      interaction = "pairwise",
      by = NULL,
      adjust = "fdr"
    ) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(network = net) %>%
      relocate(network)
  })
  
  list(
    msex_within_site = msex_within_site,
    site_differences_in_msex_effect = site_differences_in_msex_effect
  )
}

posthoc_results <- run_posthoc_tests(
  models_by_network = models_by_network,
  significant_networks = significant_networks
)

msex_within_site_posthoc <- posthoc_results$msex_within_site
site_differences_in_msex_effect_posthoc <- posthoc_results$site_differences_in_msex_effect

print(msex_within_site_posthoc)
print(site_differences_in_msex_effect_posthoc)

write.csv(
  msex_within_site_posthoc,
  "withinconn_posthoc_msex_within_each_site.csv",
  row.names = FALSE
)

write.csv(
  site_differences_in_msex_effect_posthoc,
  "withinconn_posthoc_site_differences_in_msex_effect.csv",
  row.names = FALSE
)

# ============================================================
# 8. Estimated marginal means for plotting
# ============================================================

get_emm_msex_site <- function(models_by_network) {
  
  map_dfr(models_by_network, function(x) {
    
    emmeans(x$model_interaction, ~ msex | site) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(network = x$network)
  }) %>%
    mutate(
      network = factor(network, levels = target_cols),
      site = factor(site),
      msex = factor(msex)
    )
}

emm_msex_site <- get_emm_msex_site(models_by_network)

write.csv(
  emm_msex_site,
  "withinconn_emmeans_msex_by_site.csv",
  row.names = FALSE
)

# ============================================================
# 9. Plot fitted values by site and msex
# ============================================================

plot_predicted_by_site_hue <- function(predicted_data, hue_var,
                                        y_limits = NULL,
                                        y_breaks = NULL) {
  
  dodge_width <- 0.75
  
  p <- ggplot(
    predicted_data,
    aes(
      x = site,
      y = predicted_conn,
      fill = .data[[hue_var]],
      color = .data[[hue_var]]
    )
  ) +
    geom_violin(
      position = position_dodge(width = dodge_width),
      alpha = 0.20,
      linewidth = 0.3,
      trim = FALSE,
      scale = "width"
    ) +
    geom_boxplot(
      position = position_dodge(width = dodge_width),
      width = 0.18,
      alpha = 0.65,
      outlier.shape = NA,
      linewidth = 0.35
    ) +
    geom_point(
      position = position_jitterdodge(
        jitter.width = 0.12,
        dodge.width = dodge_width
      ),
      alpha = 0.25,
      size = 0.55
    ) +
    facet_wrap(~ network, scales = "fixed") +
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0.12, 0.12))
    ) +
    coord_cartesian(ylim = y_limits) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste0(
        "LM-predicted within-network connectivity by site and ",
        hue_var,
        " from msex x site interaction model"
      ),
      x = "Site",
      y = "Predicted within-network connectivity",
      fill = hue_var,
      color = hue_var
    )
  
  print(p)
  
  ggsave(
    filename = paste0(
      "withinconn_lm_predicted_by_site_hue_",
      hue_var,
      "_msex_site_interaction.png"
    ),
    plot = p,
    width = 15,
    height = 9,
    dpi = 300
  )
  
  return(p)
}

p_sex <- plot_predicted_by_site_hue(
  predicted_data = predicted_data,
  hue_var = "msex",
  y_limits = fixed_y_limits,
  y_breaks = fixed_y_breaks
)

# ============================================================
# 10. Plot estimated marginal means
# ============================================================

plot_emm_msex_site <- function(emm_msex_site,
                               y_limits = NULL,
                               y_breaks = NULL) {
  
  p <- ggplot(
    emm_msex_site,
    aes(
      x = site,
      y = emmean,
      color = msex,
      group = msex
    )
  ) +
    geom_point(
      position = position_dodge(width = 0.4),
      size = 2
    ) +
    geom_errorbar(
      aes(
        ymin = lower.CL,
        ymax = upper.CL
      ),
      position = position_dodge(width = 0.4),
      width = 0.15,
      linewidth = 0.35
    ) +
    geom_line(
      position = position_dodge(width = 0.4),
      linewidth = 0.5
    ) +
    facet_wrap(~ network, scales = "fixed") +
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0.12, 0.12))
    ) +
    coord_cartesian(ylim = y_limits) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Estimated marginal means for msex by site",
      subtitle = "From models including msex x site interaction",
      x = "Site",
      y = "Estimated marginal mean within-network connectivity",
      color = "msex"
    )
  
  print(p)
  
  ggsave(
    filename = "withinconn_emmeans_msex_by_site_interaction.png",
    plot = p,
    width = 15,
    height = 9,
    dpi = 300
  )
  
  return(p)
}

p_emm_msex_site <- plot_emm_msex_site(
  emm_msex_site = emm_msex_site,
  y_limits = fixed_y_limits,
  y_breaks = fixed_y_breaks
)

# ============================================================
# 11. Optional: plot only significant networks
# ============================================================

if (length(significant_networks) > 0) {
  
  emm_msex_site_sig <- emm_msex_site %>%
    filter(network %in% significant_networks) %>%
    droplevels()
  
  p_emm_msex_site_sig <- plot_emm_msex_site(
    emm_msex_site = emm_msex_site_sig,
    y_limits = fixed_y_limits,
    y_breaks = fixed_y_breaks
  )
  
  ggsave(
    filename = "withinconn_emmeans_msex_by_site_significant_networks_only.png",
    plot = p_emm_msex_site_sig,
    width = 12,
    height = 7,
    dpi = 300
  )
}


interaction_terms <- map_dfr(target_cols, function(net) {
  
  df_net <- data_long %>%
    filter(network == net) %>%
    droplevels()
  
  model <- lm(
    within_conn ~ mean_FD + msex * site + age_scandate +
      eyes + dcfdx + syn_bin,
    data = df_net
  )
})
  
coef(summary(model))[grep("msex.*:siteRIRC|siteRIRC.*:msex", 
                          rownames(coef(summary(model)))), ]

print(interaction_terms, n = Inf)

df_net <- data_long %>%
    filter(network == "Limbic") %>%
    droplevels()
  
  model <- lm(
    within_conn ~ mean_FD + msex * site + age_scandate +
      eyes + dcfdx + syn_bin,
    data = df_net
  )



interaction_terms_RIRC <- map_dfr(target_cols, function(net) {
  
  df_net <- data_long %>%
    filter(network == net) %>%
    droplevels()
  
  model <- lm(
    within_conn ~ mean_FD + msex + age_scandate +
      eyes * site + dcfdx + syn_bin,
    data = df_net
  )
  
  coef_table <- coef(summary(model))
  
  rirc_rows <- coef_table[
    grep(
      "msex.*:siteRIRC|siteRIRC.*:msex",
      rownames(coef_table)
    ),
    ,
    drop = FALSE
  ]
  
  as.data.frame(rirc_rows) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(
      network = net
    ) %>%
    relocate(network)
}) %>%
  mutate(
    p_fdr = p.adjust(`Pr(>|t|)`, method = "fdr"),
    p_stars_raw = case_when(
      `Pr(>|t|)` < 0.001 ~ "***",
      `Pr(>|t|)` < 0.01  ~ "**",
      `Pr(>|t|)` < 0.05  ~ "*",
      `Pr(>|t|)` < 0.1   ~ ".",
      TRUE               ~ ""
    ),
    p_stars_fdr = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01  ~ "**",
      p_fdr < 0.05  ~ "*",
      p_fdr < 0.1   ~ ".",
      TRUE          ~ ""
    ),
    significant_raw = `Pr(>|t|)` < 0.05,
    significant_fdr = p_fdr < 0.05
  ) %>%
  arrange(p_fdr)

print(interaction_terms_RIRC, n = Inf)

print(interaction_terms_RIRC, n = Inf)


interaction_terms_all_sites <- map_dfr(target_cols, function(net) {
  
  df_net <- data_long %>%
    filter(network == net) %>%
    droplevels()
  
  model <- lm(
    within_conn ~ mean_FD + msex + age_scandate +
      eyes + dcfdx + syn_bin * site,
    data = df_net
  )
  
  coef_table <- coef(summary(model))
  
  interaction_rows <- coef_table[
    grep(
      "syn_bin.*:site|site.*:syn_bin",
      rownames(coef_table)
    ),
    ,
    drop = FALSE
  ]
  
  as.data.frame(interaction_rows) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(network = net) %>%
    relocate(network)
}) %>%
  mutate(
    site = str_extract(term, "site[^:]+"),
    site = str_remove(site, "^site"),
    p_fdr = p.adjust(`Pr(>|t|)`, method = "fdr"),
    p_stars_raw = case_when(
      `Pr(>|t|)` < 0.001 ~ "***",
      `Pr(>|t|)` < 0.01  ~ "**",
      `Pr(>|t|)` < 0.05  ~ "*",
      `Pr(>|t|)` < 0.1   ~ ".",
      TRUE               ~ ""
    ),
    p_stars_fdr = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01  ~ "**",
      p_fdr < 0.05  ~ "*",
      p_fdr < 0.1   ~ ".",
      TRUE          ~ ""
    ),
    significant_raw = `Pr(>|t|)` < 0.05,
    significant_fdr = p_fdr < 0.05
  ) %>%
  relocate(site, .after = network) %>%
  arrange(site, network)

print(interaction_terms_all_sites, n = Inf)

interaction_terms_all_sites %>%
  group_by(site) %>%
  group_walk(~ {
    cat("\n\n==============================\n")
    cat("Site:", .y$site, "\n")
    cat("==============================\n")
    print(.x, n = Inf)
  })

levels(data_long$site)
