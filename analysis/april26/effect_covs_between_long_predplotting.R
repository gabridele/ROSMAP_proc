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

emm_options(pbkrtest.limit = 10000)
emm_options(lmerTest.limit = 50000)

# ============================================================
# 1. Load and prepare data
# ============================================================

demos_betweenconn <- read.csv("demos_conn_1905.csv")

# Add numeric session
demos_betweenconn <- demos_betweenconn %>%
  mutate(
    ses_num = as.numeric(str_extract(ses_id, "\\d+"))
  )

## Add binary SyN variable
#demos_betweenconn <- demos_betweenconn %>%
#  mutate(
#    syn_bin = ifelse(distortion_correction == "SyN", 1, 0)
#  )
#
## Load diagnosis/session-specific variables
#variables <- read_excel("variables_ses_specific_may26.xlsx") %>%
#  select(sub_id, ses_id, dcfdx) %>%
#  mutate(
#    dcfdx = ifelse(dcfdx == ".", NA, dcfdx),
#    dcfdx = case_when(
#      dcfdx == "1" ~ "NCI",
#      dcfdx == "2" ~ "MCI",
#      dcfdx == "3" ~ "MCI",
#      dcfdx == "4" ~ "AD",
#      dcfdx == "5" ~ "AD",
#      dcfdx == "6" ~ "other",
#      TRUE ~ as.character(dcfdx)
#    )
#  )
#
## Merge diagnosis into main dataframe
#demos_betweenconn <- demos_betweenconn %>%
#  left_join(variables, by = c("sub_id", "ses_id"))
#
## Keep lowest/earliest session per subject
#demos_betweenconn <- demos_betweenconn %>%
#  mutate(
#    ses_num = as.numeric(str_extract(ses_id, "\\d+"))
#  ) %>%
#  group_by(sub_id) %>%
#  arrange(ses_num) %>%
#  slice(1) %>%
#  ungroup()

# Missingness check
print(colSums(is.na(demos_betweenconn)))

# Between-network columns
target_combos <- c(
  "Cont_to_Default", "Cont_to_DorsAttn",
  "Cont_to_Limbic", "Cont_to_SalVentAttn", "Cont_to_SomMot",
  "Cont_to_Vis", "Default_to_DorsAttn", "Default_to_Limbic",
  "Default_to_SalVentAttn", "Default_to_SomMot", "Default_to_Vis",
  "DorsAttn_to_Limbic", "DorsAttn_to_SalVentAttn", "DorsAttn_to_SomMot",
  "DorsAttn_to_Vis", "Limbic_to_SalVentAttn", "Limbic_to_SomMot",
  "Limbic_to_Vis", "SalVentAttn_to_SomMot", "SalVentAttn_to_Vis",
  "SomMot_to_Vis"
)

fd_threshold <- 0.25

# Type conversion
#demos_betweenconn <- demos_betweenconn %>%
#  mutate(
#    mean_FD = as.numeric(mean_FD),
#    msex = factor(
#      msex,
#      levels = c(0, 1),
#      labels = c("female", "male")
#    ),
#    site = factor(site),
#    age_scandate = as.numeric(age_scandate),
#    eyes = factor(eyes),
#    dcfdx = factor(
#      dcfdx,
#      levels = c("NCI", "MCI", "AD", "other")
#    ),
#    syn_bin = factor(
#      syn_bin,
#      levels = c(0, 1),
#      labels = c("not SyN", "SyN")
#    )
#  )
#
### add scandate
#demos_betweenconn <- read_csv("age_atscan.csv") %>%
#  separate(col = "scandate_visit_projID", into = c("scandate", "visit", "projID"), sep = "_") %>%
#  select(c("ses_id", "sub_id", "scandate")) %>%
#  right_join(demos_betweenconn, by = c("sub_id", "ses_id"))
#
## make scandate format yyyymmdd into a datee
#demos_betweenconn <- demos_betweenconn %>%
#  mutate(scandate = as.Date(as.character(scandate), format = "%Y%m%d"))
#
#demos_betweenconn <- demos_betweenconn %>%
#  mutate(
#    ses_num = as.numeric(str_extract(ses, "\\d+"))
#  ) %>%
#  mutate(sub = factor(sub))
#
## compute years from baseline for each subject
#demos_betweenconn <- demos_betweenconn %>%
#  group_by(sub_id) %>%
#  mutate(
#    baseline_date = scandate[which.min(ses_num)],
#    years_from_baseline = interval(baseline_date, scandate) / years(1)
#  ) %>%
#  ungroup()

# ============================================================
# 2. Colors for between-network combos
# ============================================================

between_network_colors <- c(
  "Cont_to_Default" = "#DC9059",
  "Cont_to_DorsAttn" = "#91A135",
  "Cont_to_Limbic" = "#D7C46A",
  "Cont_to_SalVentAttn" = "#D78185",
  "Cont_to_SomMot" = "#A9A27D",
  "Cont_to_Vis" = "#C08778",
  "Default_to_DorsAttn" = "#887C54",
  "Default_to_Limbic" = "#CE9F89",
  "Default_to_SalVentAttn" = "#CE5CA3",
  "Default_to_SomMot" = "#A07D9C",
  "Default_to_Vis" = "#B76297",
  "DorsAttn_to_Limbic" = "#83B065",
  "DorsAttn_to_SalVentAttn" = "#826D7F",
  "DorsAttn_to_SomMot" = "#548E77",
  "DorsAttn_to_Vis" = "#6B7373",
  "Limbic_to_SalVentAttn" = "#C990B4",
  "Limbic_to_SomMot" = "#9BB1AD",
  "Limbic_to_Vis" = "#B296A8",
  "SalVentAttn_to_SomMot" = "#9A6DC7",
  "SalVentAttn_to_Vis" = "#B252C3",
  "SomMot_to_Vis" = "#8474BB"
)

# Make sure color vector is in the same order as target_combos
between_network_colors <- between_network_colors[target_combos]

# ============================================================
# 3. Model formula
# ============================================================

model_formula <- between_conn ~ mean_FD + msex + site + age_scandate +
  eyes + dcfdx + syn_bin + (1 | sub_id)

# Categorical covariates to test/plot
covariates_to_run <- c(
  "msex",
  "site",
  "eyes",
  "syn_bin",
  "dcfdx"
)

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
    filename = paste0("betweenconn_predicted_", covariate, "_", file_suffix, ".png"),
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
write_csv(all_pairwise_results_table, "all_pairwise_results_betweenconn_covs_long.csv")