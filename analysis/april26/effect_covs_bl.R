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
# This script:
#   1. Loads within-network connectivity data.
#   2. Keeps only the earliest session per subject.
#   3. Reshapes the data from wide to long format.
#   4. Fits one adjusted linear model per network.
#   5. Uses emmeans to compute model-based pairwise comparisons
#      for categorical covariates such as sex, site, and eyes.
#   6. Plots raw distributions with violin + boxplot + jitter.
#   7. Adds significance brackets from model-based pairwise tests.
#
# Important:
#   The plotted points/violins/boxes are raw data.
#   The significance brackets are not raw tests; they come from
#   adjusted models using emmeans.

demos_withinconn <- read.csv("demos_withinconn.csv")
# filter df by BNK site and makee neew df containing only code and site cols
#bnks <- demos_withinconn %>%
#  filter(site == "BNK") %>%
#  select(sub_ses, site)
#write_csv(bnks, "BNK_sub_ses.csv")

# add syn_bin to demos_withinconn from column distortion_correction
demos_withinconn <- demos_withinconn %>%
  mutate(syn_bin = ifelse(distortion_correction == "SyN", 1, 0))

variables = read_excel("variables_ses_specific_may26.xlsx") %>%
  select(sub_id, ses_id, dcfdx)

# make . entry in dcfdx column to be NA
variables <- variables %>%
  mutate(dcfdx = ifelse(dcfdx == ".", NA, dcfdx))

# transform dcfdx 1 into NCI, 2 into CI, 3 into MCI, 4 into AD
variables <- variables %>%
  mutate(dcfdx = case_when(
    dcfdx == "1" ~ "NCI",
    dcfdx == "2" ~ "MCI",
    dcfdx == "3" ~ "MCI",
    dcfdx == "4" ~ "AD",
    dcfdx == "5" ~ "AD",
    dcfdx == "6" ~ "other",
    TRUE ~ as.character(dcfdx)
  ))

# merge the two dfs by sub_id and ses_id 
demos_withinconn <- demos_withinconn %>%
  left_join(variables, by = c("sub_id", "ses_id")) 

##
# ------------------------------------------------------------
# Keep only the lowest/earliest session per subject
# ------------------------------------------------------------
# ses_id is assumed to look like "ses-1", "ses-2", etc.
# We extract the numeric session value so ordering is numeric,
# not alphabetical.
#
# Example:
#   "ses-10" should come after "ses-2"
#   numeric extraction makes this behave correctly.

demos_withinconn <- demos_withinconn %>%
  mutate(
    ses_num = as.numeric(str_extract(ses_id, "\\d+"))
  ) %>%
  group_by(sub_id) %>%
  arrange(ses_num) %>%
  slice(1) %>%
  ungroup()

# Basic missingness check.
# This prints the number of missing values per column.
print(colSums(is.na(demos_withinconn)))

# Network columns to analyze.
# These columns are treated as outcome variables.
target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

# Motion threshold used for the post-filtering analysis.
fd_threshold <- 0.25

# Ensure variables have the correct type before modeling.
# Numeric variables should be numeric; categorical variables should be factors.
demos_withinconn <- demos_withinconn %>%
  mutate(
    mean_FD = as.numeric(mean_FD),
    msex = factor(
      msex,
      levels = c(0, 1),
      labels = c("female", "male")
    ),
    site = factor(site),
    age_scandate = as.numeric(age_scandate),
    distortion_correction = factor(distortion_correction),
    eyes = factor(eyes),
    dcfdx = factor(dcfdx),
    syn_bin = factor(syn_bin)
  )

# Color palette for network-specific plots.
network_colors <- c(
  "Vis" = "#9B59B6",
  "SomMot" = "#6C8EBF",
  "Default" = "#D36B78",
  "Limbic" = "#C9D39A",
  "DorsAttn" = "#3C8D2F",
  "SalVentAttn" = "#C84CCF",
  "Cont" = "#E5B53A"
)

# Categorical covariates to test and plot.
# Each of these will be analyzed with emmeans pairwise comparisons.
covariates_to_run <- c(
  "msex",
  "site",
  "eyes",
  "syn_bin",
  "dcfdx"
)

# ------------------------------------------------------------
# Adjusted model formula
# ------------------------------------------------------------
# This is the model used for emmeans pairwise comparisons.
#
# Outcome is added later after pivoting:
#   within_conn
#
# Predictors:
#   mean_FD       = motion
#   msex          = sex
#   site          = scanner/site
#   age_scandate  = age at scan
#   eyes          = eyes open/closed
#
# Current model does NOT include distortion dummy columns.
# If you want to include selected distortion dummy columns, you can
# add them dynamically using distortion_dummy_cols.

model_formula <- as.formula(
  paste(
    "within_conn ~ mean_FD + msex + site + age_scandate + eyes + dcfdx + syn_bin"
  )
)

# ============================================================
# 2. Helper functions
# ============================================================

#' Convert wide network data to long format.
#'
#' @param data A dataframe containing one row per subject/session and
#'   one column per network in `target_cols`.
#'
#' @return A long-format dataframe with:
#'   - network: network name, e.g. Vis, SomMot, Default
#'   - within_conn: within-network connectivity value
#'   - all original covariates retained
#'
#' Details:
#'   The original data have separate columns for each network.
#'   For modeling/plotting across networks, it is cleaner to reshape
#'   to one row per subject per network.
#'
#'   Missing values in required modeling variables are removed.
#'   `droplevels()` removes factor levels that are no longer present
#'   after filtering.
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
      !is.na(dcfdx)
    ) %>%
    droplevels()
}

#' Convert p-values to significance stars.
#'
#' @param p Numeric vector of p-values.
#'
#' @return Character vector with:
#'   - "***" for p < 0.001
#'   - "**"  for p < 0.01
#'   - "*"   for p < 0.05
#'   - ""    otherwise
sig_from_p <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

# ============================================================
# 3. Model-based pairwise comparisons
# ============================================================

#' Fit adjusted models and compute pairwise emmeans comparisons.
#'
#' @param data_long Long-format dataframe produced by `make_long()`.
#' @param covariate Character string naming a categorical covariate
#'   to compare, e.g. "msex", "site", or "eyes".
#'
#' @return A dataframe with pairwise contrasts for each network.
#'
#' Details:
#'   For each network, this function fits:
#'
#'     within_conn ~ mean_FD + msex + site + age_scandate + eyes
#'
#'   Then it uses:
#'
#'     emmeans(model, specs = pairwise ~ covariate)
#'
#'   to compare levels of the requested covariate while accounting for
#'   the other terms in the model.
#'
#' Interpretation:
#'   The resulting contrasts are model-based pairwise comparisons.
#'   For example, for site:
#'
#'     siteA - siteB
#'
#'   means the estimated marginal mean for siteA minus the estimated
#'   marginal mean for siteB, adjusted for the other covariates.
#'

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
# 4. Prepare bracket annotations
# ============================================================

#' Convert emmeans pairwise results into bracket annotations.
#'
#' @param data_long Long-format dataframe used for plotting.
#' @param pairwise_results Output of `fit_model_pairwise()`.
#' @param covariate Character string naming the covariate being plotted.
#'
#' @return A dataframe formatted for `ggpubr::stat_pvalue_manual()`.
#'
#' Details:
#'   `stat_pvalue_manual()` needs columns named:
#'     - group1
#'     - group2
#'     - y.position
#'     - label
#'
#'   This function:
#'     1. Keeps only significant pairwise comparisons.
#'     2. Splits contrasts like "A - B" into group1 = "A", group2 = "B".
#'     3. Cleans cases where emmeans names levels like "msex0" instead
#'        of just "0".
#'     4. Calculates y positions separately for each network facet.
#'     5. Alternates brackets above and below the data to reduce clutter.
#'
#' Note:
#'   Stars are based on raw p-values because no correction is applied.
make_pairwise_brackets <- function(data_long, pairwise_results, covariate) {
  
  x_levels <- levels(factor(data_long[[covariate]]))
  
  y_positions <- data_long %>%
    group_by(network) %>%
    summarise(
      y_max = max(within_conn, na.rm = TRUE),
      y_min = min(within_conn, na.rm = TRUE),
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
# 5. Plot function
# ============================================================

#' Plot raw distributions with model-based pairwise significance brackets.
#'
#' @param data_long Long-format dataframe.
#' @param covariate Character string naming the categorical covariate
#'   to plot on the x-axis.
#' @param pairwise_results Output of `fit_model_pairwise()`.
#' @param title Plot title.
#'
#' @return A ggplot object.
#'
#' Plot layers:
#'   - Violin plot: distribution shape
#'   - Boxplot: median and interquartile range
#'   - Jittered points: individual raw observations
#'   - Significance brackets: model-based pairwise comparisons
#'
#' Important:
#'   The raw data layers are not adjusted.
#'   The significance brackets are adjusted/model-based because they
#'   come from emmeans applied to the linear model.
plot_factor_covariate <- function(data_long, covariate, pairwise_results, title) {
  
  pairwise_annot <- make_pairwise_brackets(
    data_long = data_long,
    pairwise_results = pairwise_results,
    covariate = covariate
  )
  
  p <- ggplot(
    data_long,
    aes(
      x = .data[[covariate]],
      y = within_conn,
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
      outliers = FALSE,
      outlier.shape = NA,
      linewidth = 0.35
    ) +
    geom_jitter(
      width = 0.12,
      alpha = 0.09,
      size = 0.4
    ) +
    facet_wrap(~ network, scales = "fixed") +
    scale_color_manual(values = network_colors) +
    scale_fill_manual(values = network_colors) +
    scale_y_continuous(
      expand = expansion(mult = c(0.18, 0.18))
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
      subtitle = "Stars show HSD-corrected model-based pairwise comparisons from emmeans",
      x = covariate,
      y = "Within-network connectivity"
    )
  
  # Add brackets only when at least one significant comparison exists.
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
# 6. Print results
# ============================================================

#' Print pairwise comparison results in a readable table.
#'
#' @param pairwise_results Output of `fit_model_pairwise()`.
#' @param title Text header printed before the table.
#'
#' @return Invisibly prints a formatted tibble.
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
# 7. Wrapper
# ============================================================

#' Run the full analysis for one categorical covariate.
#'
#' @param data_long Long-format dataframe.
#' @param covariate Character string naming the covariate to analyze.
#' @param analysis_label Text label used in plot titles and printed output.
#' @param file_suffix Text suffix used in the saved PNG filename.
#'
#' @return A list containing:
#'   - pairwise: pairwise emmeans results
#'   - plot: ggplot object
#'
#' Workflow:
#'   1. Fit models and compute emmeans pairwise comparisons.
#'   2. Print pairwise comparison table.
#'   3. Plot raw distributions with model-based stars.
#'   4. Save plot to disk.
run_factor_analysis <- function(data_long, covariate, analysis_label, file_suffix) {
  
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
    title = paste0(analysis_label, ": raw distribution by ", covariate, " at baseline")
  )
  
  print(p)
  
  ggsave(
    filename = paste0("withinconn_", covariate, "_", file_suffix, ".png"),
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

# ============================================================
# 8. Pre-filtering analysis
# ============================================================

# Convert to long format after all preprocessing.
data_long_pre <- make_long(demos_withinconn)

# Run the full analysis for each covariate.
# `.x` means the current covariate name from `covariates_to_run`.
pre_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    data_long = data_long_pre,
    covariate = .x,
    analysis_label = "Pre-filtering",
    file_suffix = "preFDfilter"
  )
)

names(pre_results) <- covariates_to_run

# ============================================================
# 9. Post-filtering analysis: mean_FD < 0.25
# ============================================================

# Repeat the same analysis after excluding higher-motion scans.
data_long_post <- data_long_pre %>%
  filter(mean_FD < fd_threshold) %>%
  droplevels()

post_results <- map(
  covariates_to_run,
  ~ run_factor_analysis(
    data_long = data_long_post,
    covariate = .x,
    analysis_label = paste0("Post-filtering, mean_FD < ", fd_threshold),
    file_suffix = "postFDfilter"
  )
)

names(post_results) <- covariates_to_run

# ============================================================
# 10. Combine all pairwise results
# ============================================================

# Combine pre- and post-filtering pairwise tables into one final table.
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



#data_long_pre <- make_long(demos_withinconn)

#df_net <- data_long_pre %>%
#      filter(network == "Vis")



#model_formula <- as.formula(
#  paste(
#    "within_conn ~ mean_FD + msex + site + age_scandate + eyes"
#  )
#)
#model <- lm(
#      model_formula,
#      data = df_net
#    )   
#emm <- emmeans(
#      model,
#      specs = pairwise ~ eyes
#    )

#print(emm)

# emm$emmeans: adjusted means per site

# This gives the model-adjusted estimated mean connectivity for each site.

# emm$contrasts: pairwise comparisons between sites, with estimates, SE, t/z ratios, and p-values.

