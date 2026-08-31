# Baseline motion-sensitivity analysis for within-network connectivity.
# Keeps the earliest session per participant and compares adjusted FD effects
# before and after the prespecified mean-FD exclusion threshold.

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

# ============================================================
# 1. Load and prepare data
# ============================================================

demos_withinconn <- read.csv(require_file(demos("demos_conn.csv")))

# Missingness sanity check
print(colSums(is.na(demos_withinconn)))

# Keep the earliest numeric session per participant.
demos_min_ses <- demos_withinconn %>%
  mutate(ses_num = as.numeric(str_extract(ses_id, "\\d+"))) %>%
  group_by(sub_id) %>%
  arrange(ses_num, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    mean_FD = as.numeric(mean_FD),
    msex = factor(msex),
    site = factor(site),
    age_scandate = as.numeric(age_scandate),
    eyes = factor(eyes),
    syn_bin = factor(syn_bin),
    dcfdx = factor(dcfdx)
  )

target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

fd_threshold <- FD_THRESHOLD

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
# 2. Helper functions
# ============================================================

make_long <- function(data) {
  data %>%
    pivot_longer(
      cols = all_of(target_cols),
      names_to = "network",
      values_to = "within_conn"
    ) %>%
    mutate(network = factor(network, levels = target_cols)) %>%
    filter(
      !is.na(mean_FD),
      !is.na(within_conn),
      !is.na(msex),
      !is.na(site),
      !is.na(age_scandate),
      !is.na(distortion_correction),
      !is.na(eyes)
    )
}

fit_fd_models <- function(data_long) {
  
  map_dfr(target_cols, function(net) {
    
    df_net <- data_long %>%
      filter(network == net)
    
    model_adj <- lm(
      within_conn ~ mean_FD + msex + site + age_scandate +
        syn_bin + eyes + dcfdx,
      data = df_net
    )
    
    tibble(
      network = net,
      beta_adjusted = coef(model_adj)["mean_FD"],
      t_val_adjusted = summary(model_adj)$coefficients["mean_FD", "t value"],
      p_adjusted = summary(model_adj)$coefficients["mean_FD", "Pr(>|t|)"]
    )
  }) %>%
    mutate(
      q_adjusted = p.adjust(p_adjusted, method = "fdr"),
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
      network = factor(network, levels = target_cols)
    )
}

get_adjusted_predictions <- function(data_long) {
  
  map_dfr(target_cols, function(net) {
    
    df_net <- data_long %>%
      filter(network == net)
    
    model_adj <- lm(
      within_conn ~ mean_FD + msex + site + age_scandate +
        syn_bin + eyes + dcfdx,
      data = df_net
    )
    
    predict_response(
      model_adj,
      terms = "mean_FD [all]"
    ) %>%
      as.data.frame() %>%
      mutate(network = net)
  }) %>%
    mutate(network = factor(network, levels = target_cols))
}

plot_fd_effects <- function(data_long, pred_adjusted, model_results, title) {
  
  label_pos <- data_long %>%
    group_by(network) %>%
    summarise(
      x = max(mean_FD, na.rm = TRUE),
      y = max(within_conn, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(model_results, by = "network")
  
  ggplot(data_long, aes(mean_FD, within_conn, color = network)) +
    geom_point(alpha = 0.3, size = 1) +
    
    geom_ribbon(
      data = pred_adjusted,
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = network),
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    
    geom_line(
      data = pred_adjusted,
      aes(x = x, y = predicted, color = network),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    
    geom_text(
      data = label_pos,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 3,
      color = "black"
    ) +
    
    facet_wrap(~ network, scales = "fixed") +
    scale_color_manual(values = network_colors) +
    scale_fill_manual(values = network_colors) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = title,
      x = "mean_FD",
      y = "Within-network connectivity"
    )
}

print_model_table <- function(model_results, title) {
  
  cat("\n============================================================\n")
  cat(title, "\n")
  cat("============================================================\n")
  
  model_results %>%
    select(network, beta_adjusted, t_val_adjusted, p_adjusted, q_adjusted, sig_adjusted) %>%
    mutate(
      beta_adjusted = round(beta_adjusted, 4),
      t_val_adjusted = round(t_val_adjusted, 3),
      p_adjusted = signif(p_adjusted, 3),
      q_adjusted = signif(q_adjusted, 3)
    ) %>%
    print()
}

# ============================================================
# 3. Pre-filtering analysis
# ============================================================

data_long_pre <- make_long(demos_min_ses)

model_results_pre <- fit_fd_models(data_long_pre)
pred_adjusted_pre <- get_adjusted_predictions(data_long_pre)

print_model_table(
  model_results_pre,
  "Pre-filtering adjusted FD model results"
)
# save table
write_csv(model_results_pre, output("motion", "fd_effects_withinconn_bl_preFDfilter_modelresults.csv"))

p_pre <- plot_fd_effects(
  data_long_pre,
  pred_adjusted_pre,
  model_results_pre,
  "Adjusted FD effect by network: baseline only"
)

print(p_pre)

# save plot
ggsave(
  output("motion", "fd_effects_withinconn_bl.png"),
  plot = p_pre,
  width = 12,
  height = 8,
  dpi = 300
)

# ============================================================
# 4. Post-filtering analysis: mean_FD < 0.25
# ============================================================

data_long_post <- data_long_pre %>%
  filter(mean_FD < fd_threshold)

model_results_post <- fit_fd_models(data_long_post)
pred_adjusted_post <- get_adjusted_predictions(data_long_post)

print_model_table(
  model_results_post,
  paste0("Post-filtering adjusted FD model results: mean_FD < ", fd_threshold)
)

# save table
write_csv(model_results_post, output("motion", "fd_effects_withinconn_bl_postFDfilter_modelresults.csv"))

p_post <- plot_fd_effects(
  data_long_post,
  pred_adjusted_post,
  model_results_post,
  paste0("Adjusted FD effect by network: baseline only, mean_FD < ", fd_threshold)
)

print(p_post)

# save plot
ggsave(
  output("motion", "fd_effects_withinconn_bl_postFDfilter.png"),
  plot = p_post,
  width = 12,
  height = 8,
  dpi = 300
)
