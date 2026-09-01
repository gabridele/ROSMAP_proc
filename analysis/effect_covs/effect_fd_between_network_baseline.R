# Baseline motion-sensitivity analysis for between-network connectivity.
# Keeps the earliest session per participant and fits adjusted linear models
# relating mean framewise displacement (FD) to each network-pair outcome.

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

# Missingness sanity check
print(colSums(is.na(demos_betweenconn)))

# Keep lowest session per subject
demos_min_ses <- demos_betweenconn %>%
  group_by(sub_id) %>%
  arrange(ses_num) %>%
  slice(1) %>%
  ungroup()

fd_threshold <- FD_THRESHOLD

demos_min_ses <- demos_min_ses %>%
  mutate(
    mean_FD = as.numeric(mean_FD),
    msex = factor(msex),
    site = factor(site),
    age_scandate = as.numeric(age_scandate),
    syn_bin = factor(syn_bin),
    eyes = factor(eyes),
    dcfdx = factor(dcfdx)
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
  syn_bin + eyes + dcfdx

# ============================================================
# 4. Helper functions
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
      !is.na(mean_FD),
      !is.na(between_conn),
      !is.na(msex),
      !is.na(site),
      !is.na(age_scandate),
      !is.na(syn_bin),
      !is.na(eyes),
      !is.na(dcfdx)
    ) %>%
    droplevels()
}

fit_fd_models <- function(data_long) {
  
  map_dfr(target_combos, function(net) {
    
    df_net <- data_long %>%
      filter(network_combo == net) %>%
      droplevels()
    
    model_adj <- lm(
      model_formula,
      data = df_net
    )
    
    tibble(
      network_combo = net,
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
      network_combo = factor(network_combo, levels = target_combos)
    )
}

get_adjusted_predictions <- function(data_long) { 
  
  map_dfr(target_combos, function(net) { 
    
  df_net <- data_long %>% filter(network_combo == net) 
    
  model_adj <- lm(
      between_conn ~ mean_FD + msex + site + age_scandate +
  syn_bin + eyes + dcfdx,
      data = df_net
    )
    
  predict_response( 
    model_adj, 
    terms = "mean_FD [all]" 
    ) %>% as.data.frame() %>% mutate(network_combo = net) 
  }) %>% 
    mutate(network_combo = factor(network_combo, levels = target_combos)) 
}

plot_fd_effects <- function(data_long, pred_adjusted, model_results, title) {
  
  label_pos <- data_long %>%
    group_by(network_combo) %>%
    summarise(
      x = max(mean_FD, na.rm = TRUE),
      y = max(between_conn, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(model_results, by = "network_combo")
  
  ggplot(data_long, aes(x = mean_FD, y = between_conn, color = network_combo)) +
    geom_point(alpha = 0.30, size = 1) +
    
    geom_ribbon(
      data = pred_adjusted,
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = network_combo),
      inherit.aes = FALSE,
      alpha = 0.22
    ) +
    
    geom_line(
      data = pred_adjusted,
      aes(x = x, y = predicted, color = network_combo),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    
    geom_label(
      data = label_pos,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 2.5,
      color = "black",
      fill = "white",
      alpha = 0.85,
      label.size = 0
    ) +
    
    facet_wrap(~ network_combo, scales = "fixed") +
    scale_color_manual(values = between_network_colors) +
    scale_fill_manual(values = between_network_colors) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 8),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = title,
      x = "mean_FD",
      y = "Between-network connectivity"
    )
}

print_model_table <- function(model_results, title) {
  
  cat("\n============================================================\n")
  cat(title, "\n")
  cat("============================================================\n")
  
  model_results %>%
    select(network_combo, beta_adjusted, t_val_adjusted, p_adjusted, q_adjusted, sig_adjusted) %>%
    mutate(
      beta_adjusted = round(beta_adjusted, 4),
      t_val_adjusted = round(t_val_adjusted, 3),
      p_adjusted = signif(p_adjusted, 3),
      q_adjusted = signif(q_adjusted, 3)
    ) %>%
    as_tibble() %>%
    print(n = Inf)
}

# ============================================================
# 5. Pre-filtering analysis
# ============================================================

data_long_pre <- make_long(demos_min_ses)

model_results_pre <- fit_fd_models(data_long_pre)
pred_adjusted_pre <- get_adjusted_predictions(data_long_pre)

print_model_table(
  model_results_pre,
  "Pre-filtering adjusted FD model results"
)
# save table 
write_csv(model_results_pre, output("motion", "fd_effects_betweenconn_bl_preFDfilter_modelresults.csv"))

p_pre <- plot_fd_effects(
  data_long_pre,
  pred_adjusted_pre,
  model_results_pre,
  "Adjusted FD effect on between-network connectivity: baseline only"
)

print(p_pre)

ggsave(
  output("motion", "fd_effects_betweenconn_bl.png"),
  plot = p_pre,
  width = 14,
  height = 10,
  dpi = 300
)

# ============================================================
# 6. Post-filtering analysis: mean_FD < 0.25
# ============================================================

data_long_post <- data_long_pre %>%
  filter(mean_FD < fd_threshold) %>%
  droplevels()

model_results_post <- fit_fd_models(data_long_post)
pred_adjusted_post <- get_adjusted_predictions(data_long_post)

print_model_table(
  model_results_post,
  paste0("Post-filtering adjusted FD model results: mean_FD < ", fd_threshold)
)
# save table
write_csv(model_results_post, output("motion", "fd_effects_betweenconn_bl_postFDfilter_modelresults.csv"))

p_post <- plot_fd_effects(
  data_long_post,
  pred_adjusted_post,
  model_results_post,
  paste0("Adjusted FD effect on between-network connectivity: baseline, mean_FD < ", fd_threshold)
)

print(p_post)

ggsave(
  output("motion", "fd_effects_betweenconn_bl_postFDfilter.png"),
  plot = p_post,
  width = 14,
  height = 10,
  dpi = 300
)