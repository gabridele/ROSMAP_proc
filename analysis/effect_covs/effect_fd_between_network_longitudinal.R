# Longitudinal motion-sensitivity analysis for between-network connectivity.
# Uses repeated sessions and participant-level random intercepts to estimate
# the adjusted association between mean framewise displacement and connectivity.

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
library(ggeffects)
library(readxl)
library(grDevices)
library(lme4)
library(lmerTest)

# ============================================================
# 1. Load and prepare data
# ============================================================

demos_betweenconn <- read.csv(require_file(demos("demos_conn.csv")))

# Add numeric session 
demos_betweenconn <- demos_betweenconn %>%
  mutate(
    ses_num = as.numeric(str_extract(ses_id, "\\d+"))
  )


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

demos_betweenconn <- demos_betweenconn %>%
  mutate(
    sub_id = factor(sub_id),
    mean_FD = as.numeric(mean_FD),
    msex = factor(msex),
    site = factor(site),
    age_scandate = as.numeric(age_scandate),
    distortion_correction = factor(distortion_correction),
    syn_bin = factor(syn_bin),
    eyes = factor(eyes),
    dcfdx = factor(dcfdx)
  )

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
  syn_bin + eyes + dcfdx + (1 | sub_id)

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
    
    model_adj <- lmer(
      model_formula,
      data = df_net
    )
    
    tibble(
      network_combo = net,
      beta_adjusted = fixef(model_adj)["mean_FD"],
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
    
  model_adj <- lmer(
      between_conn ~ mean_FD + msex + site + age_scandate +
  syn_bin + eyes + dcfdx + (1 | sub_id),
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

data_long_pre <- make_long(demos_betweenconn)

model_results_pre <- fit_fd_models(data_long_pre)
pred_adjusted_pre <- get_adjusted_predictions(data_long_pre)

print_model_table(
  model_results_pre,
  "Pre-filtering adjusted FD model results"
)
# save table 
write_csv(model_results_pre, output("motion", "fd_effects_betweenconn_long_preFDfilter_modelresults.csv"))

p_pre <- plot_fd_effects(
  data_long_pre,
  pred_adjusted_pre,
  model_results_pre,
  "Adjusted FD effect on between-network connectivity: longitudinal"
)

print(p_pre)

ggsave(
  output("motion", "fd_effects_betweenconn_long_preFDfilter.png"),
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
summary(pred_adjusted_post)
print_model_table(
  model_results_post,
  paste0("Post-filtering adjusted FD model results: longitudinal, mean_FD < ", fd_threshold)
)
# save table
write_csv(model_results_post, output("motion", "fd_effects_betweenconn_long_postFDfilter_modelresults.csv"))

p_post <- plot_fd_effects(
  data_long_post,
  pred_adjusted_post,
  model_results_post,
  paste0("Adjusted FD effect on between-network connectivity: longitudinal, mean_FD < ", fd_threshold)
)

print(p_post)

ggsave(
  output("motion", "fd_effects_betweenconn_long_postFDfilter.png"),
  plot = p_post,
  width = 14,
  height = 10,
  dpi = 300
)
