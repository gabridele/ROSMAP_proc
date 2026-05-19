library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggeffects)
library(lme4)
library(lmerTest)

# ============================================================
# 1. Load and prepare data
# ============================================================

demos_withinconn <- read.csv("demos_conn_1905.csv")

demos_withinconn <- demos_withinconn %>%
  mutate(
    ses_num = as.numeric(str_extract(ses, "\\d+"))
  ) %>%
  mutate(sub = factor(sub))

target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

fd_threshold <- 0.25

demos_withinconn <- demos_withinconn %>%
  mutate(
    mean_FD = as.numeric(mean_FD),
    msex = factor(msex),
    site = factor(site),
    age_scandate = as.numeric(age_scandate),
    distortion_correction = factor(distortion_correction),
    eyes = factor(eyes)
  )

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
      !is.na(ses_num),
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
    
    model_adj <- lmer(
      within_conn ~ mean_FD + ses_num + msex + site + age_scandate + eyes + syn_bin + dcfdx + (1 | sub),
      data = df_net
    )
    
    tibble(
      network = net,
      beta_adjusted = fixef(model_adj)["mean_FD"],
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
    
    model_adj <- lmer(
      within_conn ~ mean_FD + ses_num + msex + site + age_scandate + eyes + syn_bin + dcfdx + (1 | sub),
      data = df_net
    )
    
    predict_response(
      model_adj,
      terms = "mean_FD [all]",
      type = "fixed"
    ) %>%
      as.data.frame() %>%
      mutate(network = net)
  }) %>%
    mutate(network = factor(network, levels = target_cols))
}

plot_fd_effects <- function(data_long, pred_adjusted, model_results, title) {
  
  label_pos <- model_results %>%
    mutate(
      x = Inf,
      y = Inf
    )
  
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
    
    geom_label(
      data = label_pos,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 3,
      color = "black",
      fill = "white",
      alpha = 0.8,
      linewidth = 0
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
      x = "mean FD",
      y = "Within-network connectivity"
    )
}

print_model_table <- function(model_results, title) {
  
  cat("\n============================================================\n")
  cat(title, "\n")
  cat("============================================================\n")
  
  model_results %>%
    select(network, beta_adjusted, p_adjusted, q_adjusted, sig_adjusted) %>%
    mutate(
      beta_adjusted = round(beta_adjusted, 4),
      p_adjusted = signif(p_adjusted, 3),
      q_adjusted = signif(q_adjusted, 3)
    ) %>%
    print()
}

# ============================================================
# 3. Pre-filtering analysis
# ============================================================

data_long_pre <- make_long(demos_withinconn)

model_results_pre <- fit_fd_models(data_long_pre)
pred_adjusted_pre <- get_adjusted_predictions(data_long_pre)

print_model_table(
  model_results_pre,
  "Pre-filtering adjusted FD model results"
)
# save table
write_csv(model_results_pre, "fd_effects_withinconn_longitudinal_preFDfilter_modelresults.csv")

p_pre <- plot_fd_effects(
  data_long_pre,
  pred_adjusted_pre,
  model_results_pre,
  "Adjusted FD effect by network (longitudinal)"
)

print(p_pre)

ggsave(
  "fd_effects_withinconn_longitudinal.png",
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
write_csv(model_results_post, "fd_effects_withinconn_longitudinal_postFDfilter_modelresults.csv")

p_post <- plot_fd_effects(
  data_long_post,
  pred_adjusted_post,
  model_results_post,
  paste0("Adjusted FD effect by network (longitudinal), mean_FD < ", fd_threshold)
)

print(p_post)

# save plot
ggsave(
  "fd_effects_withinconn_longitudinal_postFDfilter.png",
  plot = p_post,
  width = 12,
  height = 8,
  dpi = 300
)






######
# misc
# fixed-effect model matrix is rank deficient so dropping 2 columns / coefficients

#data_long_post %>%
#  summarise(
#    msex = n_distinct(msex),
#    site = n_distinct(site),
#    distortion = n_distinct(distortion_correction),
#    eyes = n_distinct(eyes)
#  )
#
#
#model_test <- lmer(
#  within_conn ~ mean_FD + ses_num + msex + site + age_scandate +
#    distortion_correction + eyes + (1 | sub),
#  data = data_long_post
#)
#
#alias(lm(
#  within_conn ~ mean_FD + ses_num + msex + site + age_scandate +
#    distortion_correction + eyes,
#  data = data_long_post
#))
#
#table(data_long_post$site, data_long_post$distortion_correction)
#table(data_long_post$site, data_long_post$eyes)
