library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

demos_withinconn <- read.csv("demos_withinconn.csv")
colnames(demos_withinconn)[3] <- "sub_ses"
# filter columns that are networks

# target network columns
target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

# predictors
predictors <- c("mean_FD", "msex", "site", "age_scandate", "distortion_correction", "eyes")

# ensure categorical variables are factors
demos_withinconn$mean_FD <- as.numeric(demos_withinconn$mean_FD)
demos_withinconn$msex <- factor(demos_withinconn$msex)
demos_withinconn$site <- factor(demos_withinconn$site)
demos_withinconn$age_scandate <- as.numeric(demos_withinconn$age_scandate)
demos_withinconn$distortion_correction <- factor(demos_withinconn$distortion_correction)
demos_withinconn$eyes <- factor(demos_withinconn$eyes)

####
# --- REGRESSION ANALYSIS ---
####


# container for results
results <- list()

# Loop over each network column (e.g., Vis, SomMot, etc.)
for (col in target_cols) {
  
  # # --- SAFETY CHECK ---
  # # Skip this network if:
  # # 1) all values are NA, OR
  # # 2) there is no variance (sd = 0 → constant values → regression would fail)
  # if (all(is.na(rosmap_df[[col]])) ||
  #     sd(rosmap_df[[col]], na.rm = TRUE) == 0) {
  #   next  # skip to next network
  # }
  
  # --- FIT LINEAR MODEL ---
  # Dynamically builds a formula like:
  #   Vis ~ mean_FD + msex + site
  #   SomMot ~ mean_FD + msex + site
  # depending on `col`
  model <- lm(
    reformulate(predictors, response = col),
    data = demos_withinconn
  )
  
  # --- EXTRACT EFFECT OF INTEREST ---
  # Get the regression coefficient (beta) for mean_FD
  beta <- coef(model)["mean_FD"]
  
  # Get the p-value associated with mean_FD
  pval <- summary(model)$coefficients["mean_FD", "Pr(>|t|)"]
  
  # --- STORE RESULTS ---
  # Save results for this network as a small data frame
  results[[col]] <- data.frame(
    network = col,     # network name
    beta_FD = beta,    # effect size of mean_FD
    p_FD    = pval     # p-value of mean_FD
  )
}

# --- COMBINE ALL NETWORK RESULTS ---
# Convert list of data frames into one big data frame
results_df <- do.call(rbind, results)

# --- MULTIPLE COMPARISON CORRECTION ---
# Adjust p-values across networks using FDR (Benjamini-Hochberg)
results_df$q_FD <- p.adjust(results_df$p_FD, method = "fdr")

# --- ADD SIGNIFICANCE LABELS ---
# Convert q-values into stars for plotting/interpretation:
# ***  q < 0.001
# **   q < 0.01
# *    q < 0.05
# ""   not significant
results_df$q_FD_star <- cut(
  results_df$q_FD,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)

print(results_df)
print(model)
# exclude those rows that have mean_FD > 0.3
demos_withinconn_filtered <- demos_withinconn %>%
  filter(mean_FD <= 0.25)

results_filtered <- list()

for (col in target_cols) {
  
  model <- lm(
    reformulate(predictors, response = col),
    data = demos_withinconn_filtered
  )
  
  visreg::visreg(model, "mean_FD", gg = TRUE, data = demos_withinconn_filtered)
  # store plots in list
  results_filtered[[col]] <- visreg::visreg(model, "mean_FD", gg = TRUE, data = demos_withinconn_filtered)
  beta <- coef(model)["mean_FD"]
  pval <- summary(model)$coefficients["mean_FD", "Pr(>|t|)"]
  


  results_filtered[[col]] <- data.frame(
    network = col,
    beta_FD = beta,
    p_FD    = pval
  )
}

results_filtered_df <- do.call(rbind, results_filtered)
print(results_filtered_df)

results_filtered_df$q_FD <- p.adjust(results_filtered_df$p_FD, method = "fdr")
results_filtered_df$q_FD_star <- cut(
  results_filtered_df$q_FD,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)
print(results_filtered_df)

results_df$network <- factor(results_df$network, levels = target_cols)
results_filtered_df$network <- factor(results_filtered_df$network, levels = target_cols)

# plot 1: all subjects
p1 <- ggplot(results_df, aes(x = network, y = beta_FD)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = round(beta_FD, 3)),
    vjust = ifelse(results_df$beta_FD >= 0, -0.4, 1.2),
    size = 4
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.10, 0.20)) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Beta estimates for mean_FD (all subjects)",
    x = "Network",
    y = "Beta for mean_FD"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

p2 <- ggplot(results_filtered_df, aes(x = network, y = beta_FD)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = round(beta_FD, 3)),
    vjust = ifelse(results_filtered_df$beta_FD >= 0, -0.4, 1.2),
    size = 4
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.10, 0.20)) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Beta estimates for mean_FD (mean_FD <= 0.25)",
    x = "Network",
    y = "Beta for mean_FD"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

summary(model)
alias(model)

# plot mean FD for each network and regression line. all in one figure, multiple plots
library(ggplot2)
library(ggpubr)
# reshape data to long format for ggplot
demos_long <- demos_withinconn %>%
  pivot_longer(cols = target_cols, names_to = "network", values_to = "within_conn")
ggplot(demos_long, aes(x = mean_FD, y = within_conn)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~ network, scales = "free_y") +
  theme_minimal() +
  labs(title = "Within-Network Connectivity vs. Mean FD",
       x = "Mean FD",
       y = "Within-Network Connectivity")
# add measure of correlation
ggplot(demos_long, aes(x = mean_FD, y = within_conn)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ network, scales = "free_y") +
  theme_minimal() +
  labs(title = "Within-Network Connectivity vs. Mean FD",
       x = "Mean FD",
       y = "Within-Network Connectivity") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3)

cor.test(demos_withinconn$mean_FD, demos_withinconn$Cont, method = "pearson")

library(tidytable)
cor_labels <- demos_long %>%
  group_by(network) %>%
  summarise(
    test = list(cor.test(mean_FD, within_conn, method = "pearson")),
    .groups = "drop"
  ) %>%
  mutate(
    r = map_dbl(test, ~ .x$estimate),
    ci_low = map_dbl(test, ~ .x$conf.int[1]),
    ci_high = map_dbl(test, ~ .x$conf.int[2]),
    p = map_dbl(test, ~ .x$p.value),
    label = sprintf("r = %.2f\n95%% CI [%.2f, %.2f]\np = %.3g",
                    r, ci_low, ci_high, p)
  )

ggplot(demos_long, aes(x = mean_FD, y = within_conn)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_text(
    data = cor_labels,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.1,
    size = 3
  ) +
  facet_wrap(~ network) +
  theme_minimal(base_size = 13) +
  labs(
    title = "mean_FD vs network connectivity (pre-filtering)",
    x = "mean_FD",
    y = "Connectivity"
  )

cor_labels_filt <- demos_long %>%
  group_by(network) %>%
  filter(mean_FD <= 0.25) %>%
  summarise(
    test = list(cor.test(mean_FD, within_conn, method = "pearson")),
    .groups = "drop"
  ) %>%
  mutate(
    r = map_dbl(test, ~ .x$estimate),
    ci_low = map_dbl(test, ~ .x$conf.int[1]),
    ci_high = map_dbl(test, ~ .x$conf.int[2]),
    p = map_dbl(test, ~ .x$p.value),
    label = sprintf("r = %.2f\n95%% CI [%.2f, %.2f]\np = %.3g",
                    r, ci_low, ci_high, p)
  )

ggplot(demos_long %>% filter(mean_FD <= 0.25), aes(x = mean_FD, y = within_conn)) +
  geom_point(alpha = 0.25, size = 1.5) +
  geom_smooth(method = "lm", formula = within_conn ~ mean_FD + msex + site + age_scandate + distortion_correction + eyes, se = TRUE) +
  geom_text(
    data = cor_labels_filt,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.1,
    size = 3
  ) +
  facet_wrap(~ network) +
  theme_minimal(base_size = 13) +
  labs(
    title = "mean_FD vs network connectivity (filtered: mean_FD <= 0.25)",
    x = "mean_FD",
    y = "Connectivity"
  )
library(visreg)

model_long <- lm(within_conn ~ mean_FD + msex + site + age_scandate + distortion_correction + eyes + network, data = demos_long)
print(demos_long)
visreg::visreg(model_long, "mean_FD", by = "network", gg = TRUE, data = demos_long)


model_long_filt <- lm(within_conn ~ mean_FD + msex + site + age_scandate + distortion_correction + eyes + network, data = demos_long %>% filter(mean_FD <= 0.25))
visreg::visreg(model_long_filt, "mean_FD", by = "network", gg = TRUE, data = demos_long %>% filter(mean_FD <= 0.25))

print(summary(model_long))
print(summary(model_long_filt))


# ============================================================
# 0. Setup
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# ============================================================
# 1. Load and prepare data
# ============================================================

demos_withinconn <- read.csv("demos_withinconn.csv")
colnames(demos_withinconn)[3] <- "sub_ses"

target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

predictors <- c(
  "mean_FD",
  "msex",
  "site",
  "age_scandate",
  "distortion_correction",
  "eyes"
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
    mutate(network = factor(network, levels = target_cols))
}

get_mode <- function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
}

# ============================================================
# 3. Main plotting function
# ============================================================

plot_fd_models <- function(data, title) {
  
  data_long <- make_long(data)
  
  # ----------------------------------------------------------
  # 3.1 Run adjusted and unadjusted models per network
  # ----------------------------------------------------------
  
  model_results <- map_dfr(target_cols, function(net) {
    
    df_net <- data_long %>%
      filter(network == net)
    
    model_adj <- lm(
      within_conn ~ mean_FD + msex + site + age_scandate +
        distortion_correction + eyes,
      data = df_net
    )
    
    model_raw <- lm(
      within_conn ~ mean_FD,
      data = df_net
    )
    
    tibble(
      network = net,
      beta_adjusted = coef(model_adj)["mean_FD"],
      p_adjusted = summary(model_adj)$coefficients["mean_FD", "Pr(>|t|)"],
      beta_unadjusted = coef(model_raw)["mean_FD"],
      p_unadjusted = summary(model_raw)$coefficients["mean_FD", "Pr(>|t|)"]
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
        "adj β = %.3f %s\nraw β = %.3f\nq = %.3g",
        beta_adjusted,
        sig_adjusted,
        beta_unadjusted,
        q_adjusted
      ),
      network = factor(network, levels = target_cols)
    )
  
  # ----------------------------------------------------------
  # 3.2 Create adjusted prediction lines + confidence intervals
  # ----------------------------------------------------------
  
  pred_adjusted <- map_dfr(target_cols, function(net) {
    
    df_net <- data_long %>%
      filter(network == net)
    
    model_adj <- lm(
      within_conn ~ mean_FD + msex + site + age_scandate +
        distortion_correction + eyes,
      data = df_net
    )
    
    newdat <- data.frame(
      mean_FD = seq(
        min(df_net$mean_FD, na.rm = TRUE),
        max(df_net$mean_FD, na.rm = TRUE),
        length.out = 100
      ),
      msex = get_mode(df_net$msex),
      site = get_mode(df_net$site),
      age_scandate = mean(df_net$age_scandate, na.rm = TRUE),
      distortion_correction = get_mode(df_net$distortion_correction),
      eyes = get_mode(df_net$eyes)
    )
    
    newdat$msex <- factor(newdat$msex, levels = levels(df_net$msex))
    newdat$site <- factor(newdat$site, levels = levels(df_net$site))
    newdat$distortion_correction <- factor(
      newdat$distortion_correction,
      levels = levels(df_net$distortion_correction)
    )
    newdat$eyes <- factor(newdat$eyes, levels = levels(df_net$eyes))
    
    preds <- predict(
      model_adj,
      newdata = newdat,
      interval = "confidence"
    )
    
    bind_cols(newdat, as.data.frame(preds)) %>%
      mutate(network = net)
  }) %>%
    mutate(network = factor(network, levels = target_cols))
  
  # ----------------------------------------------------------
  # 3.3 Create unadjusted prediction lines
  # ----------------------------------------------------------
  
  pred_unadjusted <- map_dfr(target_cols, function(net) {
    
    df_net <- data_long %>%
      filter(network == net)
    
    model_raw <- lm(
      within_conn ~ mean_FD,
      data = df_net
    )
    
    newdat <- data.frame(
      mean_FD = seq(
        min(df_net$mean_FD, na.rm = TRUE),
        max(df_net$mean_FD, na.rm = TRUE),
        length.out = 100
      )
    )
    
    preds <- predict(
      model_raw,
      newdata = newdat,
      interval = "confidence"
    )
    
    bind_cols(newdat, as.data.frame(preds)) %>%
      mutate(network = net)
  }) %>%
    mutate(network = factor(network, levels = target_cols))
  
  # ----------------------------------------------------------
  # 3.4 Label positions
  # ----------------------------------------------------------
  
  label_pos <- data_long %>%
    group_by(network) %>%
    summarise(
      x = max(mean_FD, na.rm = TRUE),
      y = max(within_conn, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(model_results, by = "network")
  
  # ----------------------------------------------------------
  # 3.5 Plot
  # ----------------------------------------------------------
  
  ggplot(data_long, aes(x = mean_FD, y = within_conn)) +
    geom_point(alpha = 0.25, size = 1.3) +
    
    # unadjusted/raw slope
    geom_line(
      data = pred_unadjusted,
      aes(x = mean_FD, y = fit, linetype = "Unadjusted"),
      inherit.aes = FALSE,
      linewidth = 0.8
    ) +
    
    # adjusted confidence interval
    geom_ribbon(
      data = pred_adjusted,
      aes(x = mean_FD, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      alpha = 0.20
    ) +
    
    # adjusted slope
    geom_line(
      data = pred_adjusted,
      aes(x = mean_FD, y = fit, linetype = "Adjusted"),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    
    # beta labels
    geom_text(
      data = label_pos,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 3
    ) +
    
    facet_wrap(~ network, scales = "free_y") +
    scale_linetype_manual(
      name = "Model",
      values = c(
        "Adjusted" = "solid",
        "Unadjusted" = "dashed"
      )
    ) +
    theme_minimal(base_size = 13) +
    labs(
      title = title,
      x = "mean_FD",
      y = "Within-network connectivity"
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

# ============================================================
# 4. Create filtered dataset
# ============================================================

demos_filtered <- demos_withinconn %>%
  filter(mean_FD <= fd_threshold)

# ============================================================
# 5. Generate plots
# ============================================================

p_pre <- plot_fd_models(
  demos_withinconn,
  "FD vs connectivity: pre-filtering"
)

p_post <- plot_fd_models(
  demos_filtered,
  paste0("FD vs connectivity: post-filtering, mean_FD <= ", fd_threshold)
)

print(p_pre)
print(p_post)
