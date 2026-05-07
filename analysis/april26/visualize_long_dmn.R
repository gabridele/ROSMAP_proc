library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggeffects)
library(lme4)
library(lmerTest)


demos_withinconn <- read.csv("demos_withinconn.csv")

demos_withinconn <- read.csv("demos_withinconn.csv")

demos_withinconn <- demos_withinconn %>%
  mutate(
    ses_num = as.numeric(str_extract(ses, "\\d+"))
  ) %>%
  mutate(sub = factor(sub))

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
data_long <- make_long(demos_withinconn)

network_to_plot <- "Default"

df_one_net <- data_long %>%
  filter(network == network_to_plot)

# ============================================================

df_segments <- df_one_net %>%
  arrange(sub, ses_num) %>%
  group_by(sub) %>%
  mutate(
    x_start = ses_num,
    x_end = lead(ses_num),
    y_start = within_conn,
    y_end = lead(within_conn),
    direction = case_when(
      y_end > y_start ~ "up",
      y_end < y_start ~ "down",
      TRUE ~ "flat"
    )
  ) %>%
  filter(!is.na(x_end)) %>%
  ungroup()

ggplot() +
  geom_segment(
    data = df_segments,
    aes(
      x = x_start,
      xend = x_end,
      y = y_start,
      yend = y_end,
      color = direction
    ),
    alpha = 0.6,
    linewidth = 0.6
  ) +
  geom_point(
    data = df_one_net,
    aes(x = ses_num, y = within_conn),
    alpha = 0.5,
    size = 1.5
  ) +
  scale_color_manual(values = c(
    "up" = "#2C7BB6",
    "down" = "#D7191C",
    "flat" = "grey50"
  )) +
  scale_x_continuous(breaks = sort(unique(df_one_net$ses_num))) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste(network_to_plot, "connectivity across sessions"),
    x = "Session number",
    y = "DMN",
    color = "Step direction"
  )

ggplot(df_one_net, aes(x = ses_num, y = within_conn)) +
  geom_line(aes(group = sub), alpha = 0.25) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_x_continuous(breaks = sort(unique(df_one_net$ses_num))) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste(network_to_plot, "connectivity across sessions"),
    x = "Session number",
    y = "DMN"
  )

# -------
## split into 4 groups for better visualization
set.seed(123)

sub_groups <- df_one_net %>%
  distinct(sub) %>%
  mutate(group = sample(rep(1:4, length.out = nrow(.))))
df_one_net <- df_one_net %>%
  left_join(sub_groups, by = "sub")

df_segments <- df_segments %>%
  left_join(sub_groups, by = "sub")

ggplot() +
  geom_segment(
    data = df_segments,
    aes(
      x = x_start,
      xend = x_end,
      y = y_start,
      yend = y_end,
      color = direction
    ),
    alpha = 0.6,
    linewidth = 0.6
  ) +
  
  geom_point(
    data = df_one_net,
    aes(x = ses_num, y = within_conn),
    alpha = 0.5,
    size = 1.5
  ) +
  
  facet_wrap(~ group) +
  
  scale_color_manual(values = c(
    "up" = "#2C7BB6",
    "down" = "#D7191C",
    "flat" = "grey50"
  )) +
  
  theme_minimal(base_size = 13) +
  labs(
    title = paste(network_to_plot, "connectivity across sessions (split view)"),
    x = "Session number",
    y = "DMN",
    color = "Step direction"
  )

####### =======================================================
# same but after filtering for mean_FD < 0.25
####### =======================================================



data_long_post <- data_long %>%
  filter(mean_FD < 0.25)

df_one_net_postFD <- data_long_post %>%
  filter(network == network_to_plot)

df_segments <- df_one_net_postFD %>%
  arrange(sub, ses_num) %>%
  group_by(sub) %>%
  mutate(
    x_start = ses_num,
    x_end = lead(ses_num),
    y_start = within_conn,
    y_end = lead(within_conn),
    direction = case_when(
      y_end > y_start ~ "up",
      y_end < y_start ~ "down",
      TRUE ~ "flat"
    )
  ) %>%
  filter(!is.na(x_end)) %>%
  ungroup()

post_fd <- ggplot() +
  geom_segment(
    data = df_segments,
    aes(
      x = x_start,
      xend = x_end,
      y = y_start,
      yend = y_end,
      color = direction
    ),
    alpha = 0.6,
    linewidth = 0.6
  ) +
  geom_point(
    data = df_one_net_postFD,
    aes(x = ses_num, y = within_conn),
    alpha = 0.5,
    size = 1.5
  ) +
  scale_color_manual(values = c(
    "up" = "#2C7BB6",
    "down" = "#D7191C",
    "flat" = "grey50"
  )) +
  scale_x_continuous(breaks = sort(unique(df_one_net_postFD$ses_num))) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste(network_to_plot, "connectivity across sessions"),
    x = "Session number",
    y = "DMN",
    color = "Step direction"
  )

# save
ggsave(
  "dmn_withinconn_longitudinal_postFD.png",
  plot = post_fd,
  width = 12,
  height = 8,
  dpi = 300
)

# -------
## split into 4 groups for better visualization
set.seed(123)

sub_groups <- df_one_net_postFD %>%
  distinct(sub) %>%
  mutate(group = sample(rep(1:4, length.out = nrow(.))))
df_one_net_postFD <- df_one_net_postFD %>%
  left_join(sub_groups, by = "sub")

df_segments <- df_segments %>%
  left_join(sub_groups, by = "sub") %>%
  filter(sub %in% df_one_net_postFD$sub) # keep only segments for subjects that passed the FD filter

post_fd_split <- ggplot() +
  geom_segment(
    data = df_segments,
    aes(
      x = x_start,
      xend = x_end,
      y = y_start,
      yend = y_end,
      color = direction
    ),
    alpha = 0.6,
    linewidth = 0.6
  ) +
  
  geom_point(
    data = df_one_net_postFD,
    aes(x = ses_num, y = within_conn),
    alpha = 0.5,
    size = 1.5
  ) +
  
  facet_wrap(~ group) +
  
  scale_color_manual(values = c(
    "up" = "#2C7BB6",
    "down" = "#D7191C",
    "flat" = "grey50"
  )) +
  
  theme_minimal(base_size = 13) +
  labs(
    title = paste(network_to_plot, "connectivity across sessions (split view)"),
    x = "Session number",
    y = "DMN",
    color = "Step direction"
  )

# save
ggsave(
  "dmn_withinconn_longitudinal_postFD_split.png",
  plot = post_fd_split,
  width = 12,
  height = 8,
  dpi = 300
)
