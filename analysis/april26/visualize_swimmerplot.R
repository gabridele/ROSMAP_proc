library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggeffects)
library(lme4)
library(lmerTest)
library(lubridate)
library(plotly)

# ===============================================================
# Load data and prepare long format
# ===============================================================

demos_connectivity <- read_csv("demos_conn_2505.csv")

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

data_long <- make_long(demos_connectivity)

# swimmer plot divided by diagnosis at baseline on y axis, on x axis is years from baseline. each datapoint is colored by diagnosis. connected by line whose color is the diagnosis at final visit

library(dplyr)
library(ggplot2)
library(forcats)

data_swim <- data_long %>%
  filter(
    !is.na(dcfdx),
    !is.na(age_scandate),
    !is.na(sub_id)
  ) %>%
  arrange(sub_id, age_scandate) %>%
  group_by(sub_id) %>%
  mutate(
    baseline_date = first(age_scandate),
    baseline_dx = first(dcfdx),
    final_dx = last(dcfdx)
    ) %>%
    ungroup() %>%
    filter(
      !is.na(baseline_dx),
      !is.na(final_dx)
    ) %>%
    mutate(
      baseline_dx = factor(baseline_dx, levels = c("NCI", "MCI", "AD")),
      dcfdx = factor(dcfdx, levels = c("NCI", "MCI", "AD")),
      final_dx = factor(final_dx, levels = c("NCI", "MCI", "AD"))
    ) %>%
    arrange(baseline_dx, sub_id) %>%
    mutate(
      sub_id_y = factor(sub_id, levels = unique(sub_id)),
      tooltip = paste0(
        "Subject: ", sub_id,
        "<br>Visit diagnosis: ", dcfdx,
        "<br>Baseline diagnosis: ", baseline_dx,
        "<br>Final diagnosis: ", final_dx,
        "<br>Years from baseline: ", round(years_from_baseline, 2),
        "<br>Scan date: ", age_scandate
      )
    )

ggplot(data_swim, aes(x = years_from_baseline, y = sub_id_y, group = sub_id)) +
  geom_line(aes(color = final_dx), linewidth = 0.8, alpha = 0.7) +
  geom_point(aes(color = dcfdx), size = 2.5, alpha = 0.9) +
  facet_grid(baseline_dx ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(
    values = c(
      "NCI" = "#1f77b4",
      "MCI" = "#ff7f0e",
      "AD"  = "#2ca02c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Years from Baseline",
    y = "Subjects, grouped by baseline diagnosis",
    color = "Diagnosis",
    title = "Swimmer Plot of Diagnosis Trajectories Over Time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing.y = unit(1.5, "lines"),
    strip.text.y = element_text(size = 13, face = "bold"),
    plot.margin = margin(10, 20, 10, 20)
  )

ggsave(
  "swimmer_plot.png",
  width = 10,
  height = 22,
  dpi = 300
)

p <- ggplot(
  data_swim,
  aes(
    x = years_from_baseline,
    y = sub_id_y,
    group = sub_id,
    text = tooltip
  )
) +
  geom_line(
    aes(color = final_dx),
    linewidth = 0.8,
    alpha = 0.7
  ) +
  geom_point(
    aes(color = dcfdx),
    size = 2.5,
    alpha = 0.9
  ) +
  facet_grid(
    baseline_dx ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_color_manual(
    values = c(
      "NCI" = "#1f77b4",
      "MCI" = "#ff7f0e",
      "AD"  = "#2ca02c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Years from Baseline",
    y = "Subjects, grouped by baseline diagnosis",
    color = "Diagnosis",
    title = "Interactive Swimmer Plot of Diagnosis Trajectories"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing.y = unit(1.5, "lines"),
    strip.text.y = element_text(size = 13, face = "bold"),
    plot.margin = margin(10, 20, 10, 20)
  )

p_interactive <- ggplotly(p, tooltip = "text") %>%
  layout(
    hovermode = "closest",
    legend = list(
      orientation = "h",
      x = 0.3,
      y = -0.1
    )
  )

p_interactive
library(htmlwidgets)

saveWidget(
  p_interactive,
  file = "interactive_swimmer_plot.html",
  selfcontained = TRUE
)
