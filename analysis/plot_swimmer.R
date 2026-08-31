# Visualize diagnosis trajectories over follow-up as swimmer plots.
# Produces static and interactive panels stratified by baseline diagnosis while
# preserving visit-level diagnosis changes and a shared time axis.

# Shared input/output path configuration. See README.md for environment variables.
.paths_file <- if (file.exists(file.path("analysis", "april26", "paths.R"))) {
  file.path("analysis", "april26", "paths.R")
} else {
  "paths.R"
}
source(.paths_file)
rm(.paths_file)

# ===============================================================
# Libraries
# ===============================================================

library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(plotly)
library(htmlwidgets)
library(grid)


# ===============================================================
# Load data
# ===============================================================

demos_connectivity <- read_csv(
  ap26_require_file(ap26_demos("demos_conn_2807.csv"))
)


# ===============================================================
# Prepare swimmer-plot data
# ===============================================================

# Use the original wide dataset rather than data_long.
# data_long contains one repeated row per connectivity network.

data_swim <- demos_connectivity %>%
  filter(
    !is.na(sub_id),
    !is.na(dcfdx),
    !is.na(age_scandate),
    !is.na(years_from_baseline)
  ) %>%

  # Retain one row per subject visit
  distinct(
    sub_id,
    age_scandate,
    years_from_baseline,
    dcfdx,
    .keep_all = TRUE
  ) %>%

  arrange(
    sub_id,
    years_from_baseline,
    age_scandate
  ) %>%

  group_by(sub_id) %>%
  mutate(
    baseline_dx = first(dcfdx),
    final_dx = last(dcfdx)
  ) %>%
  ungroup() %>%

  filter(
    baseline_dx %in% c("NCI", "MCI", "AD"),
    final_dx %in% c("NCI", "MCI", "AD"),
    dcfdx %in% c("NCI", "MCI", "AD")
  ) %>%

  mutate(
    baseline_dx = factor(
      baseline_dx,
      levels = c("NCI", "MCI", "AD")
    ),

    dcfdx = factor(
      dcfdx,
      levels = c("NCI", "MCI", "AD")
    ),

    final_dx = factor(
      final_dx,
      levels = c("NCI", "MCI", "AD")
    ),

    tooltip = paste0(
      "Subject: ", sub_id,
      "<br>Visit diagnosis: ", dcfdx,
      "<br>Baseline diagnosis: ", baseline_dx,
      "<br>Final diagnosis: ", final_dx,
      "<br>Years from baseline: ",
      round(years_from_baseline, 2),
      "<br>Scan date/age: ", age_scandate
    )
  )


# ===============================================================
# Diagnosis colours
# ===============================================================

diagnosis_colours <- c(
  "NCI" = "#1f77b4",
  "MCI" = "#ff7f0e",
  "AD"  = "#2ca02c"
)


# ===============================================================
# Prepare each panel
# ===============================================================

prepare_panel_data <- function(data, diagnosis) {

  panel_data <- data %>%
    filter(baseline_dx == diagnosis) %>%
    arrange(
      sub_id,
      years_from_baseline,
      age_scandate
    )

  subject_order <- panel_data %>%
    distinct(sub_id) %>%
    pull(sub_id)

  panel_data %>%
    mutate(
      # Separate subject ordering within each diagnosis panel
      sub_id_y = factor(
        sub_id,
        levels = rev(subject_order)
      ),

      # Numeric equivalent for Plotly
      subject_y = match(sub_id, subject_order)
    )
}


data_nci <- prepare_panel_data(data_swim, "NCI")
data_mci <- prepare_panel_data(data_swim, "MCI")
data_ad  <- prepare_panel_data(data_swim, "AD")


# ===============================================================
# Subject counts
# ===============================================================

n_nci <- n_distinct(data_nci$sub_id)
n_mci <- n_distinct(data_mci$sub_id)
n_ad  <- n_distinct(data_ad$sub_id)

group_counts <- tibble(
  baseline_dx = c("NCI", "MCI", "AD"),
  subjects = c(n_nci, n_mci, n_ad)
)

print(group_counts)


# ===============================================================
# Shared x-axis limits
#
# All three panels use exactly the same time range.
# Their y-axes remain independent.
# ===============================================================

x_range <- range(
  data_swim$years_from_baseline,
  na.rm = TRUE
)

x_padding <- diff(x_range) * 0.03

if (!is.finite(x_padding) || x_padding == 0) {
  x_padding <- 0.1
}

common_x_limits <- c(
  x_range[1] - x_padding,
  x_range[2] + x_padding
)


# ===============================================================
# Static ggplot panel function
# ===============================================================

make_static_panel <- function(
    panel_data,
    panel_title
) {

  ggplot(
    panel_data,
    aes(
      x = years_from_baseline,
      y = sub_id_y,
      group = sub_id
    )
  ) +

    # Very thin subject trajectories
    geom_line(
      aes(color = final_dx),
      linewidth = 0.12,
      alpha = 0.65,
      lineend = "round"
    ) +

    # Very small visit points
    geom_point(
      aes(color = dcfdx),
      size = 0.45,
      alpha = 0.85
    ) +

    scale_color_manual(
      values = diagnosis_colours,
      limits = c("NCI", "MCI", "AD"),
      breaks = c("NCI", "MCI", "AD"),
      drop = FALSE,
      na.translate = FALSE
    ) +

    # Same x-axis across NCI, MCI and AD
    scale_x_continuous(
      limits = common_x_limits,
      expand = c(0, 0)
    ) +

    # No fixed y scale:
    # each panel keeps its own subjects and vertical spacing
    labs(
      title = panel_title,
      x = "Years from Baseline",
      y = NULL,
      color = "Diagnosis"
    ) +

    theme_minimal(base_size = 10) +

    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),

      axis.title.x = element_text(
        size = 9,
        margin = margin(t = 5)
      ),

      axis.text.x = element_text(
        size = 8
      ),

      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),

      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),

      panel.border = element_rect(
        fill = NA,
        colour = "grey70",
        linewidth = 0.3
      ),

      plot.margin = margin(
        t = 6,
        r = 6,
        b = 6,
        l = 6
      )
    )
}


# ===============================================================
# Create static panels
# ===============================================================

p_nci <- make_static_panel(
  panel_data = data_nci,
  panel_title = paste0(
    "NCI at Baseline (n = ",
    n_nci,
    ")"
  )
)

p_mci <- make_static_panel(
  panel_data = data_mci,
  panel_title = paste0(
    "MCI at Baseline (n = ",
    n_mci,
    ")"
  )
)

p_ad <- make_static_panel(
  panel_data = data_ad,
  panel_title = paste0(
    "AD at Baseline (n = ",
    n_ad,
    ")"
  )
)


# ===============================================================
# Static layout
#
# NCI fills the left side.
# MCI is upper-right.
# AD is lower-right.
#
# X-axis: identical across all panels
# Y-axis: independent in each panel
# ===============================================================

panel_design <- c(
  area(t = 1, l = 1, b = 2, r = 1),
  area(t = 1, l = 2, b = 1, r = 2),
  area(t = 2, l = 2, b = 2, r = 2)
)


# Give MCI and AD vertical space in proportion to their sample sizes
right_panel_heights <- c(
  max(n_mci, 1),
  max(n_ad, 1)
)


p_static <- (
  p_nci +
    p_mci +
    p_ad +

    plot_layout(
      design = panel_design,
      widths = c(1.15, 1),
      heights = right_panel_heights,
      guides = "collect"
    ) +

    plot_annotation(
      title = "Diagnosis Trajectories Over Time",

      subtitle = paste0(
        "Lines are coloured by final diagnosis; ",
        "points are coloured by visit diagnosis"
      ),

      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 16
        ),

        plot.subtitle = element_text(
          hjust = 0.5,
          size = 10
        )
      )
    )
) &

  theme(
    legend.position = "bottom",

    legend.title = element_text(
      face = "bold"
    ),

    legend.key.width = unit(
      1.2,
      "cm"
    )
  )


print(p_static)


# ===============================================================
# Save static plot
# ===============================================================

ggsave(
  filename = ap26_output("figures", "swimmer_plot_shared_x_axis.png"),
  plot = p_static,
  width = 18,
  height = 9,
  units = "in",
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)


# Optional vector version
ggsave(
  filename = ap26_output("figures", "swimmer_plot_shared_x_axis.pdf"),
  plot = p_static,
  width = 18,
  height = 9,
  units = "in",
  limitsize = FALSE,
  bg = "white"
)


# ===============================================================
# Native Plotly panel function
#
# Native Plotly avoids the ggplot2 warning about the
# unrecognised "text" aesthetic.
# ===============================================================

make_plotly_panel <- function(
    panel_data,
    panel_title,
    add_legend = FALSE,
    height_px = 950
) {

  panel_n <- n_distinct(panel_data$sub_id)

  panel_y_max <- max(
    panel_n,
    1
  )

  # Height is specified in plot_ly(), not layout()
  p <- plot_ly(
    height = height_px
  )


  # -------------------------------------------------------------
  # Add legend entries once
  # -------------------------------------------------------------

  if (add_legend) {

    for (diagnosis in c("NCI", "MCI", "AD")) {

      p <- p %>%
        add_markers(
          x = common_x_limits[1],
          y = 1,

          name = diagnosis,
          legendgroup = diagnosis,
          visible = "legendonly",

          marker = list(
            color = unname(
              diagnosis_colours[diagnosis]
            ),
            size = 5
          ),

          hoverinfo = "skip",
          showlegend = TRUE
        )
    }
  }


  # -------------------------------------------------------------
  # Add subject trajectory lines
  # -------------------------------------------------------------

  for (diagnosis in c("NCI", "MCI", "AD")) {

    line_data <- panel_data %>%
      filter(final_dx == diagnosis)

    if (nrow(line_data) > 0) {

      p <- p %>%
        add_trace(
          data = line_data,

          x = ~years_from_baseline,
          y = ~subject_y,

          # One separate line per subject
          split = ~sub_id,

          type = "scatter",
          mode = "lines",

          text = ~tooltip,
          hoverinfo = "text",

          line = list(
            color = unname(
              diagnosis_colours[diagnosis]
            ),

            # Very thin interactive line
            width = 0.35
          ),

          opacity = 0.65,

          name = diagnosis,
          legendgroup = diagnosis,
          showlegend = FALSE
        )
    }
  }


  # -------------------------------------------------------------
  # Add visit points
  # -------------------------------------------------------------

  for (diagnosis in c("NCI", "MCI", "AD")) {

    point_data <- panel_data %>%
      filter(dcfdx == diagnosis)

    if (nrow(point_data) > 0) {

      p <- p %>%
        add_markers(
          data = point_data,

          x = ~years_from_baseline,
          y = ~subject_y,

          text = ~tooltip,
          hoverinfo = "text",

          marker = list(
            color = unname(
              diagnosis_colours[diagnosis]
            ),

            # Very small interactive points
            size = 2.5,

            opacity = 0.85,

            line = list(
              width = 0
            )
          ),

          name = diagnosis,
          legendgroup = diagnosis,
          showlegend = FALSE
        )
    }
  }


  # -------------------------------------------------------------
  # Panel layout
  # -------------------------------------------------------------

  p %>%
    layout(
      title = list(
        text = panel_title,
        x = 0.5,
        xanchor = "center",

        font = list(
          size = 14
        )
      ),

      # Same x-axis limits in every panel
      xaxis = list(
        title = "Years from Baseline",
        range = common_x_limits,
        showgrid = TRUE,
        zeroline = FALSE
      ),

      # Independent y-axis based on this panel's subjects
      yaxis = list(
        title = "",
        range = c(
          panel_y_max + 0.5,
          0.5
        ),
        showticklabels = FALSE,
        ticks = "",
        showgrid = FALSE,
        zeroline = FALSE
      ),

      hovermode = "closest",

      margin = list(
        l = 20,
        r = 15,
        t = 55,
        b = 50
      )
    )
}


# ===============================================================
# Interactive panel dimensions
# ===============================================================

interactive_height <- 950

right_total <- max(
  n_mci + n_ad,
  1
)

mci_fraction <- n_mci / right_total
ad_fraction  <- n_ad / right_total


# Prevent either panel from becoming extremely short
mci_fraction <- max(
  mci_fraction,
  0.25
)

ad_fraction <- max(
  ad_fraction,
  0.25
)


# Normalise so the two fractions sum to one
fraction_total <- mci_fraction + ad_fraction

right_heights <- c(
  mci_fraction / fraction_total,
  ad_fraction / fraction_total
)


# ===============================================================
# Create interactive panels
# ===============================================================

plotly_nci <- make_plotly_panel(
  panel_data = data_nci,

  panel_title = paste0(
    "NCI at Baseline (n = ",
    n_nci,
    ")"
  ),

  add_legend = TRUE,
  height_px = interactive_height
)


plotly_mci <- make_plotly_panel(
  panel_data = data_mci,

  panel_title = paste0(
    "MCI at Baseline (n = ",
    n_mci,
    ")"
  ),

  add_legend = FALSE,

  height_px = round(
    interactive_height * right_heights[1]
  )
)


plotly_ad <- make_plotly_panel(
  panel_data = data_ad,

  panel_title = paste0(
    "AD at Baseline (n = ",
    n_ad,
    ")"
  ),

  add_legend = FALSE,

  height_px = round(
    interactive_height * right_heights[2]
  )
)


# ===============================================================
# Stack MCI and AD on the right
#
# shareX = TRUE
# shareY = FALSE
# ===============================================================

plotly_right <- subplot(
  plotly_mci,
  plotly_ad,

  nrows = 2,

  heights = right_heights,

  shareX = TRUE,
  shareY = FALSE,

  titleX = TRUE,
  titleY = FALSE,

  margin = 0.04
)


# ===============================================================
# Combine NCI with the right-hand panels
#
# The x limits are explicitly identical in every panel.
# The y axes remain independent.
# ===============================================================

p_interactive <- subplot(
  plotly_nci,
  plotly_right,

  nrows = 1,

  widths = c(
    0.56,
    0.44
  ),

  shareX = TRUE,
  shareY = FALSE,

  titleX = TRUE,
  titleY = FALSE,

  margin = 0.04
) %>%

  layout(
    title = list(
      text = paste0(
        "Interactive Diagnosis Trajectories Over Time",
        "<br>",
        "<sup>",
        "Lines: final diagnosis; ",
        "points: visit diagnosis",
        "</sup>"
      ),

      x = 0.5,
      xanchor = "center"
    ),

    hovermode = "closest",

    legend = list(
      orientation = "h",
      x = 0.5,
      xanchor = "center",
      y = -0.07
    ),

    margin = list(
      l = 30,
      r = 20,
      t = 90,
      b = 90
    )
  ) %>%

  config(
    displaylogo = FALSE,
    responsive = TRUE
  )


p_interactive


# ===============================================================
# Save interactive HTML
# ===============================================================

saveWidget(
  widget = p_interactive,
  file = ap26_output("figures", "interactive_swimmer_plot_shared_x_axis.html"),
  selfcontained = TRUE
)