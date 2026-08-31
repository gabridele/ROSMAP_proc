# Visualize diagnosis trajectories over follow-up as swimmer plots.
# Produces static and interactive panels stratified by baseline diagnosis while
# preserving visit-level diagnosis changes and a shared time axis.

# Shared input/output path configuration. See README.md for environment variables.
.paths_file <- if (file.exists(file.path("analysis", "paths.R"))) {
  file.path("analysis", "paths.R")
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
  require_file(demos("demos_conn.csv"))
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
  filename = output("figures", "swimmer_plot_shared_x_axis.png"),
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
  filename = output("figures", "swimmer_plot_shared_x_axis.pdf"),
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
# Interactive Plotly version
#
# Use one Plotly figure with explicit axis domains instead of
# nested subplot() calls.
#
# Each panel contains at most:
#   - 3 trajectory traces (one per final diagnosis)
#   - 3 visit-point traces (one per visit diagnosis)
#
# Individual subjects within a trajectory trace are separated by NA.
# This is much faster than creating one Plotly trace per subject.
# ===============================================================


# ---------------------------------------------------------------
# Convert multiple subject trajectories into one Plotly trace
# ---------------------------------------------------------------

make_line_vectors <- function(data) {

  subject_ids <- unique(data$sub_id)

  x <- unlist(
    lapply(
      subject_ids,
      function(id) {
        c(
          data$years_from_baseline[data$sub_id == id],
          NA_real_
        )
      }
    ),
    use.names = FALSE
  )

  y <- unlist(
    lapply(
      subject_ids,
      function(id) {
        c(
          data$subject_y[data$sub_id == id],
          NA_real_
        )
      }
    ),
    use.names = FALSE
  )

  text <- unlist(
    lapply(
      subject_ids,
      function(id) {
        c(
          data$tooltip[data$sub_id == id],
          NA_character_
        )
      }
    ),
    use.names = FALSE
  )

  list(
    x = x,
    y = y,
    text = text
  )
}


# ---------------------------------------------------------------
# Add one swimmer panel to an existing Plotly object
# ---------------------------------------------------------------

add_swimmer_panel <- function(
    p,
    panel_data,
    xaxis_ref,
    yaxis_ref
) {

  # -------------------------------------------------------------
  # Subject trajectory lines
  # -------------------------------------------------------------

  for (diagnosis in c("NCI", "MCI", "AD")) {

    line_data <- panel_data %>%
      filter(final_dx == diagnosis) %>%
      arrange(
        sub_id,
        years_from_baseline,
        age_scandate
      )

    if (nrow(line_data) > 0) {

      line_vectors <- make_line_vectors(line_data)

      p <- p %>%
        add_trace(
          x = line_vectors$x,
          y = line_vectors$y,
          text = line_vectors$text,
          type = "scatter",
          mode = "lines",

          xaxis = xaxis_ref,
          yaxis = yaxis_ref,

          line = list(
            color = unname(
              diagnosis_colours[diagnosis]
            ),
            width = 0.7
          ),

          opacity = 0.60,

          hovertemplate = paste0(
            "%{text}",
            "<extra></extra>"
          ),

          connectgaps = FALSE,

          name = diagnosis,
          legendgroup = diagnosis,
          showlegend = FALSE
        )
    }
  }


  # -------------------------------------------------------------
  # Visit points
  # -------------------------------------------------------------

  for (diagnosis in c("NCI", "MCI", "AD")) {

    point_data <- panel_data %>%
      filter(dcfdx == diagnosis)

    if (nrow(point_data) > 0) {

      p <- p %>%
        add_trace(
          data = point_data,

          x = ~years_from_baseline,
          y = ~subject_y,
          text = ~tooltip,

          type = "scatter",
          mode = "markers",

          xaxis = xaxis_ref,
          yaxis = yaxis_ref,

          marker = list(
            color = unname(
              diagnosis_colours[diagnosis]
            ),
            size = 4,
            opacity = 0.85,
            line = list(
              width = 0
            )
          ),

          hovertemplate = paste0(
            "%{text}",
            "<extra></extra>"
          ),

          name = diagnosis,
          legendgroup = diagnosis,
          showlegend = FALSE
        )
    }
  }

  p
}


# ===============================================================
# Right-hand panel heights
# ===============================================================

right_total <- max(
  n_mci + n_ad,
  1
)

mci_fraction <- n_mci / right_total
ad_fraction  <- n_ad / right_total


# Do not allow either panel to become too small

mci_fraction <- max(
  mci_fraction,
  0.25
)

ad_fraction <- max(
  ad_fraction,
  0.25
)


# Re-normalise

fraction_total <- (
  mci_fraction +
    ad_fraction
)

mci_fraction <- (
  mci_fraction /
    fraction_total
)

ad_fraction <- (
  ad_fraction /
    fraction_total
)


# Small vertical gap between MCI and AD

right_gap <- 0.06

available_right_height <- (
  1 - right_gap
)

ad_height <- (
  available_right_height *
    ad_fraction
)

ad_domain <- c(
  0,
  ad_height
)

mci_domain <- c(
  ad_height + right_gap,
  1
)


# ===============================================================
# Y-axis ranges
# ===============================================================

nci_y_range <- c(
  max(n_nci, 1) + 0.5,
  0.5
)

mci_y_range <- c(
  max(n_mci, 1) + 0.5,
  0.5
)

ad_y_range <- c(
  max(n_ad, 1) + 0.5,
  0.5
)


# ===============================================================
# Create one Plotly figure
# ===============================================================

p_interactive <- plot_ly(
  height = 950
)


# ---------------------------------------------------------------
# Add three dummy traces for one shared legend
# ---------------------------------------------------------------

for (diagnosis in c("NCI", "MCI", "AD")) {

  p_interactive <- p_interactive %>%
    add_markers(
      x = common_x_limits[1],
      y = 1,

      marker = list(
        color = unname(
          diagnosis_colours[diagnosis]
        ),
        size = 7
      ),

      name = diagnosis,
      legendgroup = diagnosis,

      visible = "legendonly",
      hoverinfo = "skip",
      showlegend = TRUE
    )
}


# ---------------------------------------------------------------
# NCI panel
#
# x / y
# ---------------------------------------------------------------

p_interactive <- add_swimmer_panel(
  p = p_interactive,
  panel_data = data_nci,
  xaxis_ref = "x",
  yaxis_ref = "y"
)


# ---------------------------------------------------------------
# MCI panel
#
# x2 / y2
# ---------------------------------------------------------------

p_interactive <- add_swimmer_panel(
  p = p_interactive,
  panel_data = data_mci,
  xaxis_ref = "x2",
  yaxis_ref = "y2"
)


# ---------------------------------------------------------------
# AD panel
#
# x3 / y3
# ---------------------------------------------------------------

p_interactive <- add_swimmer_panel(
  p = p_interactive,
  panel_data = data_ad,
  xaxis_ref = "x3",
  yaxis_ref = "y3"
)


# ===============================================================
# Explicit panel layout
# ===============================================================

p_interactive <- p_interactive %>%
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


    # -----------------------------------------------------------
    # NCI — full-height left panel
    # -----------------------------------------------------------

    xaxis = list(
      domain = c(
        0,
        0.55
      ),

      range = common_x_limits,

      title = list(
        text = "Years from Baseline"
      ),

      showgrid = TRUE,
      zeroline = FALSE
    ),

    yaxis = list(
      domain = c(
        0,
        1
      ),

      range = nci_y_range,

      title = "",

      showticklabels = FALSE,
      ticks = "",
      showgrid = FALSE,
      zeroline = FALSE
    ),


    # -----------------------------------------------------------
    # MCI — upper-right
    # -----------------------------------------------------------

    xaxis2 = list(
      domain = c(
        0.60,
        1
      ),

      range = common_x_limits,

      anchor = "y2",

      showticklabels = FALSE,

      showgrid = TRUE,
      zeroline = FALSE
    ),

    yaxis2 = list(
      domain = mci_domain,

      range = mci_y_range,

      anchor = "x2",

      title = "",

      showticklabels = FALSE,
      ticks = "",
      showgrid = FALSE,
      zeroline = FALSE
    ),


    # -----------------------------------------------------------
    # AD — lower-right
    # -----------------------------------------------------------

    xaxis3 = list(
      domain = c(
        0.60,
        1
      ),

      range = common_x_limits,

      anchor = "y3",

      title = list(
        text = "Years from Baseline"
      ),

      showgrid = TRUE,
      zeroline = FALSE
    ),

    yaxis3 = list(
      domain = ad_domain,

      range = ad_y_range,

      anchor = "x3",

      title = "",

      showticklabels = FALSE,
      ticks = "",
      showgrid = FALSE,
      zeroline = FALSE
    ),


    # -----------------------------------------------------------
    # Panel titles
    # -----------------------------------------------------------

    annotations = list(

      list(
        text = paste0(
          "<b>NCI at Baseline</b> (n = ",
          n_nci,
          ")"
        ),

        x = 0.275,
        y = 1.025,

        xref = "paper",
        yref = "paper",

        showarrow = FALSE,

        xanchor = "center",
        yanchor = "bottom",

        font = list(
          size = 14
        )
      ),

      list(
        text = paste0(
          "<b>MCI at Baseline</b> (n = ",
          n_mci,
          ")"
        ),

        x = 0.80,
        y = 1.025,

        xref = "paper",
        yref = "paper",

        showarrow = FALSE,

        xanchor = "center",
        yanchor = "bottom",

        font = list(
          size = 14
        )
      ),

      list(
        text = paste0(
          "<b>AD at Baseline</b> (n = ",
          n_ad,
          ")"
        ),

        x = 0.80,
        y = ad_domain[2] + 0.012,

        xref = "paper",
        yref = "paper",

        showarrow = FALSE,

        xanchor = "center",
        yanchor = "bottom",

        font = list(
          size = 14
        )
      )
    ),


    # -----------------------------------------------------------
    # General layout
    # -----------------------------------------------------------

    hovermode = "closest",

    legend = list(
      orientation = "h",

      x = 0.5,
      xanchor = "center",

      y = -0.09,

      title = list(
        text = "Diagnosis"
      )
    ),

    margin = list(
      l = 25,
      r = 25,
      t = 120,
      b = 100
    )
  ) %>%

  config(
    displaylogo = FALSE,
    responsive = TRUE,

    modeBarButtonsToRemove = c(
      "lasso2d",
      "select2d"
    )
  )


# Display
p_interactive


# ===============================================================
# Save interactive HTML
# ===============================================================

saveWidget(
  widget = p_interactive,

  file = output(
    "figures",
    "interactive_swimmer_plot_shared_x_axis.html"
  ),

  selfcontained = TRUE
)