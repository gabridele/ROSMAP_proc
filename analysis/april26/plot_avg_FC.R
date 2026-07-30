library(reticulate)
library(ggplot2)
library(dplyr)
library(readr)

# Load NumPy
np <- import("numpy", convert = TRUE)

# ------------------------------------------------------------
# 1. Paths and matrix names
# ------------------------------------------------------------

files <- c(
  "/Users/ga0034de/Desktop/BNKBBR_ts/FC/avg/avg_fc_bnkbbr_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/bnk/avg/avg_fc_bnk_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/uc/avg/avg_fc_uc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/mg/avg/avg_fc_mg_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/rirc/avg/avg_fc_rirc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/fc_matrices_456/avg/avg_fc_matrix.npy"
)

matrix_names <- c(
  "BNKBBR",
  "BNK",
  "UC",
  "MG",
  "RIRC",
  "whole_dataset"
)

if (length(files) != length(matrix_names)) {
  stop("The number of files does not match the number of matrix names.")
}

# ------------------------------------------------------------
# 2. Read network labels
# ------------------------------------------------------------

dseg_file <- paste0(
  "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/",
  "atlas-4S456Parcels/atlas-4S456Parcels_dseg.tsv"
)

dseg <- read_tsv(
  dseg_file,
  na = c("", "NA", "N/A", "n/a", "NaN"),
  show_col_types = FALSE
)

required_columns <- c("network_label", "atlas_name")

if (!all(required_columns %in% names(dseg))) {
  stop(
    "The dseg TSV must contain these columns: ",
    paste(required_columns, collapse = ", ")
  )
}

network <- dseg %>%
  transmute(
    network_label = na_if(
      trimws(as.character(network_label)),
      ""
    ),
    atlas_name = na_if(
      trimws(as.character(atlas_name)),
      ""
    ),

    # Use atlas_name when network_label is missing
    plot_label = coalesce(
      network_label,
      atlas_name,
      "Unlabelled"
    )
  ) %>%
  pull(plot_label)

# ------------------------------------------------------------
# 3. Load and validate matrices
# ------------------------------------------------------------

matrices <- lapply(files, function(file) {

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  x <- np$load(
    file,
    allow_pickle = FALSE
  )

  x <- as.matrix(x)

  if (nrow(x) != ncol(x)) {
    stop(file, " does not contain a square matrix.")
  }

  x
})

# Confirm all matrices have identical dimensions
matrix_dimensions <- vapply(
  matrices,
  function(x) paste(dim(x), collapse = "x"),
  character(1)
)

if (length(unique(matrix_dimensions)) != 1) {
  stop(
    "Matrices do not have matching dimensions: ",
    paste(unique(matrix_dimensions), collapse = ", ")
  )
}

number_of_rois <- nrow(matrices[[1]])

if (length(network) != number_of_rois) {
  stop(
    "The number of dseg rows does not match the matrix dimensions. ",
    "Matrix has ",
    number_of_rois,
    " ROIs, but dseg contains ",
    length(network),
    " rows."
  )
}

# ------------------------------------------------------------
# 4. Optionally group ROIs by network
# ------------------------------------------------------------

# TRUE:
#   ROIs are reordered so that members of each network are adjacent.
#
# FALSE:
#   Original atlas/dseg order is preserved.
#
# The same ordering is applied to rows and columns of every matrix.

reorder_by_network <- FALSE

if (reorder_by_network) {

  # Preserve network order based on first appearance in dseg
  network_factor <- factor(
    network,
    levels = unique(network)
  )

  roi_order <- order(network_factor)

} else {

  roi_order <- seq_len(number_of_rois)
}

network_plot_order <- network[roi_order]

# Apply identical ordering to every connectivity matrix
matrices <- lapply(
  matrices,
  function(x) x[roi_order, roi_order, drop = FALSE]
)

# ------------------------------------------------------------
# 5. Find network block locations
# ------------------------------------------------------------

# rle() identifies consecutive blocks of identical labels
network_runs <- rle(network_plot_order)

network_blocks <- data.frame(
  network = network_runs$values,
  block_length = network_runs$lengths,
  stringsAsFactors = FALSE
)

network_blocks$end <- cumsum(
  network_blocks$block_length
)

network_blocks$start <- c(
  1,
  head(network_blocks$end, -1) + 1
)

# Position at which the network name is displayed
network_blocks$midpoint <- (
  network_blocks$start +
    network_blocks$end
) / 2

# Lines are placed between adjacent network blocks
network_boundaries <- head(
  network_blocks$end,
  -1
) + 0.5

# ------------------------------------------------------------
# 6. Convert matrices to long format
# ------------------------------------------------------------

matrix_to_dataframe <- function(x, matrix_name) {

  data.frame(
    row = rep(
      seq_len(nrow(x)),
      times = ncol(x)
    ),
    column = rep(
      seq_len(ncol(x)),
      each = nrow(x)
    ),
    connectivity = as.vector(x),
    matrix = matrix_name
  )
}

plot_data <- do.call(
  rbind,
  Map(
    matrix_to_dataframe,
    matrices,
    matrix_names
  )
)

plot_data$matrix <- factor(
  plot_data$matrix,
  levels = matrix_names
)

# Symmetric colour limits for every matrix
colour_limit <- 1

# ------------------------------------------------------------
# 7. Function for creating a labelled FC plot
# ------------------------------------------------------------

make_fc_plot <- function(
    data,
    title,
    subtitle = NULL,
    facet = FALSE,
    axis_text_size = 7
) {

  p <- ggplot(
    data,
    aes(
      x = column,
      y = row,
      fill = connectivity
    )
  ) +
    geom_tile() +

    # Vertical boundaries between networks
    geom_vline(
      xintercept = network_boundaries,
      linewidth = 0.25,
      colour = "black",
      alpha = 0.7
    ) +

    # Horizontal boundaries between networks
    geom_hline(
      yintercept = network_boundaries,
      linewidth = 0.25,
      colour = "black",
      alpha = 0.7
    ) +

    # Network labels along the top
    scale_x_continuous(
      position = "top",
      breaks = network_blocks$midpoint,
      labels = network_blocks$network,
      expand = c(0, 0)
    ) +

    # Network labels along the left
    scale_y_reverse(
      breaks = network_blocks$midpoint,
      labels = network_blocks$network,
      expand = c(0, 0)
    ) +

    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(
        -colour_limit,
        colour_limit
      ),
      oob = scales::squish,
      na.value = "grey90",
      name = "Connectivity"
    ) +

    coord_fixed(
      clip = "off"
    ) +

    labs(
      title = title,
      #subtitle = subtitle,
      x = NULL,
      y = NULL
    ) +

    theme_minimal(
      base_size = 12
    ) +

    theme(
      panel.grid = element_blank(),

      strip.text = element_text(
        face = "bold"
      ),

      axis.ticks = element_blank(),

      # Top network labels
      axis.text.x.top = element_text(
        size = axis_text_size,
        angle = 45,
        hjust = 0,
        vjust = 0
      ),

      # Left network labels
      axis.text.y.left = element_text(
        size = axis_text_size,
        hjust = 1
      ),

      plot.title = element_text(
        face = "bold",
        size = 14
      ),

      plot.subtitle = element_text(
        size = 10,
        margin = margin(b = 8)
      ),

      plot.margin = margin(
        t = 15,
        r = 15,
        b = 10,
        l = 15
      )
    )

  if (facet) {
    p <- p +
      facet_wrap(
        ~matrix,
        ncol = 2
      )
  }

  p
}

# ------------------------------------------------------------
# 8. Faceted plot containing all matrices
# ------------------------------------------------------------

fc_plot <- make_fc_plot(
  data = plot_data,
  title = "Average Functional Connectivity Matrices",
  facet = TRUE,
  axis_text_size = 6
)

print(fc_plot)

ggsave(
  filename = "all_average_fc_matrices.svg",
  plot = fc_plot,
  width = 14,
  height = 18,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 9. Save each matrix separately
# ------------------------------------------------------------

for (i in seq_along(matrices)) {

  plot_data_i <- matrix_to_dataframe(
    matrices[[i]],
    matrix_names[i]
  )

  single_plot <- make_fc_plot(
    data = plot_data_i,
    title = paste(
      "Average Functional Connectivity Matrix -",
      matrix_names[i]
    ),
    facet = FALSE,
    axis_text_size = 7
  )

  print(single_plot)

  safe_name <- gsub(
    pattern = "[^A-Za-z0-9_-]+",
    replacement = "_",
    x = matrix_names[i]
  )

  ggsave(
    filename = paste0(
      "average_fc_matrix_",
      safe_name,
      ".svg"
    ),
    plot = single_plot,
    width = 9,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# ------------------------------------------------------------
# 10. Individual matrices with summary statistics
# ------------------------------------------------------------

# Set TRUE to exclude diagonal values from statistics
exclude_diagonal <- TRUE

for (i in seq_along(matrices)) {

  current_matrix <- matrices[[i]]

  # Values used for summary statistics
  stats_values <- current_matrix

  if (exclude_diagonal) {
    diag(stats_values) <- NA_real_
  }

  stats_values <- as.numeric(stats_values)

  stats_values <- stats_values[
    is.finite(stats_values)
  ]

  if (length(stats_values) == 0) {
    warning(
      "No finite values available for ",
      matrix_names[i]
    )

    matrix_stats <- c(
      Mean = NA_real_,
      SD = NA_real_,
      Min = NA_real_,
      Max = NA_real_,
      Median = NA_real_
    )

  } else {

    matrix_stats <- c(
      Mean = mean(stats_values),
      SD = sd(stats_values),
      Min = min(stats_values),
      Max = max(stats_values),
      Median = median(stats_values)
    )
  }

  # stats_text <- sprintf(
  #   paste0(
  #     "Mean: %.3f   SD: %.3f   Min: %.3f   ",
  #     "Max: %.3f   Median: %.3f"
  #   ),
  #   matrix_stats["Mean"],
  #   matrix_stats["SD"],
  #   matrix_stats["Min"],
  #   matrix_stats["Max"],
  #   matrix_stats["Median"]
  # )

  plot_data_i <- matrix_to_dataframe(
    current_matrix,
    matrix_names[i]
  )

  single_plot <- make_fc_plot(
    data = plot_data_i,
    title = paste(
      "Average Functional Connectivity Matrix -",
      matrix_names[i]
    ),
    #subtitle = stats_text,
    facet = FALSE,
    axis_text_size = 7
  )

  print(single_plot)

  safe_name <- gsub(
    pattern = "[^A-Za-z0-9_-]+",
    replacement = "_",
    x = matrix_names[i]
  )

  # ggsave(
  #   filename = paste0(
  #     "average_fc_matrix_",
  #     safe_name,
  #     "_with_statistics.svg"
  #   ),
  #   plot = single_plot,
  #   width = 9,
  #   height = 8,
  #   dpi = 300,
  #   bg = "white"
  # )
}

# ------------------------------------------------------------
# Save each matrix as a clean png
# ------------------------------------------------------------

for (i in seq_along(matrices)) {

  plot_data_i <- matrix_to_dataframe(
    matrices[[i]],
    matrix_names[i]
  )

  matrix_plot <- ggplot(
    plot_data_i,
    aes(
      x = column,
      y = row,
      fill = connectivity
    )
  ) +
    geom_raster(
      interpolate = FALSE
    ) +

    # Lines separating network blocks
    geom_vline(
      xintercept = network_boundaries,
      colour = "black",
      linewidth = 0.25
    ) +
    geom_hline(
      yintercept = network_boundaries,
      colour = "black",
      linewidth = 0.25
    ) +

    scale_x_continuous(
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      expand = c(0, 0)
    ) +

    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-colour_limit, colour_limit),
      oob = scales::squish,
      na.value = "grey90"
    ) +

    coord_fixed(
      expand = FALSE,
      clip = "off"
    ) +

    # Remove legend, labels, axes, title, and margins
    guides(fill = "none") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0, unit = "pt")
    )

  safe_name <- gsub(
    pattern = "[^A-Za-z0-9_-]+",
    replacement = "_",
    x = matrix_names[i]
  )

  ggsave(
    filename = paste0(
      "average_fc_matrix_",
      safe_name,
      ".png"
    ),
    plot = matrix_plot,
    width = 10,
    height = 10,
    units = "in",
    bg = "transparent"
  )
}
