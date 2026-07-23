# Install once
install.packages(c("reticulate", "ggplot2", "scales"))

library(reticulate)
library(ggplot2)

# Load NumPy
np <- import("numpy", convert = TRUE)

# ------------------------------------------------------------
# 1. Paths to four averaged connectivity matrices
# ------------------------------------------------------------

files <- c(
  "/Users/ga0034de/Desktop/BNKBBR_ts/bnkbbr_avg_fc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/bnk/avg/avg_bnk_fc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/uc/avg/avg_uc_fc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/mg/avg/avg_mg_fc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/rirc/avg/avg_rirc_fc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/fc_matrices_456/avg/avg_fc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/fc_matrices_456/FD_filtered/avg_fc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/bnk/FD_filtered/avg_fc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/uc/FD_filtered/avg_fc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/mg/FD_filtered/avg_fc_matrix.npy",
  "/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/rirc/FD_filtered/avg_fc_matrix.npy"

)

matrix_names <- c(
  "BNKBBR",
  "BNK",
  "UC",
  "MG",
  "RIRC",
  "whole_dataset",
  "FD_filtered",
  "FD_filtered_BNK",
  "FD_filtered_UC",
  "FD_filtered_MG",
  "FD_filtered_RIRC"
)

# ------------------------------------------------------------
# 2. Load and validate matrices
# ------------------------------------------------------------

matrices <- lapply(files, function(file) {
  x <- np$load(file, allow_pickle = FALSE)
  x <- as.matrix(x)

  if (nrow(x) != ncol(x)) {
    stop(file, " does not contain a square matrix.")
  }

  x
})

# Confirm that all matrices have identical dimensions
matrix_dimensions <- vapply(
  matrices,
  function(x) paste(dim(x), collapse = "x"),
  character(1)
)

if (length(unique(matrix_dimensions)) != 1) {
  stop("The four matrices do not have matching dimensions.")
}

# ------------------------------------------------------------
# 3. Convert matrices to long-format data
# ------------------------------------------------------------

matrix_to_dataframe <- function(x, matrix_name) {
  data.frame(
    row = rep(seq_len(nrow(x)), times = ncol(x)),
    column = rep(seq_len(ncol(x)), each = nrow(x)),
    connectivity = as.vector(x),
    matrix = matrix_name
  )
}

plot_data <- do.call(
  rbind,
  Map(matrix_to_dataframe, matrices, matrix_names)
)

plot_data$matrix <- factor(
  plot_data$matrix,
  levels = matrix_names
)

# Use the same symmetric colour limits for every matrix
colour_limit <- 0.6

# For correlation matrices, you can instead use:
# colour_limit <- 1

# ------------------------------------------------------------
# 4. Plot the four matrices
# ------------------------------------------------------------

fc_plot <- ggplot(
  plot_data,
  aes(x = column, y = row, fill = connectivity)
) +
  geom_tile() +
  facet_wrap(~matrix, ncol = 2) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-colour_limit, colour_limit),
    oob = scales::squish,
    na.value = "grey90",
    name = "Connectivity"
  ) +
  coord_fixed() +
  labs(
    title = "Average Functional Connectivity Matrices",
    x = "ROI",
    y = "ROI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

print(fc_plot)

# # Save to file
ggsave(
  filename = "4sites_average_fc_matrices.png",
  plot = fc_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# each figure one plot
for (i in seq_along(matrices)) {
  single_plot <- ggplot(
    matrix_to_dataframe(matrices[[i]], matrix_names[i]),
    aes(x = column, y = row, fill = connectivity)
  ) +
    geom_tile() +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-colour_limit, colour_limit),
      oob = scales::squish,
      na.value = "grey90",
      name = "Connectivity"
    ) +
    coord_fixed() +
    labs(
      title = paste("Average Functional Connectivity Matrix -", matrix_names[i]),
      x = "ROI",
      y = "ROI"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
    print(single_plot)
  ggsave(
    filename = paste0("average_fc_matrix_", matrix_names[i], ".png"),
    plot = single_plot,
    width = 6,
    height = 5,
    dpi = 300
  )
} 

# Set TRUE to exclude diagonal values from the summary statistics
exclude_diagonal <- TRUE

for (i in seq_along(matrices)) {

  current_matrix <- matrices[[i]]

  # Values used for summary statistics
  stats_values <- current_matrix

  if (exclude_diagonal) {
    diag(stats_values) <- NA
  }

  stats_values <- as.numeric(stats_values)
  stats_values <- stats_values[is.finite(stats_values)]

  # Calculate summary statistics
  matrix_stats <- c(
    Mean   = mean(stats_values),
    SD     = sd(stats_values),
    Min    = min(stats_values),
    Max    = max(stats_values),
    Median = median(stats_values)
  )

  # Text displayed on the plot
  stats_text <- sprintf(
    "Mean: %.3f   SD: %.3f   Min: %.3f   Max: %.3f   Median: %.3f",
    matrix_stats["Mean"],
    matrix_stats["SD"],
    matrix_stats["Min"],
    matrix_stats["Max"],
    matrix_stats["Median"]
  )

  plot_data_i <- matrix_to_dataframe(
    current_matrix,
    matrix_names[i]
  )

  single_plot <- ggplot(
    plot_data_i,
    aes(x = column, y = row, fill = connectivity)
  ) +
    geom_tile() +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-colour_limit, colour_limit),
      oob = scales::squish,
      na.value = "grey90",
      name = "Connectivity"
    ) +
    coord_fixed() +
    labs(
      title = paste(
        "Average Functional Connectivity Matrix -",
        matrix_names[i]
      ),
      subtitle = stats_text,
      x = "ROI",
      y = "ROI",
      caption = if (exclude_diagonal) {
        "Summary statistics exclude diagonal values."
      } else {
        "Summary statistics include diagonal values."
      }
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = element_text(
        size = 10,
        margin = margin(b = 8)
      ),
      plot.caption = element_text(
        size = 8,
        hjust = 0
      )
    )

  print(single_plot)

  # Make filename safe if matrix names contain spaces
  safe_name <- gsub(
    pattern = "[^A-Za-z0-9_-]+",
    replacement = "_",
    x = matrix_names[i]
  )

  ggsave(
    filename = paste0("average_fc_matrix_", safe_name, ".png"),
    plot = single_plot,
    width = 6,
    height = 5,
    dpi = 300
  )
}
