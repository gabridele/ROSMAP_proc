# Visualize average parcelwise functional-connectivity matrices.
#
# Matrix locations are supplied through a two-column manifest. Atlas labels are
# read from the unmodified AtlasPack 4S456 TSV so workstation-specific paths and
# redistributed atlas derivatives are not required.

# -----------------------------------------------------------------------------
# Shared paths
# -----------------------------------------------------------------------------
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
rm(.script_arg, .script_dir, .paths_candidates, .paths_file)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(reticulate)
})

np <- import("numpy", convert = TRUE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
# Set TRUE only *if the atlas TSV is not already ordered by network* and grouped
# network blocks are desired in the figure.
REORDER_BY_NETWORK <- FALSE

# All matrices share the same fixed connectivity scale for direct comparison.
COLOUR_LIMIT <- 1

# Exclude diagonal/self-connections from the descriptive matrix statistics.
EXCLUDE_DIAGONAL_FROM_STATS <- TRUE

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
# Matrix paths are supplied through a small two-column CSV file. The file
# must contain columns named `name` and `path`; see fc_matrix_manifest.example.csv.
# Path to the manifest should be set with ROSMAP_FC_MATRIX_MANIFEST in paths.R

manifest_file <- require_file(
  FC_MATRIX_MANIFEST
)

manifest <- readr::read_csv(
  manifest_file,
  show_col_types = FALSE
)

required_manifest_columns <- c("name", "path")
if (!all(required_manifest_columns %in% names(manifest))) {
  stop("FC matrix manifest must contain columns: name, path")
}

manifest <- manifest %>%
  transmute(
    name = trimws(as.character(name)),
    path = trimws(as.character(path))
  )

if (nrow(manifest) == 0 || any(!nzchar(manifest$name)) || any(!nzchar(manifest$path))) {
  stop("FC matrix manifest contains empty names or paths.")
}
if (anyDuplicated(manifest$name)) {
  stop("FC matrix manifest contains duplicate matrix names.")
}

# Resolve relative paths relative to the manifest itself.
manifest_dir <- dirname(normalizePath(manifest_file, mustWork = TRUE))
is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\])", manifest$path)
manifest$resolved_path <- manifest$path
manifest$resolved_path[!is_absolute] <- file.path(
  manifest_dir,
  manifest$path[!is_absolute]
)

matrix_names <- manifest$name
matrix_files <- vapply(
  manifest$resolved_path,
  require_file,
  character(1),
  label = "FC matrix"
)

# Path to the atlas TSV file. Included in the repository, but can be overidden in paths.R

dseg_file <- require_file(SCHAEFER4S456_ATLAS_TSV, "4S456 atlas TSV")

dseg <- read_tsv(
  dseg_file,
  na = c("", "NA", "N/A", "n/a", "NaN"),
  show_col_types = FALSE
)

required_dseg_columns <- c("network_label", "atlas_name")
if (!all(required_dseg_columns %in% names(dseg))) {
  stop(
    "The atlas TSV must contain columns: ",
    paste(required_dseg_columns, collapse = ", ")
  )
}

network <- dseg %>%
  transmute(
    network_label = na_if(trimws(as.character(network_label)), ""),
    atlas_name = na_if(trimws(as.character(atlas_name)), ""),
    plot_label = coalesce(network_label, atlas_name, "Unlabelled")
  ) %>%
  pull(plot_label)

# -----------------------------------------------------------------------------
# Load and validate matrices
# -----------------------------------------------------------------------------
matrices <- lapply(matrix_files, function(file) {
  x <- as.matrix(np$load(file, allow_pickle = FALSE))

  if (length(dim(x)) != 2 || nrow(x) != ncol(x)) {
    stop(file, " does not contain a square two-dimensional matrix.")
  }
  if (any(!is.finite(x))) {
    warning(file, " contains non-finite values; they will be plotted as missing.")
  }

  x
})

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
    "Atlas/matrix size mismatch: matrix has ", number_of_rois,
    " ROIs, but the atlas TSV contains ", length(network), " rows."
  )
}

# -----------------------------------------------------------------------------
# Common ROI ordering and network blocks
# -----------------------------------------------------------------------------
# TRUE:
#   ROIs are reordered so that members of each network are adjacent.
#
# FALSE:
#   Original atlas/dseg order is preserved.
#
# The same ordering is applied to rows and columns of every matrix.

if (REORDER_BY_NETWORK) {
  network_factor <- factor(network, levels = unique(network))
  roi_order <- order(network_factor)
} else {
  roi_order <- seq_len(number_of_rois)
}

network_plot_order <- network[roi_order]
matrices <- lapply(
  matrices,
  function(x) x[roi_order, roi_order, drop = FALSE]
)

network_runs <- rle(network_plot_order)
network_blocks <- tibble(
  network = network_runs$values,
  block_length = network_runs$lengths
) %>%
  mutate(
    end = cumsum(block_length),
    start = lag(end, default = 0) + 1,
    midpoint = (start + end) / 2
  )

network_boundaries <- head(network_blocks$end, -1) + 0.5

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
matrix_to_dataframe <- function(x, matrix_name) {
  tibble(
    row = rep(seq_len(nrow(x)), times = ncol(x)),
    column = rep(seq_len(ncol(x)), each = nrow(x)),
    connectivity = as.vector(x),
    matrix = matrix_name
  )
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

matrix_summary <- function(x, matrix_name, exclude_diagonal = TRUE) {
  values <- x
  if (exclude_diagonal) diag(values) <- NA_real_
  values <- as.numeric(values)
  values <- values[is.finite(values)]

  if (length(values) == 0) {
    return(tibble(
      matrix = matrix_name,
      n_values = 0L,
      mean = NA_real_,
      sd = NA_real_,
      median = NA_real_,
      min = NA_real_,
      max = NA_real_
    ))
  }

  tibble(
    matrix = matrix_name,
    n_values = length(values),
    mean = mean(values),
    sd = sd(values),
    median = median(values),
    min = min(values),
    max = max(values)
  )
}

make_fc_plot <- function(data, title, facet = FALSE, axis_text_size = 7) {
  p <- ggplot(
    data,
    aes(x = column, y = row, fill = connectivity)
  ) +
    geom_tile() +
    geom_vline(
      xintercept = network_boundaries,
      linewidth = 0.25,
      colour = "black",
      alpha = 0.7
    ) +
    geom_hline(
      yintercept = network_boundaries,
      linewidth = 0.25,
      colour = "black",
      alpha = 0.7
    ) +
    scale_x_continuous(
      position = "top",
      breaks = network_blocks$midpoint,
      labels = network_blocks$network,
      expand = c(0, 0)
    ) +
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
      limits = c(-COLOUR_LIMIT, COLOUR_LIMIT),
      oob = scales::squish,
      na.value = "grey90",
      name = "Connectivity"
    ) +
    coord_fixed(clip = "off") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.ticks = element_blank(),
      axis.text.x.top = element_text(
        size = axis_text_size,
        angle = 45,
        hjust = 0,
        vjust = 0
      ),
      axis.text.y.left = element_text(size = axis_text_size, hjust = 1),
      plot.title = element_text(face = "bold", size = 18),
      plot.margin = margin(t = 15, r = 15, b = 10, l = 15)
    )

  if (facet) {
    p <- p + facet_wrap(~matrix, ncol = 2)
  }

  p
}

make_clean_fc_plot <- function(data) {
  ggplot(
    data,
    aes(x = column, y = row, fill = connectivity)
  ) +
    geom_raster(interpolate = FALSE) +
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
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-COLOUR_LIMIT, COLOUR_LIMIT),
      oob = scales::squish,
      na.value = "grey90"
    ) +
    coord_fixed(expand = FALSE, clip = "off") +
    guides(fill = "none") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0, unit = "pt")
    )
}

# -----------------------------------------------------------------------------
# Individual figures and descriptive stats
# -----------------------------------------------------------------------------
summary_table <- bind_rows(Map(
  function(x, name) {
    matrix_summary(
      x,
      matrix_name = name,
      exclude_diagonal = EXCLUDE_DIAGONAL_FROM_STATS
    )
  },
  matrices,
  matrix_names
))

write_csv(
  summary_table,
  output("fc_matrices", "average_fc_matrix_summary.csv")
)
print(summary_table)

for (i in seq_along(matrices)) {
  plot_data_i <- matrix_to_dataframe(matrices[[i]], matrix_names[i])
  safe_name <- safe_filename(matrix_names[i])

  labelled_plot <- make_fc_plot(
    data = plot_data_i,
    title = paste("Average Functional Connectivity Matrix -", matrix_names[i]),
    facet = FALSE,
    axis_text_size = 12
  )

  clean_plot <- make_clean_fc_plot(plot_data_i)

  print(labelled_plot)

  ggsave(
    filename = output(
      "fc_matrices",
      paste0("average_fc_matrix_", safe_name, ".svg")
    ),
    plot = labelled_plot,
    width = 9,
    height = 8,
    bg = "white"
  )

  ggsave(
    filename = output(
      "fc_matrices",
      paste0("average_fc_matrix_", safe_name, ".png")
    ),
    plot = clean_plot,
    width = 10,
    height = 10,
    units = "in",
    bg = "transparent"
  )
}
