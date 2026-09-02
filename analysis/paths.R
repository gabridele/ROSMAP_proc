# Shared configuration for ROSMAP analysis scripts.
#
# Private ROSMAP-derived inputs are intentionally not distributed with this
# repository. Set ROSMAP_DATA to their directory and ROSMAP_OUTPUT to the
# directory where generated tables and figures should be written.
#
# Publication motion-exclusion rule:
#   retain scans only when mean_FD < 0.25 mm.

.find_analysis_dir <- function() {
  explicit <- Sys.getenv("ROSMAP_ANALYSIS_DIR", unset = "")
  if (nzchar(explicit)) {
    return(normalizePath(explicit, mustWork = FALSE))
  }

  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(source_file) && nzchar(source_file)) {
    return(dirname(normalizePath(source_file, mustWork = FALSE)))
  }

  candidates <- c("analysis", ".", "..")
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "paths.R"))) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
  }

  stop(
    "Could not locate analysis/paths.R. Run from the repository root or ",
    "analysis directory, or set ROSMAP_ANALYSIS_DIR."
  )
}

ANALYSIS_DIR <- .find_analysis_dir()
DATA_DIR <- normalizePath(
  Sys.getenv("ROSMAP_DATA", unset = file.path(ANALYSIS_DIR, "sheets")),
  mustWork = FALSE
)
OUTPUT_DIR <- normalizePath(
  Sys.getenv("ROSMAP_OUTPUT", unset = file.path(ANALYSIS_DIR, "outputs")),
  mustWork = FALSE
)

# Canonical publication threshold. All FD exclusion filters use strict '<'.
FD_THRESHOLD <- 0.25

# Build a path to a private input file.
data <- function(...) {
  file.path(DATA_DIR, ...)
}

# Resolve the prepared connectivity/demographics table. An explicit override
# takes precedence; otherwise use analysis/outputs/prepared/<filename>.
demos <- function(default_filename = "demos_conn.csv") {
  override <- Sys.getenv("ROSMAP_DEMOS_CSV", unset = "")
  if (nzchar(override)) {
    return(normalizePath(override, mustWork = FALSE))
  }

  prepared <- file.path(OUTPUT_DIR, "prepared", default_filename)
  if (file.exists(prepared)) {
    return(prepared)
  }

  # Backward-compatible fallbacks for historical private input layouts.
  direct_private <- file.path(DATA_DIR, default_filename)
  if (file.exists(direct_private)) {
    return(direct_private)
  }

  legacy_private <- file.path(DATA_DIR, "v1.3", default_filename)
  if (file.exists(legacy_private)) {
    return(legacy_private)
  }

  prepared
}

# Build an output path and create its parent directory on demand.
output <- function(...) {
  path <- file.path(OUTPUT_DIR, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

# Build an output directory and create it immediately.
output_dir <- function(...) {
  path <- file.path(OUTPUT_DIR, ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# Fail early with a clear message for required files.
require_file <- function(path, label = "input file") {
  if (!file.exists(path)) {
    stop(
      label, " not found: ", path, "\n",
      "Configure inputs as documented in analysis/README.md."
    )
  }
  path
}

## path to csv containing the site-wise fc matrices for plotting
FC_MATRIX_MANIFEST <- Sys.getenv(
  "ROSMAP_FC_MATRIX_MANIFEST",
  unset = file.path(
    ANALYSIS_DIR,
    "fc_related",
    "fc_matrix_manifest.csv"
  )
)

SCHAEFER4S456_ATLAS_TSV <- Sys.getenv(
  "SCHAEFER4S456_ATLAS_TSV",
  unset = file.path(
    ANALYSIS_DIR,
    "fc_related",
    "atlas-4S456Parcels",
    "atlas-4S456Parcels_dseg.tsv"
  )
)