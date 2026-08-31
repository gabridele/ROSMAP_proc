# Shared path configuration for analysis/effect_covs.
#
# The ROSMAP-derived input sheets are intentionally not distributed with this
# repository. Set ROSMAP_DATA to the directory that contains the
# private input input sheets. Generated figures
# and input sheets are written under ROSMAP_OUTPUT.

.find_code_dir <- function() {
  explicit <- Sys.getenv("ROSMAP_DIR", unset = "")
  if (nzchar(explicit)) {
    return(normalizePath(explicit, mustWork = FALSE))
  }

  if (dir.exists(file.path("analysis", "effect_covs"))) {
    return(normalizePath(file.path("analysis", "effect_covs"), mustWork = FALSE))
  }

  normalizePath(getwd(), mustWork = FALSE)
}

DIR <- .find_code_dir()
DATA_DIR <- Sys.getenv(
  "ROSMAP_DATA",
  unset = file.path(DIR, "sheets")
)
OUTPUT_DIR <- Sys.getenv(
  "ROSMAP_OUTPUT",
  unset = file.path(DIR, "outputs")
)

# Build a path to a private input table. The helper does not require the file
# to exist immediately so scripts can define paths before validating inputs.
data <- function(...) {
  file.path(DATA_DIR, ...)
}

# Most analyses historically used a dated demos_conn_*.csv snapshot. Set
# ROSMAP_DEMOS_CSV to override that snapshot without editing code.
demos <- function(default_filename) {
  override <- Sys.getenv("ROSMAP_DEMOS_CSV", unset = "")
  if (nzchar(override)) {
    return(override)
  }
  data("v1.3", default_filename)
}

# Build an output path and create its parent directory on demand.
output <- function(...) {
  path <- file.path(OUTPUT_DIR, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

# Fail early with a clear message for required files.
require_file <- function(path, label = "input file") {
  if (!file.exists(path)) {
    stop(
      label, " not found: ", path, "\n",
      "Set ROSMAP_DATA/ROSMAP_DEMOS_CSV as documented in ",
      "analysis/effect_covs/README.md."
    )
  }
  path
}

# Build an output directory and create it immediately.
output_dir <- function(...) {
  path <- file.path(OUTPUT_DIR, ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}
