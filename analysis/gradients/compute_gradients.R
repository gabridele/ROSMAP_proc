### Scripts for computing gradients are adapted from:
### https://github.com/DeMONLab-BioFINDER/fc_changes_follow_gradients
### Original implementation by @jorittmo.

# ============================================================
# Shared paths
# ============================================================

.script_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

.script_dir <- if (length(.script_arg)) {
  dirname(
    normalizePath(
      sub("^--file=", "", .script_arg[[1]]),
      mustWork = FALSE
    )
  )
} else {
  getwd()
}

.paths_candidates <- unique(
  c(
    file.path("analysis", "paths.R"),
    file.path(.script_dir, "paths.R"),
    file.path(.script_dir, "..", "paths.R"),
    "paths.R",
    file.path("..", "paths.R")
  )
)

.paths_file <- .paths_candidates[
  file.exists(.paths_candidates)
][1]

if (is.na(.paths_file)) {
  stop("Could not locate analysis/paths.R")
}

source(.paths_file)

rm(
  .script_arg,
  .script_dir,
  .paths_candidates,
  .paths_file
)


# ============================================================
# Libraries
# ============================================================

library(tidyverse)
library(SCORPIUS)
library(conflicted)
library(pbapply)
library(reticulate)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::lag)

np <- import("numpy", convert = TRUE)

# ============================================================
#
# ============================================================

source(
  file.path(
    ANALYSIS_DIR,
    "gradients",
    "util_gradients.R"
  )
)

# ============================================================
# Atlas setup
# ============================================================

rois <- read.csv2(
  require_file(
    file.path(
      ANALYSIS_DIR,
      "gradients",
      "atlas_data",
      "Schaefer2018_400Parcels_order.txt"
    )
  ),
  header = TRUE
)

rois <- rois %>%
  rename(
    region = V2
  )

rois$region <- as.character(
  rois$region
)


# Yeo 7-network colours and names
net_names <- data.frame(
  name = c(
    "Vis",
    "SomMot",
    "DorsAttn",
    "SalVentAttn",
    "Limbic",
    "Cont",
    "Default"
  ),

  col = c(
    "#781286",
    "#4682B4",
    "#00760E",
    "#C43AFA",
    "#c7cc7a",
    "#E69422",
    "#CD3E4E"
  ),

  label = 1:7
)


# Yeo-network membership of Schaefer 400 parcels
yeo_msk <- as.numeric(
  readr::read_lines(
    require_file(
      file.path(
        ANALYSIS_DIR,
        "gradients",
        "atlas_data",
        "org_mask_yeo_400.txt"
      )
    ),
    skip = 1
  )
)


# Parcel metadata
roi_data <- data.frame(
  region = rois$region,
  yeo_label = yeo_msk
) %>%
  inner_join(
    net_names,
    join_by(
      yeo_label == label
    )
  )



# Colors for the different ends of the gradients
gradient_cols <- data.frame(gradient1 = c("#3F596D", "#D38A4E"), 
                            gradient2 = c("#4682B4", "#781286"),  
                            gradient3 = c("#8A6081", "#738518"))

# ============================================================
# Python / NumPy
# ============================================================




# ============================================================
# Margulies reference gradients
# ============================================================

marg_gradients <- read_csv(
  require_file(
    MARGULIES_GRADIENTS_400
  ),
  show_col_types = FALSE
) %>%
  select(
    gradient1,
    gradient2,
    gradient3
  )


# ============================================================
# Load site-average connectomes
# ============================================================

gradient_sheet <- read_csv(
  require_file(
    GRADIENT_FC_SHEET
  ),
  show_col_types = FALSE
)

if (!all(c("name", "path") %in% names(gradient_sheet))) {
  stop(
    "Gradient FC sheet must contain columns: name, path"
  )
}


get_matrix_path <- function(site_name) {

  path <- gradient_sheet %>%
    filter(
      name == site_name
    ) %>%
    pull(path)

  if (length(path) != 1) {
    stop(
      "Expected exactly one matrix path for site: ",
      site_name
    )
  }

  require_file(
    path,
    paste0(site_name, " FC matrix")
  )
}


bnk400xcpd <- np$load(
  get_matrix_path("BNK"),
  allow_pickle = FALSE
)

uc400xcpd <- np$load(
  get_matrix_path("UC"),
  allow_pickle = FALSE
)

mg400xcpd <- np$load(
  get_matrix_path("MG"),
  allow_pickle = FALSE
)

rirc400xcpd <- np$load(
  get_matrix_path("RIRC"),
  allow_pickle = FALSE
)


# ============================================================
# Crop to cortical parcels only
# ============================================================

bnk400xcpd <- bnk400xcpd[1:400,1:400]
uc400xcpd <- uc400xcpd[1:400,1:400]
mg400xcpd <- mg400xcpd[1:400,1:400]
rirc400xcpd <- rirc400xcpd[1:400,1:400]

# ============================================================
# Gradient extraction
# ============================================================

grad_list_individ <- get_gradients(
  connectome_ests = list(
    bnk400xcpd = bnk400xcpd,
    uc400xcpd = uc400xcpd,
    mg400xcpd = mg400xcpd,
    rirc400xcpd = rirc400xcpd
  ),

  reference_gradients = marg_gradients,
  n_gradients = c(
    1,
    2,
    3
  ),
  threshold = 0.0,
  similarity_method = "cosine",
  on_affinity = FALSE,
  method = "pca",
  visualize = TRUE
)

# ============================================================
# Save primary results
# ============================================================

write_csv(
  grad_list_individ$gradients,
  output(
    "gradients",
    "primary_gradients.csv"
  )
)

write_csv(
  grad_list_individ$varexp,
  output(
    "gradients",
    "primary_variance_explained.csv"
  )
)


# ============================================================
# Parameter sensitivity analysis
#
# Uncomment to reproduce the exploratory parameter grid.
# ============================================================

# params <- expand_grid(
#   method = c(
#     "pca",
#     "diffusion"
#   ),
#
#   affinity = c(
#     FALSE,
#     TRUE
#   ),
#
#   sim_method = "cosine",
#
#   threshold = c(
#     0.0,
#     0.25,
#     0.5,
#     0.75
#   )
# ) %>%
#
#   filter(
#     !(
#       method == "diffusion" &
#         !affinity
#     )
#   ) %>%
#
#   filter(
#     !(
#       method == "pca" &
#         affinity
#     )
#   ) %>%
#
#   mutate(
#     sim_method = ifelse(
#       !affinity,
#       NA,
#       sim_method
#     )
#   )
#
#
# print(
#   "Calculating gradients over various parameters"
# )
#
#
# gradient_data <- NULL
# varexp_df <- NULL
#
#
# for (i in seq_len(nrow(params))) {
#
#   param_i <- params[i, ]
#
#   grad_list <- get_gradients(
#
#     connectome_ests = list(
#       bnk400xcpd = bnk400xcpd,
#       uc400xcpd = uc400xcpd,
#       mg400xcpd = mg400xcpd,
#       rirc400xcpd = rirc400xcpd
#     ),
#
#     n_gradients = c(
#       1,
#       2,
#       3
#     ),
#
#     threshold = param_i$threshold,
#
#     similarity_method =
#       param_i$sim_method,
#
#     on_affinity =
#       param_i$affinity,
#
#     method =
#       param_i$method,
#
#     visualize = TRUE,
#
#     side_density = FALSE
#   )
#
#
#   gradient_data <- bind_rows(
#     gradient_data,
#     grad_list$gradients
#   )
#
#   varexp_df <- bind_rows(
#     varexp_df,
#     grad_list$varexp
#   )
# }
#
#
# gradient_data <- gradient_data %>%
#   distinct()