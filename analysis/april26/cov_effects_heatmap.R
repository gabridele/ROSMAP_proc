library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

# ============================================================
# 1. Load baseline CSVs
# ============================================================

csvs_bl <- list.files(
  path = "analysis/april26",
  pattern = "all_pairwise_results_.*conn_covs_long.csv",
  full.names = TRUE
)

covs_bl_list <- lapply(csvs_bl, read_csv)

names(covs_bl_list) <- csvs_bl %>%
  str_replace("all_pairwise_results_", "") %>%
  str_replace("conn_covs_bl.csv", "")

covs_bl_df <- bind_rows(covs_bl_list, .id = "comparison")

# Check what comparison values actually are
print(unique(covs_bl_df$comparison))
print(colnames(covs_bl_df))
print(unique(covs_bl_df$contrast))

# ============================================================
# 2. Prepare data
# ============================================================
get_shared_t_limits <- function(
    df,
    covariate_to_plot,
    filter_statuses = c("pre", "post"),
    padding = 0
) {

  # Allow the function to work before or after t_value is created
  value_column <- if ("t_value" %in% names(df)) {
    "t_value"
  } else if ("statistic" %in% names(df)) {
    "statistic"
  } else {
    stop("The data frame must contain either 't_value' or 'statistic'.")
  }

  t_values <- df %>%
    filter(
      covariate == covariate_to_plot,
      filter_status %in% filter_statuses
    ) %>%
    pull(.data[[value_column]]) %>%
    as.numeric()

  # Remove NA, Inf and -Inf values
  t_values <- t_values[is.finite(t_values)]

  if (length(t_values) == 0) {
    stop(
      paste0(
        "No finite t-values found for covariate = '",
        covariate_to_plot,
        "'."
      )
    )
  }

  # Symmetric limits centred on zero
  max_abs_t <- max(abs(t_values))

  # Avoid identical limits if every t-value is zero
  if (max_abs_t == 0) {
    max_abs_t <- 1
  }

  max_abs_t <- max_abs_t * (1 + padding)

  c(-max_abs_t, max_abs_t)
}

stats_df <- covs_bl_df %>%
  mutate(
    comparison = case_when(
      str_detect(comparison, regex("within", ignore_case = TRUE)) ~ "within",
      str_detect(comparison, regex("between", ignore_case = TRUE)) ~ "between",
      TRUE ~ comparison
    ),
    
    covariate = as.character(covariate),
    contrast = as.character(contrast),
    sig = ifelse(is.na(sig), "", sig),
    t_value = as.numeric(statistic)
  ) %>%
  separate(
    contrast,
    into = c("contrast1", "contrast2"),
    sep = " - ",
    remove = FALSE,
    fill = "right",
    extra = "merge"
  )

stats_df <- stats_df %>%
  separate(
    network_combo,
    into = c("network1_from_combo", "network2_from_combo"),
    sep = "_to_",
    remove = FALSE,
    fill = "right",
    extra = "merge"
  ) %>%
  mutate(
    network1 = case_when(
      comparison == "between" ~ network1_from_combo,
      comparison == "within"  ~ as.character(network),
      TRUE ~ NA_character_
    ),
    network2 = case_when(
      comparison == "between" ~ network2_from_combo,
      comparison == "within"  ~ as.character(network),
      TRUE ~ NA_character_
    ),
    network1 = str_replace(network1, "Network", ""),
    network2 = str_replace(network2, "Network", "")
  ) %>%
  select(-network1_from_combo, -network2_from_combo)

stats_df <- stats_df %>%
  mutate(
    comparison = as.character(comparison),
    covariate = as.character(covariate),
    contrast = as.character(contrast),
    sig = ifelse(is.na(sig), "", sig),
    t_value = as.numeric(statistic)
  ) %>%
  filter(
    !is.na(covariate),
    !is.na(contrast),
    !is.na(t_value),
    !is.na(network1),
    !is.na(network2)
  )
# ============================================================
# 3. Network order
# ============================================================

network_order <- c(
  "Vis",
  "SomMot",
  "DorsAttn",
  "SalVentAttn",
  "Limbic",
  "Cont",
  "Default"
)


# function to make heatmaps for a given covariate and filter status

## try
make_covariate_heatmap_grid_lower <- function(
    df,
    covariate_to_plot,
    filter_status_to_plot = "pre",
    network_order = c(
      "Vis", "SomMot", "DorsAttn",
      "SalVentAttn", "Limbic", "Cont", "Default"
    ),
    text_size = 2.2,
    fill_limits = NULL
) {
  
  # ============================================================
  # 1. Filter to one covariate and one filter status
  # ============================================================
  
  effects_df <- df %>%
    filter(
      covariate == covariate_to_plot,
      filter_status == filter_status_to_plot
    ) %>%
    mutate(
      sig = ifelse(is.na(sig), "", sig),
      t_value = as.numeric(t_value),
      network1 = as.character(network1),
      network2 = as.character(network2),
      contrast = as.character(contrast),
      contrast1 = as.character(contrast1),
      contrast2 = as.character(contrast2)
    )
  
  if (nrow(effects_df) == 0) {
    stop(
      paste0(
        "No rows found for covariate = '", covariate_to_plot,
        "' and filter_status = '", filter_status_to_plot, "'."
      )
    )
  }
  
  # ============================================================
  # 2. Define covariate-level order for the outer grid
  # ============================================================
  
  contrast_levels <- unique(c(effects_df$contrast1, effects_df$contrast2))

  if ("NCI" %in% contrast_levels) {
    contrast_levels <- c(
      "NCI",
      setdiff(contrast_levels, "NCI")
    )
  }  
  effects_df <- effects_df %>%
    mutate(
      contrast1 = factor(contrast1, levels = contrast_levels),
      contrast2 = factor(contrast2, levels = contrast_levels),
      contrast1_id = as.integer(contrast1),
      contrast2_id = as.integer(contrast2)
    )
  
  # ============================================================
  # 3. Force all contrasts into lower-triangle layout
  # ============================================================
  # Displayed contrast is always:
  #
  #   panel_row - panel_col
  #
  # If original contrast is opposite of that, flip the sign of t.
  #
  # Red means panel_row > panel_col.
  # Blue means panel_row < panel_col.
  
  effects_df <- effects_df %>%
    mutate(
      flip_contrast = contrast1_id < contrast2_id,
      
      panel_row = ifelse(
        flip_contrast,
        as.character(contrast2),
        as.character(contrast1)
      ),
      panel_col = ifelse(
        flip_contrast,
        as.character(contrast1),
        as.character(contrast2)
      ),
      
      panel_row_id = pmax(contrast1_id, contrast2_id),
      panel_col_id = pmin(contrast1_id, contrast2_id),
      
      t_value_plot = ifelse(flip_contrast, -t_value, t_value),
      sig_plot = sig
    ) %>%
    filter(panel_row_id > panel_col_id) %>%
    mutate(
      panel_row = factor(panel_row, levels = contrast_levels),
      panel_col = factor(panel_col, levels = contrast_levels)
    )
  
  # ============================================================
  # 4. Build symmetric network x network heatmap inside each panel
  # ============================================================
  
  heatmap_df <- bind_rows(
    effects_df %>%
      transmute(
        panel_row,
        panel_col,
        row_net = network1,
        col_net = network2,
        t_value = t_value_plot,
        sig = sig_plot
      ),
    effects_df %>%
      filter(network1 != network2) %>%
      transmute(
        panel_row,
        panel_col,
        row_net = network2,
        col_net = network1,
        t_value = t_value_plot,
        sig = sig_plot
      )
  ) %>%
    group_by(panel_row, panel_col, row_net, col_net) %>%
    summarise(
      t_value = first(t_value),
      sig = first(sig),
      .groups = "drop"
    ) %>%
    mutate(
      cell_label = ifelse(
        is.na(t_value),
        "",
        paste0(sprintf("%.2f", t_value), sig)
      ),
      row_net = factor(row_net, levels = rev(network_order)),
      col_net = factor(col_net, levels = network_order)
    )
  
  # ============================================================
  # 5. Create full network x network grid for every lower-triangle panel
  # ============================================================
  
  panel_grid <- effects_df %>%
    distinct(panel_row, panel_col)
  
  heatmap_full <- panel_grid %>%
    crossing(
      row_net = factor(rev(network_order), levels = rev(network_order)),
      col_net = factor(network_order, levels = network_order)
    ) %>%
    left_join(
      heatmap_df,
      by = c("panel_row", "panel_col", "row_net", "col_net")
    )
  
  # ============================================================
  # 6. Remove empty first row and empty last column from outer grid
  # ============================================================
  # For lower triangle:
  #   first covariate level has no row comparisons
  #   last covariate level has no column comparisons
  
  first_level <- contrast_levels[1]
  last_level <- contrast_levels[length(contrast_levels)]
  
  heatmap_full <- heatmap_full %>%
    filter(
      panel_row != first_level,
      panel_col != last_level
    ) %>%
    mutate(
      panel_row = factor(
        panel_row,
        levels = contrast_levels[contrast_levels != first_level]
      ),
      panel_col = factor(
        panel_col,
        levels = contrast_levels[contrast_levels != last_level]
      )
    )
  
  # ============================================================
  # 7. Plot
  # ============================================================
  
  p <- ggplot(
    heatmap_full,
    aes(x = col_net, y = row_net, fill = t_value)
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(
      aes(label = cell_label),
      size = text_size,
      color = "black",
      na.rm = TRUE
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = fill_limits,
      na.value = "grey90",
      name = "t-value"
    ) +
    facet_grid(
      rows = vars(panel_row),
      cols = vars(panel_col),
      drop = FALSE
    ) +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 6
      ),
      axis.text.y = element_text(
        size = 6
      ),
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      legend.position = "right",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    ) +
    labs(
      title = paste0(covariate_to_plot, " effect heatmap grid"),
      subtitle = paste0(
        "Lower triangle only. Each panel is row level - column level. ",
        "Red = row > column, blue = row < column."
      )
    )
  
  return(p)
}


make_covariate_heatmaps <- function(
    df,
    covariate_to_plot,
    filter_status_to_plot = "pre",
    network_order = c(
      "Vis", "SomMot", "DorsAttn",
      "SalVentAttn", "Limbic", "Cont", "Default"
    ),
    text_size = 6,
    fill_limits = NULL,
    save_plots = FALSE,
    output_dir = "."
) {
  
  # ------------------------------------------------------------
  # Filter and clean data
  # ------------------------------------------------------------
  
  effects_df <- df %>%
    filter(
      covariate == covariate_to_plot,
      filter_status == filter_status_to_plot
    ) %>%
    mutate(
      sig = ifelse(is.na(sig), "", sig),
      t_value = as.numeric(t_value),
      network1 = as.character(network1),
      network2 = as.character(network2),
      contrast = as.character(contrast),
      contrast1 = as.character(contrast1),
      contrast2 = as.character(contrast2)
    )
  
  if (nrow(effects_df) == 0) {
    stop(
      paste0(
        "No rows found for covariate = '", covariate_to_plot,
        "' and filter_status = '", filter_status_to_plot, "'."
      )
    )
  }
  
  contrast_combinations <- unique(effects_df$contrast)
  
  heatmaps <- list()
  
  # ------------------------------------------------------------
  # Make one heatmap per contrast
  # ------------------------------------------------------------
  
  for (this_contrast in contrast_combinations) {
    
    contrast_df <- effects_df %>%
      filter(contrast == this_contrast)
    
    contrast1_name <- unique(contrast_df$contrast1)[1]
    contrast2_name <- unique(contrast_df$contrast2)[1]
    
    # ----------------------------------------------------------
    # Make symmetric heatmap data
    # ----------------------------------------------------------
    
    heatmap_df <- bind_rows(
      contrast_df %>%
        transmute(
          row_net = network1,
          col_net = network2,
          t_value = t_value,
          sig = sig
        ),
      contrast_df %>%
        filter(network1 != network2) %>%
        transmute(
          row_net = network2,
          col_net = network1,
          t_value = t_value,
          sig = sig
        )
    ) %>%
      group_by(row_net, col_net) %>%
      summarise(
        t_value = first(t_value),
        sig = first(sig),
        .groups = "drop"
      ) %>%
      mutate(
        cell_label = ifelse(
          is.na(t_value),
          "",
          paste0(sprintf("%.2f", t_value), sig)
        ),
        row_net = factor(row_net, levels = rev(network_order)),
        col_net = factor(col_net, levels = network_order)
      )
    
    # ----------------------------------------------------------
    # Add missing cells to force full network x network matrix
    # ----------------------------------------------------------
    
    heatmap_full <- expand_grid(
      row_net = factor(rev(network_order), levels = rev(network_order)),
      col_net = factor(network_order, levels = network_order)
    ) %>%
      left_join(
        heatmap_df,
        by = c("row_net", "col_net")
      )
    
    # ----------------------------------------------------------
    # Plot
    # ----------------------------------------------------------
    
    p <- ggplot(
      heatmap_full,
      aes(x = col_net, y = row_net, fill = t_value)
    ) +
      geom_tile(color = "white", linewidth = 0.7) +
      geom_text(
        aes(label = cell_label),
        size = text_size,
        color = "black",
        na.rm = TRUE
      ) +
      scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = fill_limits,
        na.value = "grey90",
        name = paste0(contrast1_name, " - ", contrast2_name, "\nt-value")
      ) +
      coord_fixed() +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          face = "bold"
        ),
        axis.text.y = element_text(
          face = "bold"
        ),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right"
      ) +
      labs(
        title = paste0(covariate_to_plot, " effect heatmap, ", filter_status_to_plot, " FD filtering"),
        subtitle = paste0(
          "Contrast: ", contrast1_name, " - ", contrast2_name,
          " | red = ", contrast1_name, " > ", contrast2_name,
          ", blue = ", contrast1_name, " < ", contrast2_name
        )
      )
    
    print(p)
    
    heatmaps[[this_contrast]] <- p
    
    # ----------------------------------------------------------
    # Optional save
    # ----------------------------------------------------------
    
    if (save_plots) {
      
      safe_covariate <- str_replace_all(covariate_to_plot, "[^A-Za-z0-9]+", "_")
      safe_contrast <- str_replace_all(this_contrast, "[^A-Za-z0-9]+", "_")
      safe_filter <- str_replace_all(filter_status_to_plot, "[^A-Za-z0-9]+", "_")
      
      ggsave(
        filename = file.path(
          output_dir,
          paste0(
            "heatmap_",
            safe_covariate,
            "_",
            safe_contrast,
            "_",
            safe_filter,
            ".png"
          )
        ),
        plot = p,
        width = 8,
        height = 7,
        dpi = 300
      )
    }
  }
  
  return(heatmaps)
}

# e.g. 
# sex 
msex_t_limits <- get_shared_t_limits(stats_df, "msex")

msex_heatmaps_pre <- make_covariate_heatmaps(
  df = stats_df,
  covariate_to_plot = "msex",
  filter_status_to_plot = "pre",
  fill_limits = msex_t_limits
)

msex_heatmaps_post <- make_covariate_heatmaps(
  df = stats_df,
  covariate_to_plot = "msex",
  filter_status_to_plot = "post",
  fill_limits = msex_t_limits
)
# eyes
eyes_t_limits <- get_shared_t_limits(stats_df, "eyes")

eyes_heatmaps_pre <- make_covariate_heatmaps(
  df = stats_df,
  covariate_to_plot = "eyes",
  filter_status_to_plot = "pre",
  fill_limits = eyes_t_limits
)

eyes_heatmaps_post <- make_covariate_heatmaps(
  df = stats_df,
  covariate_to_plot = "eyes",
  filter_status_to_plot = "post",
  fill_limits = eyes_t_limits
)

syn_t_limits <- get_shared_t_limits(stats_df, "syn_bin")
syn_heatmaps_pre <- make_covariate_heatmaps(
  df = stats_df,
  covariate_to_plot = "syn_bin",
  filter_status_to_plot = "pre",
  fill_limits = syn_t_limits
)

syn_heatmaps_post <- make_covariate_heatmaps(
  df = stats_df,
  covariate_to_plot = "syn_bin",
  filter_status_to_plot = "post",
  fill_limits = syn_t_limits
)


# 2+ levels
site_t_limits <- get_shared_t_limits(stats_df, "site")
dcfdx_t_limits <- get_shared_t_limits(stats_df, "dcfdx")


p_site_pre_grid <- make_covariate_heatmap_grid_lower(
  df = stats_df,
  covariate_to_plot = "site",
  filter_status_to_plot = "pre",
  fill_limits = site_t_limits
)

p_site_post_grid <- make_covariate_heatmap_grid_lower(
  df = stats_df,
  covariate_to_plot = "site",
  filter_status_to_plot = "post",
  fill_limits = site_t_limits
)

print(p_site_pre_grid)
print(p_site_post_grid)

p_dx_pre_grid <- make_covariate_heatmap_grid_lower(
  df = stats_df,
  covariate_to_plot = "dcfdx",
  filter_status_to_plot = "pre"
)

p_dx_post_grid <- make_covariate_heatmap_grid_lower(
  df = stats_df,
  covariate_to_plot = "dcfdx",
  filter_status_to_plot = "post"
)

print(p_dx_pre_grid)
print(p_dx_post_grid)

ggsave(
  filename = "analysis/april26/heatmap_grid_site_pre.pdf",
  plot = p_site_pre_grid,
  width = 18,
  height = 10,
  dpi = 300
)
