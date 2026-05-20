library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

# ============================================================
# 1. Load baseline CSVs
# ============================================================

csvs_bl <- list.files(pattern = "covs_bl.csv")

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

make_covariate_heatmaps <- function(
    df,
    covariate_to_plot,
    filter_status_to_plot = "pre",
    network_order = c(
      "Vis", "SomMot", "DorsAttn",
      "SalVentAttn", "Limbic", "Cont", "Default"
    ),
    text_size = 6,
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
        title = paste0(covariate_to_plot, " effect heatmap"),
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
site_heatmaps <- make_covariate_heatmaps(
  df = stats_df,
  covariate_to_plot = "site",
  filter_status_to_plot = "post"
)
