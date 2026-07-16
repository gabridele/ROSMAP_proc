library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)

source("~/github_dir/ROSMAP_proc/analysis/gradients/plot_gams.R")
source("~/github_dir/ROSMAP_proc/analysis/gradients/plot_grad_rels.R")

pad_plot <- function(plot, pad = 1) {
  library(grid)
  library(gtable)
  
  gt <- if (inherits(plot, "patchwork")) {
    patchwork::patchworkGrob(plot)
  } else {
    ggplot2::ggplotGrob(plot)
  }
  
  gtable::gtable_add_padding(gt, grid::unit(pad, "pt"))
}

make_gradient_plots <- function(gradient_data,
                                gradient_colors = data.frame(gradient1 = c("#3F596D", "#D38A4E"), 
                                                             gradient2 = c("#4682B4", "#781286"),  
                                                             gradient3 = c("#8A6081", "#738518")),
                                atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds"),
                                region_col = "region",
                                grad_char = c("gradient1", "gradient2", "gradient3"),
                                include_gradient_plots = TRUE,
                                add_shade = FALSE,
                                shade_size = 0.1,
                                shade_alpha = 0.3,
                                grad_name = TRUE,
                                base_size_ = 10) {
  
  gradient_plots <- list()
  
  for (grad in grad_char) {
    
    if (include_gradient_plots) {
      
      plot_data <- gradient_data %>%
        right_join(atlas_geometry$atlas)
      
      p <- ggplot(plot_data) +
        geom_sf(aes(
          fill = .data[[grad]],
          geometry = geometry
        ),
        linewidth = 0.1,
        show.legend = FALSE) +
        (if (add_shade)
          geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
         else NULL) +
        theme_void(base_size = base_size_) +
        labs(
          fill = "",
          title = str_to_title(str_replace(paste0("_", grad), "_", " "))
        ) +
        (if (grad_name)
          labs(title = dplyr::case_when(
            grad == "gradient1" ~ "Sens-Assoc",
            grad == "gradient2" ~ "Vis-Mot",
            grad == "gradient3" ~ "Rep-Exec",
            TRUE ~ grad
          )) else NULL) +
        theme(
          legend.position = "",
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(color = "black", hjust = 0.5)
        ) +
        scale_fill_gradient2(
          low = gradient_colors[[grad]][1],
          mid = "white",
          high = gradient_colors[[grad]][2]
        )
      
      # store plot with simple gradient name
      gradient_plots[[grad]] <- p
    }
  }
  
  return(gradient_plots)
}

prepare_longitudinal_and_window_analysis <- function(long_df,
                                                     run_windowing = FALSE,
                                                     processed_dir = "data/processed_and_cleaned",
                                                     mod_formula = formula(paste0("FC ~ age_bl + time + pathΔ + path_bl + sex + rsqa__MeanFD + (1 | sid)")),
                                                     cache_path = file.path(processed_dir, "longitudinal_window_analysis_data.rds"),
                                                     use_cache = FALSE,
                                                     save_cache = FALSE,
                                                     b_size = 7) {
  library(tidyverse)
  
  if (use_cache && file.exists(cache_path)) {
    return(readRDS(cache_path))
  }
  
  long_bf_ <- long_df
  
  bf_longitudinal <- plot_gradient_relationships(
    long_bf_,
    gradient_data = grad_df %>% filter(study == "biofinder"),
    gradients = c(1, 3),
    empty_row_height = -0.2,
    brain_title_vjust = -10,
    r_spin_size = 0.9,
    gradient_colors = gradient_cols,
    add_shade = TRUE,
    shade_alpha = 0.1,
    shade_size = 0.1,
    point_size = 0.2,
    point_alpha = 0.3,
    rasterize = TRUE,
    ggrastr_dpi = 300,
    ggrastr_dev = "ragg",
    r2_size = rel(3.2),
    base_size_ = b_size,
    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
    covariates = c("time", "sex", "rsqa__MeanFD"),
    filter_criteria = quo(),
    show_networks = FALSE,
    plt_title = NULL,
    tag_prefix = "",
    tag_sep = "",
    layout_construction = "horizontal",
    include_gradient_plots = TRUE,
    right_term_side = FALSE,
    cache_runs = FALSE,
    longitudinal = TRUE,
    sub_id = "sid",
    longitudinal_formula = formula(paste0("FC ~ time  + age_bl + path_bl + pathΔ + sex + rsqa__MeanFD + (1 | sid)"))
  )
  
  num_bins <- 30
  hist_map_data <- long_bf_ %>%
    filter(fmri_bl) %>%
    select(age_bl, path_bl) %>%
    mutate(
      age_bin = cut(age_bl, breaks = num_bins, include.lowest = TRUE, labels = FALSE),
      path_bin = cut(path_bl, breaks = num_bins, include.lowest = TRUE, labels = FALSE)
    ) %>%
    mutate(
      age_mid = min(age_bl) + (age_bin - 0.5) * (max(age_bl) - min(age_bl)) / num_bins,
      path_mid = min(path_bl) + (path_bin - 0.5) * (max(path_bl) - min(path_bl)) / num_bins
    ) %>%
    group_by(age_mid, path_mid) %>%
    summarise(n = n(), .groups = "drop")
  
  age_df <- long_bf_ %>% filter(fmri_bl) %>% select(age_bl)
  age_win_size <- 25
  frames <- data.frame(
    win_n = seq(min(age_df$age_bl) %>% round(), max(age_df$age_bl) %>% round() - age_win_size, by = 1)
  ) %>%
    mutate(
      window_min = win_n,
      window_max = win_n + age_win_size
    ) %>%
    rowwise() %>%
    mutate(
      sample_size = sum(age_df$age_bl >= window_min & age_df$age_bl <= window_max),
      mean_age = mean(age_df$age_bl[age_df$age_bl >= window_min & age_df$age_bl <= window_max])
    ) %>%
    ungroup() %>%
    mutate(window_range = paste0("[", window_min, ", ", window_max, "]"))
  
  if (run_windowing) {
    sliding_window_results <- data.frame()
    for (win_n in frames$win_n) {
      window_df <- long_bf_ %>% filter((age_bl >= win_n), (age_bl <= age_win_size + win_n))
      ests_window <- nodal_lmm_ests(window_df, fc_measures_bf$affinity, roi_names = rois, model_formula = mod_formula)
      ests_window$window <- win_n
      sliding_window_results <- rbind(sliding_window_results, ests_window)
    }
    write_rds(sliding_window_results, file.path(processed_dir, "sliding_window_age_res.rds"))
  } else {
    sliding_window_results <- readRDS(file.path(processed_dir, "sliding_window_age_res.rds"))
  }
  
  grad_cor_age <- sliding_window_results %>%
    filter(!(term %in% c("(Intercept)", "time", "rsqa__MeanFD", "sex"))) %>%
    group_by(window, term, n) %>%
    summarise(
      Gradient1 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient1), statistic),
      Gradient3 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient3), statistic),
      .groups = "drop"
    ) %>%
    inner_join(frames, join_by(window == win_n))
  
  path_df <- long_bf_ %>% filter(fmri_bl) %>% select(path_bl)
  path_win_size <- 0.35
  steps <- seq(min(path_df$path_bl) %>% round(), (max(path_df$path_bl) %>% round() - path_win_size) + 0.01, by = 0.015)
  frames_path <- data.frame(win_n = steps) %>%
    mutate(
      window_min = win_n,
      window_max = win_n + path_win_size
    ) %>%
    rowwise() %>%
    mutate(
      sample_size = sum(path_df$path_bl >= window_min & path_df$path_bl <= window_max),
      mean_path = mean(path_df$path_bl[path_df$path_bl >= window_min & path_df$path_bl <= window_max])
    ) %>%
    ungroup() %>%
    mutate(window_range = paste0("[", window_min, ", ", window_max, "]"))
  
  if (run_windowing) {
    sliding_window_results_path <- data.frame()
    for (win_n in frames_path$win_n) {
      window_df <- long_bf_ %>% filter((path_bl >= win_n), (path_bl <= path_win_size + win_n))
      ests_window <- nodal_lmm_ests(window_df, fc_measures_bf$affinity, roi_names = rois, model_formula = mod_formula)
      ests_window$window <- win_n
      sliding_window_results_path <- rbind(sliding_window_results_path, ests_window)
    }
    write_rds(sliding_window_results_path, file.path(processed_dir, "sliding_window_path_res.rds"))
  } else {
    sliding_window_results_path <- readRDS(file.path(processed_dir, "sliding_window_path_res.rds"))
  }
  
  grad_cor_path <- sliding_window_results_path %>%
    filter(!(term %in% c("(Intercept)", "time", "rsqa__MeanFD", "sex"))) %>%
    group_by(window, term, n) %>%
    summarise(
      Gradient1 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient1), statistic),
      Gradient3 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient3), statistic),
      .groups = "drop"
    ) %>%
    inner_join(frames_path, join_by(window == win_n))
  
  analysis_data <- list(
    bf_longitudinal = bf_longitudinal,
    n_long = bf_longitudinal$n %>% as.numeric(),
    hist_map_data = hist_map_data,
    frames = frames,
    frames_path = frames_path,
    grad_cor_age = grad_cor_age,
    grad_cor_path = grad_cor_path
  )
  
  if (save_cache) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    write_rds(analysis_data, cache_path)
  }
  
  analysis_data
}

plot_longitudinal_and_window_analysis <- function(analysis_data,
                                                  b_size = 7,
                                                  label_size = b_size * 1.2,
                                                  title_size = 8,
                                                  subtitle_size = 6,
                                                  annotation_size = 5,
                                                  axis_title_rel = 1,
                                                  axis_text_rel = 1,
                                                  strip_size_rel = 1,
                                                  plot_linewidth = 0.75,
                                                  struct_linewidth = 0.5,
                                                  point_size = 0.7) {
  library(cowplot)
  library(tidyverse)
  library(ggside)
  
  old <- theme_set(theme_bw(base_size = b_size))
  on.exit(theme_set(old), add = TRUE)
  theme_update(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )
  
  compact_theme <- theme(
    axis.title = element_text(size = rel(axis_title_rel)),
    axis.text = element_text(size = rel(axis_text_rel)),
    strip.text = element_text(size = rel(1)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1))
  )
  
  bf_longitudinal <- analysis_data$bf_longitudinal
  n_long <- analysis_data$n_long
  hist_map_data <- analysis_data$hist_map_data
  frames <- analysis_data$frames
  frames_path <- analysis_data$frames_path
  grad_cor_age <- analysis_data$grad_cor_age
  grad_cor_path <- analysis_data$grad_cor_path
  
  p_bf_long <- bf_longitudinal$plot +
    plot_annotation(
      title = waiver(),
      subtitle = waiver(),
      theme = theme(
        plot.background = element_rect(color = "black", linewidth = struct_linewidth),
        plot.subtitle = element_text(size = rel(0.7), hjust = 0, vjust = -0.05, margin = margin(l = 0.05, unit = "npc")),
        plot.title = element_text(size = rel(0.8), hjust = 0, margin = margin(l = 0.05, unit = "npc"))
      )
    ) &
    #compact_theme &
    theme(plot.tag = element_blank(), plot.title = element_text(size = rel(1)))
  
  p_bf_long[[1]] <- p_bf_long[[1]] + theme(plot.title = element_text(vjust = -2))
  p_bf_long[[2]] <- p_bf_long[[2]] + theme(plot.title = element_text(vjust = -2))
  p_bf_long[[3]] <- p_bf_long[[3]] + labs(title = "Age BL") + theme(plot.title = element_text(vjust = 0, margin = margin(b = -5)))
  p_bf_long[[4]] <- p_bf_long[[4]] + labs(title = "Pathology BL") + theme(plot.title = element_text(vjust = 0, margin = margin(b = -5)))
  p_bf_long[[5]] <- p_bf_long[[5]] + labs(title = bquote(Delta * "Pathology" ~ (t[i] - t[0]))) + theme(plot.title = element_text(vjust = -0.3, margin = margin(b = -5)))
  
  hist_map <- hist_map_data %>%
    ggplot() +
    geom_tile(aes(x = age_mid, y = path_mid, fill = n), color = "black", linewidth = 0.05) +
    annotate("rect", xmin = 18.5, xmax = 88, ymin = 0, ymax = 0.35, fill = NA, color = "black", linewidth = struct_linewidth*0.5) +
    annotate("rect", xmin = 19, xmax = 88.5, ymin = 0.015, ymax = 0.365, fill = NA, linetype = 2, color = "black", linewidth = struct_linewidth*0.5) +
    annotate("rect", xmin = 20, xmax = 45, ymin = -0.025, ymax = 1, fill = NA, color = "black", linewidth = struct_linewidth*0.5) +
    annotate("rect", xmin = 21, xmax = 46, ymin = -0.015, ymax = 1.01, fill = NA, linetype = 2, color = "black", linewidth = struct_linewidth*0.5) +
    annotate("curve", x = 90, y = 0.20, xend = 90, yend = 0.365, arrow = arrow(length = unit(0.05, "inches")), linewidth = struct_linewidth) +
    annotate("curve", x = 35, y = 1.025, xend = 46, yend = 1.025, arrow = arrow(length = unit(0.05, "inches")), curvature = -0.3, linewidth = struct_linewidth) +
    annotate("text", x = 40, y = 1.1, label = "1 Year", size = annotation_size) +
    annotate("text", x = 95, y = 0.305, label = "0.015 Path", angle = -90, size = annotation_size) +
    theme_gray(base_size = b_size) +
    scale_x_continuous("Age BL") +
    scale_y_continuous("Pathology BL", breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
    scale_fill_continuous("Count", breaks = c(2, 5, 8)) +
    compact_theme +
    theme(
      legend.position = "top",
      legend.box.margin = margin(-7, -0, -7, -7),
      legend.key.height = unit(0.01, "npc"),
      legend.key.width = unit(0.03, "npc"),
      legend.justification = "right",
      legend.title = element_text(vjust = 0),
      legend.text.position = "top",
      legend.text = element_text(vjust = -1),
      legend.background = element_rect(fill = NA),
      plot.background = element_blank(),
      panel.border = element_rect(linewidth = struct_linewidth, color = "black", fill = NA)
    )
  
  make_window_plot <- function(grad_cor, frames_df, x_lab, label_step, include_legend = FALSE, supp = FALSE) {
    grad_cor %>%
      pivot_longer(starts_with("Grad"), names_to = "gradient", values_to = "grad_corr") %>%
      { if (!supp) filter(., (gradient == "Gradient3" & term == "age_bl") | (gradient == "Gradient1" & term == "path_bl") | (gradient == "Gradient1" & term == "pathΔ")) else . } %>%
      mutate(term = case_when(
        term == "age_bl" ~ "Age~at~baseline",
        term == "path_bl" ~ "Pathology~at~baseline",
        term == "pathΔ" ~ "Delta*Pathology~(t[i]-t[0])",
        TRUE ~ term
      ) |> factor(levels = c("Age~at~baseline", "Pathology~at~baseline", "Delta*Pathology~(t[i]-t[0])"))) %>%
      ggplot(aes(window, grad_corr, color = gradient, fill = gradient)) +
      geom_point(aes(alpha = n), size = point_size, show.legend = FALSE) +
      geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), linewidth = plot_linewidth) +
      geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5, linewidth = struct_linewidth) +
      facet_wrap(~term, nrow = 1, labeller = label_parsed) +
      theme_bw(base_size = b_size) +
      scale_x_continuous(
        labels = frames_df$window_range[seq(1, length(frames_df$win_n), by = label_step)],
        breaks = frames_df$win_n[seq(1, length(frames_df$win_n), by = label_step)],
        guide = guide_axis(check.overlap = TRUE)
      ) +
      labs(y = "Gradient corr (r)", x = x_lab) +
      geom_xsidecol(aes(y = if (supp && x_lab == "Baseline Pathology Window") sample_size / 2 else sample_size, fill = NULL, color = NULL), linewidth = 0.05) +
      ggside(x.pos = "bottom") +
      scale_xsidey_continuous(breaks = c(0, 250)) +
      scale_color_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
      scale_fill_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
      scale_alpha_continuous(range = c(0.3, 1)) +
     # compact_theme +
      theme(
        legend.position = if (include_legend) "inside" else "",
        legend.position.inside = if (include_legend) c(0.12, -0.36) else NULL,
        legend.direction = "horizontal",
        legend.text = element_text(size = rel(1), margin = margin(r = 8, l = 4, unit = "pt")),
        legend.key.size = unit(0.35, "cm"),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text = element_text(size = rel(strip_size_rel)),
        axis.title = element_text(size = rel(axis_title_rel)),
        legend.title = element_blank(),
        ggside.panel.scale = 0.2,
        axis.text.x = element_text(angle = -30, hjust = 0, size = rel(1))
      )
  }
  
  grad_cor_p_age <- make_window_plot(grad_cor_age, frames, "Baseline Age Window", label_step = 5)
  grad_cor_p_path <- make_window_plot(grad_cor_path, frames_path, "Baseline Pathology Window", label_step = 10, include_legend = TRUE)
  
  window_plots <- wrap_plots(list(grad_cor_p_age, grad_cor_p_path), nrow = 2) +
    plot_annotation(theme = theme(plot.background = element_rect(color = "black", linewidth = struct_linewidth)))
  
  get_net_legend <- function() {
    x <- grad_df %>%
      filter(study == "biofinder") %>%
      ggplot(aes(gradient1, gradient3, color = name)) +
      geom_point(alpha = 0.5, size = point_size) +
      labs(color = "Yeo Network") +
      guides(color = guide_legend(label.hjust = 0, byrow = TRUE, nrow = 2, reverse = TRUE, override.aes = list(size = point_size * 1.5))) +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe) +
      theme(
        legend.position = "bottom",
        legend.key.spacing.x = unit(0, "cm"),
        legend.key.spacing.y = unit(-0.8, "cm"),
        legend.direction = "horizontal",
        legend.text.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.8), margin = margin(l = -4, r = 0, unit = "pt")),
        legend.background = element_blank()
      )
    ggpubr::as_ggplot(ggpubr::get_legend(x))
  }
  
  net_legend <- get_net_legend()
  
  p_bf_long_cp <- ggdraw() +
    draw_plot(p_bf_long) +
    draw_plot(net_legend, x = 0.1, y = 0.01, width = 0.25, height = 0.065) +
    draw_plot_label("A", size = label_size) +
    draw_label(paste0("BioFINDER Longitudinal", " (N=", n_long, ")"), x = 0.045, y = 0.965, hjust = 0, size = title_size) +
    draw_label("FCS ~ age_bl + path_bl + Delta path +\ntime + sex + motion + (1|sub)", x = 0.015, y = 0.89, hjust = 0, size = subtitle_size)
  
  window_plots_cp <- ggdraw() +
    draw_plot(window_plots) +
    draw_plot_label("C", size = label_size)
  
  lmm_tab <- magick::image_read("paper/figures/conceptual_plot/LMM_table.png")
  
  slide_meth <- ggdraw() +
    draw_plot(hist_map, x = 0.01, y = 0.3, width = 0.98, height = 0.7) +
    draw_plot(p_bf_long[[5]] + labs(title = expression(Delta * Path)) + theme(plot.title = element_text(vjust = 0, margin = margin(b = 0))), x = 0.39, y = 0.017, width = 0.29, height = 0.29) +
    draw_plot(p_bf_long[[1]] + theme(plot.title = element_text(vjust = 0)), x = 0.69, y = 0.025, width = 0.29, height = 0.29) +
    draw_image(lmm_tab, x = 0.02, y = 0.00, width = 0.35, height = 0.25) +
    draw_plot_label("B", size = label_size) +
    annotate("curve", x = 0.54, y = 0.09, xend = 0.83, yend = 0.09, linewidth = struct_linewidth, arrow = arrow(length = unit(0.06, "inches"), ends = "both")) +
    annotate("segment", x = 0.33, y = 0.425, xend = 0.225, yend = 0.22, linewidth = struct_linewidth, arrow = arrow(length = unit(0.06, "inches"))) +
    annotate("label", label = "Parcel-wise LMM\nin age window", x = 0.25, y = 0.3, size = annotation_size) +
    annotate("text", label = "Pearson Correlation", x = 0.7, y = 0.02, size = annotation_size) +
    annotate("segment", x = 0.38, y = 0.147, xend = 0.43, yend = 0.147, arrow = arrow(length = unit(0.06, "inches")), linewidth = struct_linewidth) +
    theme(plot.background = element_rect(color = "black", linewidth = struct_linewidth))
  
  final_fig <- ggdraw() +
    draw_plot(p_bf_long_cp, x = 0, y = 0.51, width = 0.645, height = 0.49) +
    draw_plot(slide_meth, x = 0.655, y = 0.51, width = 0.345, height = 0.49) +
    draw_plot(window_plots_cp, x = 0, y = 0.0, width = 1, height = 0.5) +
    annotate("curve", x = 0.89, y = 0.5105, xend = 0.702, yend = 0.408, linewidth = 1.2, color = "white", curvature = 0.29) +
    annotate("curve", x = 0.89, y = 0.5105, xend = 0.702, yend = 0.408, linewidth = struct_linewidth, arrow = arrow(length = unit(0.06, "inches")), curvature = 0.29)
  
  grad_cor_p_path_supp <- make_window_plot(grad_cor_path, frames_path, "Baseline Pathology Window", label_step = 10, include_legend = TRUE, supp = TRUE)
  grad_cor_p_age_supp <- make_window_plot(grad_cor_age, frames, "Baseline Age Window", label_step = 5, supp = TRUE)
  
  window_plots_both_G <- wrap_plots(list(grad_cor_p_age_supp, grad_cor_p_path_supp), nrow = 2) +
    plot_annotation(theme = theme(plot.background = element_rect(color = "black", linewidth = struct_linewidth)))
  
  list(main_fig = final_fig, supp_fig = window_plots_both_G)
}

longitudinal_and_window_analysis <- function(long_df,
                                             b_size = 7,
                                             run_windowing = FALSE,
                                             processed_dir = "data/processed_and_cleaned",
                                             mod_formula = formula(paste0("FC ~ age_bl + time + pathΔ + path_bl + sex + rsqa__MeanFD + (1 | sid)")),
                                             analysis_data = NULL,
                                             use_analysis_cache = FALSE,
                                             save_analysis_cache = FALSE,
                                             analysis_cache_path = file.path(processed_dir, "longitudinal_window_analysis_data.rds"),
                                             ...) {
  if (is.null(analysis_data)) {
    analysis_data <- prepare_longitudinal_and_window_analysis(
      long_df = long_df,
      run_windowing = run_windowing,
      processed_dir = processed_dir,
      mod_formula = mod_formula,
      cache_path = analysis_cache_path,
      use_cache = use_analysis_cache,
      save_cache = save_analysis_cache,
      b_size = b_size
    )
  }
  
  figs <- plot_longitudinal_and_window_analysis(
    analysis_data = analysis_data,
    b_size = b_size,
    ...
  )
  
  c(figs, list(analysis_data = analysis_data))
}


figure_one <- function(subject_data,
                       subject_data_replication,
                       measures_list, measures_list_replication, gradients_df = grad_df %>% filter(study=="biofinder"),
                       gradients_df_replication = grad_df %>% filter(study=="adni"),
                       fig1_formula_bf = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"),
                       fig1_formula_ad = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"),
                       brain_plot_names_f1 = NULL,
                       raster = FALSE,
                       ggrastr_dpi = NULL,
                       rast_scale_atlas = 1,
                       rast_scale_shade = 1,
                       tag_size = 21,
                       p_size = 1,
                       draw_size = 18,
                       shade = FALSE,
                       shade_a = 0.01,
                       shade_s = 0.01,
                       b_size = 11,
                       plot_title_size = 1,
                       brain_title_size = rel(plot_title_size),
                       brain_title_vjust = 0,
                       axes_title_size = 1,
                       ax_txt_size = rel(1),
                       r2_sizing1 = 1,
                       r2_sizing2 = 1,
                       boxed = FALSE,
                       empt_row_height = 0,
                       selected_gradients = c(1, 3), 
                       split = FALSE,
                       fig = NULL) {
  library(cowplot)
  # library(showtext)
  # showtext_opts(dpi = 300)
  
  if(!is.list(measures_list)) stop("measures should be in a list")
  if(length(measures_list)>1) stop("list of measures should only contain a single metric")
  if(is.null(names(measures_list))) stop("measures list should be named")
  
  if(!is.list(measures_list_replication)) stop("measures should be in a list")
  if(length(measures_list_replication)>1) stop("list of measures should only contain a single metric")
  if(is.null(names(measures_list_replication))) stop("measures list should be named")
  
  
  l_marg = 0.05
  #text_size = 28
  
  # if (split & fig = "1") {
  #   
  # }
  
  fig1 <- function(){ 
    bf_p <-  plot_gradient_relationships(subject_data %>% filter(fmri_bl), 
                                         gradient_data = gradients_df, 
                                         gradients = selected_gradients,
                                         gradient_colors = gradient_cols,
                                         list_of_parcel_data = measures_list,
                                         empty_row_height = empt_row_height,
                                         base_size_ = b_size,
                                         add_shade = shade,
                                         shade_alpha = shade_a,
                                         shade_size = shade_s,
                                         point_size = p_size,
                                         vect = TRUE,
                                         mod_formula = fig1_formula_bf,
                                         covariates = c("sex", "rsqa__MeanFD"),
                                         #r2_size = rel(r2_sizing1),
                                         r_spin_size = r2_sizing1,
                                         filter_criteria = quo(),
                                         show_networks = FALSE,
                                         axis_text_size = ax_txt_size,
                                         brain_title_size = brain_title_size,
                                         brain_title_vjust = brain_title_vjust,
                                         rasterize = raster,
                                         ggrastr_dpi = ggrastr_dpi,
                                         raster_scale_atlas = rast_scale_atlas,
                                         raster_scale_shade = rast_scale_shade,
                                         tag_prefix = "",
                                         tag_sep = "",
                                         layout_construction = "horizontal",
                                         plt_subtitle = TRUE,
                                         include_gradient_plots = TRUE,
                                         right_term_side = FALSE,
                                         plt_title = "",
                                         cache_runs = FALSE)
    
    n_cs <- bf_p$n
    p_bf <- bf_p$plot &
      theme(#text = element_text(size = text_size),
        plot.tag = element_blank(),
        axis.title = element_text(size = rel(axes_title_size))) 
    
    
    # f <- fig1_formula_bf
    # rhs <- deparse(f[[2]])
    # expr_str <- paste0("italic(FCS[parcel] ~ '~' ~ ", rhs, ")")
    # subtit_expr <- parse(text = expr_str)[[1]]
    
    p_bf <- p_bf + plot_annotation(title = paste0("BioFINDER", " (N=", n_cs, ")"),
                                   #subtitle = subtit_expr,
                                   theme = theme(
                                     plot.subtitle = element_text(size = rel(0.9),
                                                                              hjust = 0,
                                                                              vjust = -0.05,
                                                                              family = "mono",
                                                                              face = "italic",
                                                                              margin = margin(l = l_marg, unit = "npc")),
                                                 plot.title.position = "plot",
                                                 plot.title = element_text(size = rel(1), hjust =0, margin = margin(l = l_marg, unit = "npc")))
    )
    
    if (!is.null(brain_plot_names_f1)) {
      plt_idx_start <- length(selected_gradients) 
      
      for (i in seq_along(brain_plot_names_f1)) {
        p_bf[[plt_idx_start + i]] <- p_bf[[plt_idx_start + i]] + (if (!is.na(brain_plot_names_f1[i])) labs(title = brain_plot_names_f1[i]) else NULL)
      }
    }

    
    p_bf <- ggdraw() + 
      draw_plot(p_bf) + 
      draw_plot_label("A", x = l_marg-l_marg, size = tag_size) 
    #draw_label("BioFINDER", x = (1/3.1/2), y = 0.75, hjust = 0.5, size =  draw_size)
    
    rect_linewidth = 0.5
    
    if (boxed) p_bf <- p_bf + theme(plot.background = element_rect(color = "black", linewidth = rect_linewidth))
    
    
    adni_p <-  plot_gradient_relationships(subject_data_replication, 
                                           gradient_data = gradients_df_replication, 
                                           gradients = selected_gradients,
                                           gradient_colors = gradient_cols,
                                           list_of_parcel_data = measures_list_replication,
                                           mod_formula = fig1_formula_ad,
                                           empty_row_height = empt_row_height,
                                           base_size_ = b_size,
                                           vect = TRUE,
                                           add_shade = shade,
                                           shade_alpha = shade_a,
                                           shade_size = shade_s,
                                           point_size = p_size,
                                           #r2_size = rel(r2_sizing1),
                                           r_spin_size = r2_sizing1,
                                           covariates = c("sex", "rsqa__MeanFD"),
                                           id_var = "file_func",
                                           filter_criteria = quo(),
                                           show_networks = FALSE,
                                           axis_text_size = ax_txt_size,
                                           brain_title_size = brain_title_size,
                                           brain_title_vjust = brain_title_vjust,
                                           rasterize = raster,
                                           ggrastr_dpi = ggrastr_dpi,
                                           raster_scale_atlas = rast_scale_atlas,
                                           raster_scale_shade = rast_scale_shade,
                                           tag_prefix = "",
                                           layout_construction = "horizontal",
                                           include_gradient_plots = TRUE,
                                           right_term_side = FALSE,
                                           plt_title = "",
                                           plt_subtitle = TRUE,
                                           cache_runs = FALSE)
    
    
    n_adni <- adni_p$n
    p_a <- adni_p$plot &
      theme(#text = element_text(size = text_size),
        plot.tag = element_blank(),
        axis.title = element_text(size = rel(axes_title_size)))
    
    # f <- fig1_formula_ad
    # rhs <- deparse(f[[2]])
    # expr_str <- paste0("italic(FCS[parcel] ~ '~' ~ ", rhs, ")")
    # subtit_expr_ad <- parse(text = expr_str)[[1]]
    
    p_a <- p_a + plot_annotation(title = paste0("ADNI", " (N=", n_adni, ")"), 
                                 #subtitle = subtit_expr_ad,
                                 theme = theme(
                                   plot.subtitle = element_text(size = rel(0.9),
                                                                            hjust = 0,
                                                                            vjust = -0.05,
                                                                            family = "mono",
                                                                            face = "italic",
                                                                            margin = margin(l = l_marg, unit = "npc")),
                                               plot.title.position = "plot",
                                               plot.title = element_text(size = rel(1), hjust =0, margin = margin(l = l_marg, unit = "npc")))
    )
    
    
    if (!is.null(brain_plot_names_f1)) {
      plt_idx_start <- length(selected_gradients) 
      
      for (i in seq_along(brain_plot_names_f1)) {
        p_a[[plt_idx_start + i]] <- p_a[[plt_idx_start + i]] + (if (!is.na(brain_plot_names_f1[i])) labs(title = brain_plot_names_f1[i]) else NULL)
      }
    }
    
    p_a <- ggdraw() + 
      draw_plot(p_a) + 
      draw_plot_label("B", x = l_marg-l_marg, size = tag_size) 
    #draw_label("ADNI", x = (1/3.1/2), y = 0.75, hjust = 0.5, size =  draw_size)
    
    if (boxed) p_a <- p_a + theme(plot.background = element_rect(color = "black", linewidth = rect_linewidth))
    return(list(p_bf = p_bf, p_a = p_a, tmaps = list(bf = bf_p$tmaps, adni = adni_p$tmaps)))
  }
  

  
  # overlay <- ggdraw() +
  #   draw_plot(p_bf, x = 0, y = 0, width = 0.35, height = 1) +
  #   draw_plot(p_bf_long, x = 0.27, y = 0, width = 0.35, height = 1) +
  #   draw_plot(p_a, x = 0.64, y = 0, width = 0.35, height = 1)  
  
  
  ######
  # Cognition
  #######
  
  fig2 <- function() {
    health_cog <-  plot_gradient_relationships(subject_data %>% filter(fmri_bl, diagnosis=="Normal" | diagnosis=="SCD", abnorm_ab==0, !apoe4
                                                                       ),
                                               gradient_data = gradients_df, 
                                               gradients = selected_gradients,
                                               gradient_colors = gradient_cols,
                                               list_of_parcel_data = measures_list,
                                               empty_row_height = empt_row_height,
                                               base_size_ = b_size,
                                               #r2_size = rel(r2_sizing2),
                                               r_spin_size = r2_sizing2,
                                               mod_formula = formula(paste0("~ scale(age) * scale(-mPACC_v1) + pathology_ad + sex + rsqa__MeanFD")),
                                               logistic_fit = FALSE,
                                               vect = TRUE,
                                               add_shade = shade,
                                               shade_alpha = shade_a,
                                               shade_size = shade_s,
                                               point_size = p_size,
                                               covariates = c("sex", "rsqa__MeanFD"),
                                               filter_criteria = quo(),
                                               show_networks = FALSE,
                                               axis_text_size = ax_txt_size,
                                               brain_title_size = brain_title_size,
                                               brain_title_vjust = brain_title_vjust,
                                               rasterize = raster,
                                               ggrastr_dpi = ggrastr_dpi,
                                               raster_scale_atlas = rast_scale_atlas,
                                               raster_scale_shade = rast_scale_shade,
                                               tag_prefix = "",
                                               tag_sep = "",
                                               layout_construction = "horizontal",
                                               plot_spacing = 0.2,
                                               include_gradient_plots = TRUE,
                                               right_term_side = FALSE,
                                               plt_title = "",
                                               plt_subtitle = TRUE,
                                               cache_runs = FALSE)
    
    n_health <- health_cog$n
    health_l_marg <- l_marg - 0.016666
    p_health_cog <- health_cog$plot &
      theme(plot.tag = element_blank(),
            axis.title = element_text(size = rel(axes_title_size))
      )
    
    p_health_cog <- p_health_cog + plot_annotation(title = bquote(
      "Cognitively unimpaired w/o APOE " * epsilon * 4 ~
        "(A" * beta^"-" * ") (N=" * .(n_health) * ")"
    ), 
    subtitle = expression(italic(FCS[parcel] ~ "~" ~ age * "×" * "-mPACC" ~ + ~ "AD pathology" ~ + ~ sex ~ + ~ motion)),
    theme = theme(plot.subtitle = element_text(size = rel(0.9),
                                               hjust = 0, 
                                               vjust = -0.05, 
                                               family = "mono",
                                               face = "italic",
                                               margin = margin(l = health_l_marg, unit = "npc")), 
                  plot.title.position = "plot",
                  plot.title = element_text(size = rel(1), hjust =0, margin = margin(l = health_l_marg, unit = "npc")))
    )
    
    plt_idx <- 3:6
    if (length(selected_gradients) < 2) plt_idx <- plt_idx-1
    p_health_cog[[plt_idx[1]]] <- p_health_cog[[plt_idx[1]]] + labs(title = "Age") + theme(plot.title = element_text(vjust = brain_title_vjust))
    p_health_cog[[plt_idx[2]]] <- p_health_cog[[plt_idx[2]]] + labs(title = "-mPACC", subtitle = "(Inverted cognition)") + theme(plot.subtitle = element_text(hjust = 0.5, size = rel(0.6)))
    p_health_cog[[plt_idx[3]]] <- p_health_cog[[plt_idx[3]]] + labs(title = "AD Pathology") + theme(plot.title = element_text(vjust = brain_title_vjust))
    p_health_cog[[plt_idx[4]]] <- p_health_cog[[plt_idx[4]]] + labs(title = "-mPACC×Age") + theme(plot.title = element_text(vjust = brain_title_vjust))
    
    p_health_cog <- ggdraw() + draw_plot(p_health_cog) + 
      draw_plot_label(ifelse(split, "A", "C"), x = health_l_marg-health_l_marg, size = tag_size) + 
      draw_label("BioFINDER", x = (1/5/2), y = 0.75, hjust = 0.5, size =  draw_size)
    
    
    
    clinical_cog <-  plot_gradient_relationships(subject_data %>% filter(fmri_bl, diagnosis=="MCI" | diagnosis=="AD", !is.na(mPACC_v1)) %>% 
                                                   mutate(`-mPACC_v1` = -mPACC_v1), 
                                                 gradient_data = gradients_df, 
                                                 gradients = selected_gradients,
                                                 gradient_colors = gradient_cols,
                                                 list_of_parcel_data = measures_list,
                                                 empty_row_height = empt_row_height,
                                                 base_size_ = b_size,
                                                 #r2_size = rel(r2_sizing2),
                                                 r_spin_size = r2_sizing2,
                                                 vect = TRUE,
                                                 add_shade = shade,
                                                 shade_alpha = shade_a,
                                                 shade_size = shade_s,
                                                 point_size = p_size,
                                                 mod_formula = formula(paste0(" ~ age + pathology_ad + `-mPACC_v1` +  sex + rsqa__MeanFD")),
                                                 logistic_fit = FALSE,
                                                 covariates = c("sex", "rsqa__MeanFD"),
                                                 filter_criteria = quo(),
                                                 show_networks = FALSE,
                                                 axis_text_size = ax_txt_size,
                                                 brain_title_size = brain_title_size,
                                                 brain_title_vjust = brain_title_vjust,
                                                 rasterize = raster,
                                                 ggrastr_dpi = ggrastr_dpi,
                                                 raster_scale_atlas = rast_scale_atlas,
                                                 raster_scale_shade = rast_scale_shade,
                                                 tag_prefix = "",
                                                 tag_sep = "",
                                                 layout_construction = "horizontal",
                                                 include_gradient_plots = FALSE,
                                                 right_term_side = FALSE,
                                                 plt_title = "",
                                                 cache_runs = FALSE)
    
    n_clin <- clinical_cog$n
    clin_l_marg = 0.3
    p_clinical_cog <-
      clinical_cog$plot  &
      theme(plot.tag = element_blank(),
            axis.title = element_text(size = rel(axes_title_size))
      )
    
    p_clinical_cog <- p_clinical_cog + plot_annotation(title = bquote("Diagnosed MCI/AD (A" * beta^"+" * ") (N=" * .(n_clin) * ")"),#paste0("Diagnosed MCI/AD (Aβ+)", " (N=", n_clin, ")"), 
                                                       subtitle = expression(italic(FCS[parcel] ~ "~" ~ age ~ + ~ "AD pathology" ~ + ~ "-mPACC" ~ + ~ sex ~ + ~ motion)),
                                                       theme = theme(plot.subtitle = element_text(size = rel(0.9), 
                                                                                                  hjust = 0, 
                                                                                                  vjust = -0.05, 
                                                                                                  family = "mono",
                                                                                                  face = "italic",
                                                                                                  margin = margin(l = clin_l_marg, unit = "npc")), 
                                                                     plot.title.position = "plot",
                                                                     plot.title = element_text(size = rel(1), hjust =0, margin = margin(l = clin_l_marg, unit = "npc")))
    )
    
    plt_idx <- 3:5
    if (length(selected_gradients) < 2) plt_idx <- plt_idx-1
    p_clinical_cog[[plt_idx[1]]] <- p_clinical_cog[[plt_idx[1]]] + labs(title = "Age") + theme(plot.title = element_text(vjust = brain_title_vjust))
    p_clinical_cog[[plt_idx[2]]] <- p_clinical_cog[[plt_idx[2]]] + labs(title = "AD Pathology") + theme(plot.title = element_text(vjust = brain_title_vjust))
    p_clinical_cog[[plt_idx[3]]] <- p_clinical_cog[[plt_idx[3]]] + labs(title = "-mPACC", subtitle = "(Inverted cognition)") + theme(plot.subtitle = element_text(hjust = 0.5, size = rel(0.6)))
    
    p_clinical_cog <- ggdraw() + draw_plot(p_clinical_cog) + draw_plot_label(ifelse(split, "B", "D"), x = clin_l_marg-l_marg, size = tag_size)
    return(list(p_health_cog = p_health_cog, p_clinical_cog = p_clinical_cog, tmaps = list(health = health_cog$tmaps, clin = clinical_cog$tmaps)))
  }
  
  
  
  
  
  # Everyting, everywhere all at once
  get_net_legend <- function(){
    #scale_factor <- 5
    x <- grad_df %>% filter(study == "biofinder") %>% 
      ggplot(aes(gradient1, gradient3, color = name)) +
      geom_point(alpha = 0.5) +
      labs(color = "Yeo Network") +
      theme_bw(base_size = b_size) +
      guides(color = guide_legend(
        label.hjust=0,
        byrow = TRUE,
        nrow = 1, 
        override.aes = list(size = rel(1.3))
      )) +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe) +
      theme(
        legend.position = "bottom",
        #legend.key.size = unit(0.0, "cm"),
        legend.key.spacing.x = unit(0.15, "cm"),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.direction = "horizontal",
        #legend.title.position = "",
        legend.text.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.8), margin = margin(l = 0.5, r = 0.75, unit = "pt")),
        legend.background = element_blank()
      )
    
    leg <- ggpubr::get_legend(x)
    leg <- ggpubr::as_ggplot(leg)
    leg 
  }
  
  net_legend1 <- get_net_legend()
  
  get_net_legend <- function(){
    #scale_factor <- 5
    x <- grad_df %>% filter(study == "biofinder") %>% 
      mutate(name = factor(name, levels = c("DorsAttn", "SomMot", "SalVentAttn", "Default", "Limbic", "Cont", "Vis"))) %>% 
      ggplot(aes(gradient1, gradient3, color = name)) +
      geom_point(alpha = 0.5) +
      theme_bw(base_size = b_size) +
      labs(color = "Yeo Network") +
      guides(color = guide_legend(
        label.hjust=0,
        byrow = TRUE,
        nrow = 2, 
        reverse = FALSE,
        override.aes = list(size = rel(1))
      )) +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe) +
      theme(
        legend.position = "bottom",
        #legend.key.size = unit(0.0, "cm"),
        legend.key.spacing.x = unit(0, "cm"),
        legend.key.spacing.y = unit(-0.8, "cm"),
        legend.direction = "horizontal",
        #legend.title.position = "",
        legend.text.position = "right",
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = rel(0.7), margin = margin(l = -3, r = 0, unit = "pt")),
        legend.background = element_blank()
      )
    
    leg <- ggpubr::get_legend(x)
    leg <- ggpubr::as_ggplot(leg)
    leg 
  }
  net_legend2 <- get_net_legend()
  
  
  if (split & is.null(fig)) {
    fig_one_two <- list()
    figure1_org <- fig1()
    figure2_org <- fig2()
    
    if (boxed) {
      
      fig_one_two[[1]] <- ggdraw() +
        draw_plot(figure1_org$p_bf, x = 0, y = 0.05, width = 0.49, height = 0.95) +
        #draw_plot(p_bf_long, x = 0.26, y = 0.05, width = 0.36, height = 0.95) +
        draw_plot(figure1_org$p_a, x = 0.51, y = 0.05, width = 0.49, height = 0.95) +
        draw_plot(net_legend1, x = 0.3,  y = 0.015, width = 0.4, height = 0.015)
      
    } else {
      
      fig_one_two[[1]] <- ggdraw() +
        draw_plot(figure1_org$p_bf, x = 0, y = 0.05, width = 0.5, height = 0.95) +
        #draw_plot(p_bf_long, x = 0.26, y = 0.05, width = 0.36, height = 0.95) +
        draw_plot(figure1_org$p_a, x = 0.5, y = 0.05, width = 0.5, height = 0.95) +
        draw_plot(net_legend1, x = 0.4,  y = 0.025, width = 0.3, height = 0.015)
    }
    
    
    if (boxed) {
      fig_one_two[[2]] <- ggdraw() +
        draw_plot(figure2_org$p_health_cog, x = 0.0, y = 0.05, width = 0.6, height = 0.95) +
        draw_plot(figure2_org$p_clinical_cog, x = 0.5, y = 0.05, width = 0.49, height = 0.95) +
        draw_plot(net_legend2, x = 0.2,  y = 0.025, width = 0.4, height = 0.015) +
        theme(plot.background = element_rect(color = "black"))
    } else {
      fig_one_two[[2]] <- ggdraw() +
        draw_plot(figure2_org$p_health_cog, x = 0.0, y = 0.00, width = 0.6, height = 1) +
        draw_plot(figure2_org$p_clinical_cog, x = 0.5, y = 0.00, width = 0.49, height = 1) +
        draw_plot(net_legend2, x = 0.01,  y = 0.035, width = 0.2, height = 0.015) 
    }
    
    return(list(fig_one_two, tmaps = list(fig1 = figure1_org$tmaps, fig2 = figure2_org$tmaps)))
    
  } else if (split & !is.null(fig)) { 
    
    if (fig == 1) {
      figure1_org <- fig1()
      
      if (boxed) {
        
        x <- ggdraw() +
          draw_plot(figure1_org$p_bf, x = 0, y = 0.05, width = 0.49, height = 0.95) +
          #draw_plot(p_bf_long, x = 0.26, y = 0.05, width = 0.36, height = 0.95) +
          draw_plot(figure1_org$p_a, x = 0.51, y = 0.05, width = 0.49, height = 0.95) +
          draw_plot(net_legend1, x = 0.3,  y = 0.015, width = 0.4, height = 0.015)
        
      } else {
        
        x <- ggdraw() +
          draw_plot(figure1_org$p_bf, x = 0, y = 0.05, width = 0.5, height = 0.95) +
          #draw_plot(p_bf_long, x = 0.26, y = 0.05, width = 0.36, height = 0.95) +
          draw_plot(figure1_org$p_a, x = 0.5, y = 0.05, width = 0.5, height = 0.95) +
          draw_plot(net_legend1, x = 0.4,  y = 0.025, width = 0.3, height = 0.015)
      }
      
      return(list(x, tmaps = figure1_org$tmaps))
    }
    
    if (fig == 2) {
      figure2_org <- fig2()
      
      if (boxed) {
        x <- ggdraw() +
          draw_plot(figure2_org$p_health_cog, x = 0.0, y = 0.05, width = 0.6, height = 0.95) +
          draw_plot(figure2_org$p_clinical_cog, x = 0.5, y = 0.05, width = 0.49, height = 0.95) +
          draw_plot(net_legend2, x = 0.4,  y = 0.025, width = 0.4, height = 0.015) +
          theme(plot.background = element_rect(color = "black"))
      } else {
        x <- ggdraw() +
          draw_plot(figure2_org$p_health_cog, x = 0.0, y = 0.00, width = 0.6, height = 1) +
          draw_plot(figure2_org$p_clinical_cog, x = 0.5, y = 0.00, width = 0.49, height = 1) +
          draw_plot(net_legend2, x = 0.03,  y = 0.035, width = 0.2, height = 0.015) 
      }
      
      return(list(x, tmaps = figure2_org$tmaps))
    }
    
  }  else {
    
    figure1_org <- fig1()
    figure2_org <- fig2()
    
    if (boxed) {
      
      
      bottom_plots <- ggdraw() +
        draw_plot(figure2_org$p_health_cog, x = 0.0, y = 0.0, width = 0.6, height = 1) +
        draw_plot(figure2_org$p_clinical_cog, x = 0.5, y = 0.0, width = 0.49, height = 1) +
        draw_plot(net_legend2, x = 0.03,  y = 0.025, width = 0.2, height = 0.015) +
        theme(plot.background = element_rect(color = "black"))
      
      figure1 <- ggdraw() +
        draw_plot(figure1_org$p_bf, x = 0, y = 0.505, width = 0.495, height = 0.495) +
        #draw_plot(p_bf_long, x = 0.26, y = 0.50, width = 0.36, height = 0.5) +
        draw_plot(figure1_org$p_a, x = 0.505, y = 0.505, width = 0.495, height = 0.495) +  
        draw_plot(bottom_plots, x = 0.0, y = 0.0, width = 1, height = 0.49) 

    } else {
      figure1 <- ggdraw() +
        draw_plot(figure1_org$p_bf, x = 0, y = 0.50, width = 0.5, height = 0.5) +
        #draw_plot(p_bf_long, x = 0.26, y = 0.50, width = 0.36, height = 0.5) +
        draw_plot(figure1_org$p_a, x = 0.5, y = 0.50, width = 0.5, height = 0.5) +  
        draw_plot(figure2_org$p_health_cog, x = 0.0, y = 0.01, width = 0.6, height = 0.5) +
        draw_plot(figure2_org$p_clinical_cog, x = 0.5, y = 0.01, width = 0.49, height = 0.5) +
        draw_plot(net_legend, x = 0.4,  y = 0.00, width = 0.4, height = 0.015)
    }
  
    
    
    list(figure1, tmaps = list(fig1 =figure1_org$tmaps, fig2 = figure2_org$tmaps))
    
  }
  
}




net_names_plot <- function(net_names, text_size = 3.8, vertical = TRUE, tile_height = 0.9, alpha_ = 0.8) {
  get_text_color <- function(bg_color) {
    rgb <- col2rgb(bg_color)
    luminance <- (0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]) / 255
    ifelse(luminance > 0.5, "black", "white")
  }
  
  net_names$text_color <- get_text_color(net_names$col)
  if (vertical) {
    net_names$x <- factor(net_names$name, levels = net_names$name)
    net_names$y <- 1
  } else {
    net_names$y <- factor(net_names$name, levels = net_names$name)
    net_names$x <- 1
  }
  
  
  ggplot(net_names, aes(x = x, y = y)) +
    geom_tile(aes(fill = col), alpha = alpha_, width = 0.95, height = tile_height, show.legend = FALSE) +
    geom_text(aes(label = name, color = text_color), size = text_size) +
    scale_fill_identity() +
    scale_color_identity() +
    theme_void() +
    #coord_fixed(ratio = 1 / spacing) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    )
}



overlaid_main_results <- function(subject_data, fc_matrix, 
                                  atlas_geometry = readRDS("data/atlas_data/Schaefer2018_1000Parcels_geometry.rds"),
                                  network_geometry = readRDS("data/atlas_data/Yeo2011_7_geometry.rds"),
                                  base_size_ = 7,
                                  net_width = 0.5,
                                  parc_width = 0.1,
                                  leg_size = 1,
                                  leg_text_size = 1){
  require(scales)
  require(ggside)
  
  x <- nodal_regression_fits(subject_data %>% filter(fmri_bl),
                             fc_matrix,
                             roi_names = rois,
                             logistic = FALSE,
                             vectorised = TRUE,
                             id_var = "image_file",
                             dep_var = "FC",
                             model_formula = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"))
  
  est <- get_nodal_ests(x)
  
  p_age <-  est %>%
    filter(term == "age") %>%
    inner_join(atlas_geometry, by = "region") %>%
    ggplot() +
    geom_sf(aes(
      fill = statistic,
      geometry = geometry), linewidth = parc_width,
      show.legend = FALSE)+
    theme_void(base_size = base_size_)+
    labs(fill = "T-value", title = "Age"
    ) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.title.position = "top",
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(color = "black", hjust = 0.5),
          #legend.key.width = unit(1, "null"),
          plot.margin = margin(4,4,4,4, "mm")
    ) +
    #guides(color = guide_legend(override.aes = list(size = 1))) +
    scale_fill_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red") 
    ) +
    geom_sf(data = network_geometry %>% drop_na(),
            aes(
              #fill = region,
              color = name,
              geometry = geometry
            ), alpha = 0, linewidth = net_width, fill = NA,
            show.legend = FALSE) +
    scale_color_manual(
      values = net_names %>% select(name, col) %>% deframe(),
      guide = "none"
    )
  
  
  range_stat <- est %>% filter(term == "age") %>% pull(statistic) %>% range
  y_range_df <- data.frame(statistic = seq(range_stat[1], range_stat[2], length.out = 100))
  
  hist_age <- est %>% filter(term == "age") %>% 
    inner_join(grad_df %>% filter(study=="biofinder")) %>% 
    select(statistic, region, gradient3) %>% 
    inner_join(roi_data) %>% 
    ggplot(aes(x = fct_rev(fct_reorder(region, yeo_label)) #region
               , y = statistic, fill = name)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = net_names %>% select(name, col) %>% deframe()) +
    labs(y = "T-value", x = "Parcels") +
    theme_gray(base_size = base_size_) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      ggside.panel.scale = 0.05,
      panel.background = element_rect(fill = "gray80", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35)  
    ) +
    ggnewscale::new_scale_fill() +
    geom_xsidetile(aes(x = region, y = 0, fill = gradient3), show.legend = FALSE) +
    scale_fill_gradient2(
      low = gradient_cols[1, 3],
      mid = "white",
      high = gradient_cols[2, 3]
    ) +
    ggnewscale::new_scale_fill() +
    geom_ysidetile(data = y_range_df,
                   aes(x = 0, y = statistic, fill = statistic), show.legend = FALSE, inherit.aes = FALSE) +
    scale_fill_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red")
    ) +
    theme_ggside_void() +
    ggside(x.pos = "bottom", y.pos = "left")
  
  
  p_path <- est %>%
    filter(term == "pathology_ad") %>%
    inner_join(atlas_geometry, by = "region") %>%
    ggplot() +
    geom_sf(aes(
      fill = statistic,
      geometry = geometry), linewidth= parc_width,
      show.legend = FALSE)+
    theme_void(base_size = base_size_)+
    labs(fill = "T-value", title = "AD pathology"
    ) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.title.position = "top",
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(color = "black", hjust = 0.5),
          #legend.key.width = unit(1, "null"),
          plot.margin = margin(4,4,4,4, "mm")
    )  +
    scale_fill_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red")
    )  + geom_sf(data = network_geometry %>% drop_na(),
                 aes(
                   #fill = region,
                   color = name,
                   geometry = geometry
                 ), alpha = 0, linewidth = net_width,
                 show.legend = FALSE) +
    scale_color_manual(
      values = net_names %>% select(name, col) %>% deframe(),
      guide = "none"
    )
  
  range_stat <- est %>% filter(term == "pathology_ad") %>% pull(statistic) %>% range
  y_range_df <- data.frame(statistic = seq(range_stat[1], range_stat[2], length.out = 100))
  
  hist_path <- est %>% filter(term == "pathology_ad") %>% 
    inner_join(grad_df %>% filter(study=="biofinder")) %>% 
    select(statistic, region, gradient1) %>% 
    inner_join(roi_data) %>% 
    arrange(yeo_label) %>% 
    ggplot(aes(x = fct_rev(fct_reorder(region, yeo_label))  #region#fct_rev(fct_reorder(region, gradient1))
               , y = statistic, fill = name)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = net_names %>% select(name, col) %>% deframe()) +
    labs(y = "", x = "Parcels") +
    theme_gray(base_size = base_size_) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      ggside.panel.scale = 0.05,
      panel.background = element_rect(fill = "gray80", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35)  
    ) +
    ggnewscale::new_scale_fill() +
    geom_xsidetile(aes(x = region, y = 0, fill = gradient1), show.legend = FALSE) +
    scale_fill_gradient2(
      low = gradient_cols[1, 1],
      mid = "white",
      high = gradient_cols[2, 1]
    ) +
    ggnewscale::new_scale_fill() +
    geom_ysidetile(data = y_range_df,
                   aes(x = 0, y = statistic, fill = statistic), show.legend = FALSE, inherit.aes = FALSE) +
    scale_fill_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red")
    ) +
    theme_ggside_void() +
    ggside(x.pos = "bottom", y.pos = "left")
  
  x <- ggplot() +
    geom_sf(data = network_geometry %>% drop_na(),
            aes(
              #fill = region,
              color = name,
              geometry = geometry
            ), alpha = 0, linewidth = 0.5,
            show.legend = TRUE) +
    scale_color_manual(
      values = net_names %>% select(name, col) %>% deframe()
    ) +  
    theme_bw(base_size = base_size_) +
    guides(color = guide_legend(
      label.hjust=0,
      byrow = TRUE,
      nrow = 1, 
      override.aes = list(size = rel(leg_size))
    )) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(0.2, "cm"),
      legend.key.spacing.x = unit(0.1, "cm"),
      legend.direction = "horizontal",
      #legend.title.position = "",
      legend.text.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(leg_text_size), margin = margin(l = 2, r = 4, unit = "pt")),
      legend.background = element_blank()
    )
  
  x <- ggpubr::get_legend(x)
  
  library(cowplot)
  
  p <- ggdraw() +
    draw_plot(p_age, height = 0.6, width = 0.5, x = 0.025, y = 0.425) +
    draw_plot(p_path, height = 0.6, width = 0.5, x = 0.525, y = 0.425) + 
    draw_plot(hist_age, height = 0.35, width = 0.475, x = 0.01, y = 0.1) +
    draw_plot(hist_path, height = 0.35, width = 0.475, x = 0.51, y = 0.1) +
    draw_plot(x, height = 0.05, width = 0.33, x = 0.33, y = 0.025)
  
  p <- p + theme(plot.background = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  
  source_data <- est %>% filter(term %in% c("age", "pathology_ad")) %>% 
    inner_join(grad_df %>% filter(study=="biofinder")) %>% 
    select(statistic, region, gradient1) %>% 
    inner_join(roi_data) %>% 
    arrange(yeo_label)
  
  list(plot = p, source_data = source_data)
}



plot_grads_over_params <- function(connectome_list, 
                                   atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds"),
                                   param_grid = NULL, n_gradients = 1:3, base_size_=18) {

  require(sf)
  require(ggside)
  require(ggpmisc)
  require(patchwork)
  
  old <- theme_set(theme_bw(base_size = base_size_))
  theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA),
               legend.background = element_rect(fill = "transparent", colour = NA),
               legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                          col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), #"#DCF8A4"
                          label = c(1:7))

  
  connectome_est_list <- connectome_list
  
  if (is.null(param_grid)){
    params <- expand_grid(method = c("pca", "diffusion"), affinity = c(FALSE, TRUE), sim_method = "cosine", threshold = c(0.0, 0.25, 0.5, 0.75)) %>% 
      filter(!(method == "diffusion" & !affinity)) %>% 
      filter(!(method == "pca" & affinity)) %>% 
      mutate(sim_method = ifelse(!affinity, NA, sim_method))
  } else {
    params <- param_grid
  }
  
  gradient_data <- c()
  varexp_df <- c()
  for (i in 1:nrow(params)) {
    param_i <- params[i, ]
    grad_list <- get_gradients(connectome_ests = connectome_est_list,
                               n_gradients = n_gradients,
                               threshold = param_i$threshold,
                               similarity_method = param_i$sim_method,
                               on_affinity = param_i$affinity,
                               method = param_i$method,
                               visualize = FALSE,
                               side_density = FALSE)
    
    gradient_data <- rbind(gradient_data, grad_list$gradients)
    varexp_df <- rbind(varexp_df, grad_list$varexp)
  }
  
  gradient_data <- gradient_data %>% distinct()
  grad_char <- paste0("gradient", n_gradients)
  
  std_grad_plots <- list()
  
  for (grad in grad_char) {
    
    low_col <- gradient_cols[1, grad]
    hi_col <- gradient_cols[2, grad]
    
    std_grad_plots[[grad]] <- gradient_data %>% filter(study=="margulies") %>% 
      mutate(region = rois) %>% 
      right_join(atlas_geometry$atlas, by = "region") %>%
      ggplot() +
      geom_sf(aes(
        fill = .data[[grad]],
        geometry = geometry), linewidth= 0.2,
        show.legend = FALSE)+
      geom_sf(data = atlas_geometry$shade, size = 0.01, alpha = 0.01) +
      theme_void(base_size = base_size_)+
      labs(fill = "", title = "Margulies") +
      theme(legend.position = "",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(color = "black", hjust = 0.5),
            plot.margin = margin(t = 0,  
                                 r = 0, 
                                 b = 0,  
                                 l = 0)
      ) +
      scale_x_continuous(position = "top")+
      scale_fill_gradient2(
        low = low_col,
        mid = "white",
        high = hi_col
      ) 
    
    
  }
  
  grad_text <- list()
  
  for(g_txt in paste0("G", n_gradients)) {
    grad_text[[g_txt]] <- 
      ggplot() +
      theme_void() +  
      annotate("text", x = 0.5, y = 0.5, label = g_txt, size = rel(20), 
               fontface = "bold",
               color = "#323232",
               hjust = 0.5, vjust = 0.5)
  }
  
  
  
  gradient_plots <- list()
  for (stud in names(connectome_est_list)){
    plt_idx <- 0
    for (grad in grad_char) {
      
      low_col <- gradient_cols[1, grad]
      hi_col <- gradient_cols[2, grad]
      
      for (params_row in 1:nrow(params) ){
        
        pars <- params[params_row, ]
        
        plot_pars <- pars %>% mutate(affinity = ifelse(is.na(sim_method), "No", sim_method), 
                                     method = ifelse(method=="pca", "PCA", "Diffusion")) %>% select(-sim_method) %>% 
          rename(thresh = threshold)
        par_char <- paste(names(plot_pars), plot_pars[1, ], sep = ": ") %>% str_to_sentence() 
        
        
        variance_explained <- varexp_df %>% semi_join(pars, by = join_by(method, affinity, sim_method, threshold)) %>% 
          filter(study == stud)
        
        
        plt_idx <- plt_idx + 1
        p <- gradient_data %>% 
          filter(study == stud) %>% 
          semi_join(pars, by = join_by(method, affinity, sim_method, threshold)) %>% 
          right_join(atlas_geometry$atlas, by = "region") %>%
          ggplot() +
          geom_sf(aes(
            fill = .data[[grad]],
            geometry = geometry), linewidth= 0.2,
            show.legend = FALSE)+
          geom_sf(data = atlas_geometry$shade, size = 0.01, alpha = 0.01) +
          theme_void(base_size = base_size_) +
          labs(fill = "", title = paste0(plot_pars[, "method"],", ", par_char[3]),
               subtitle = ifelse(grad == grad_char[1], pars[1, "threshold"], NA),
               y = ifelse(plot_pars[, "method"] == "PCA" & plot_pars[, "thresh"] == 0, "BioFINDER", ""),
               caption = paste0(round(variance_explained %>% pull(grad)*100), "% var exp")
          ) +
          theme(legend.position = "",
                plot.caption = element_text(vjust = 8),
                panel.background = element_rect(fill = "transparent", colour = NA),
                plot.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                axis.title.y = element_text(angle=90),
                plot.title = element_blank(),
                plot.subtitle = element_text(hjust = 0.5),
                plot.margin = margin(t = 0,  
                                     r = 0, 
                                     b = -10,  
                                     l = 0)
          ) +
          scale_fill_gradient2(
            low = low_col,
            mid = "white",
            high = hi_col 
          ) 
        
        if (grad != grad_char[1]) p <- p + theme(plot.title = element_blank(),
                                                 plot.subtitle = element_blank())
        gradient_plots[[stud]][[plt_idx]] <- p
        
      }
    }
  }
  
  
  plots <- list()
  
  ref_grads <- gradient_data %>% filter(study == "margulies")
  
  for (stud in names(connectome_est_list)){
    plt_idx <- 0
    plot_data_stud <- gradient_data %>% filter(study %in% c(stud))
    for (grad in grad_char) {
      
      low_col <- gradient_cols[1, grad]
      hi_col <- gradient_cols[2, grad]
      
      plot_data_grad <- plot_data_stud %>% pivot_longer(starts_with("gradient"), names_to = "gradient", values_to = "value") %>% 
        filter(gradient == grad) %>% 
        pivot_wider(names_from = "study", values_from = "value")
      
      for (params_row in 1:nrow(params) ){
        pars <- params[params_row, ]
        
        plot_data <- plot_data_grad %>% semi_join(pars) %>% 
          mutate(margulies = ref_grads %>% pull(grad))
        
        x_range_df <- data.frame(x = seq(min(plot_data[, stud]), max(plot_data[, stud]), length.out = 100))
        
        y_range_df <- data.frame(margulies = seq(min(plot_data[, "margulies"]), max(plot_data[, "margulies"]), length.out = 100))
        
        p <- plot_data %>% 
          ggplot(aes(x = .data[[stud]], y = margulies,
                     color = name)) +
          geom_point(alpha = 0.2) +
          stat_poly_eq(color = "#323232", label.x = "left", label.y = "top", size = 6.5) +
          stat_poly_line(se = FALSE, color = "#323232") +
          labs(
            #title = str_to_title(grad),
            y = "",
            x = str_to_title(stud),
            #tag = tag_labs[count],
            color = "Network") +
          theme(
            legend.position = "",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            ggside.panel.scale = 0.05,
            plot.margin = margin(t = 0,  
                                 r = 5, 
                                 b = 0,  
                                 l = 5),
            axis.text.x = element_text(angle =45, hjust = 1)) +
          scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
          scale_y_continuous(limits = c(min(plot_data["margulies"]), max(plot_data["margulies"])))+
          scale_x_continuous(limits = c(min(plot_data[stud]), max(plot_data[stud])), guide = guide_axis(check.overlap = TRUE))
        
        y_side_switch <-  ifelse((pars %>% pull(method) == "pca" & pars %>% pull(threshold) == 0), 1, 0)
        
        p <- p +
          ggnewscale::new_scale_fill() +
          geom_xsidetile(data = x_range_df, 
                         aes(x = x, y = 0, fill = x), 
                         show.legend = FALSE, inherit.aes = FALSE) +
          scale_fill_gradient2(
            low = low_col,
            mid = "white",
            high = hi_col
          ) + ggnewscale::new_scale_fill() +
          geom_ysidetile(data = y_range_df,
                         aes(x = 0, y = margulies, fill = margulies), 
                         alpha = y_side_switch,
                         show.legend = FALSE, 
                         inherit.aes = FALSE) +
          scale_fill_gradient2(
            low = low_col,
            mid = "white",
            high = hi_col
          )
        
        
        p <-  p +
          theme_ggside_void() +
          ggside(x.pos = "bottom", y.pos = "left")
        
        plt_idx <- plt_idx + 1
        plots[[stud]][[plt_idx]] <- p
      }
    }
  }
  
  
  no_of_grads = length(n_gradients)
  n_study <- length(connectome_est_list)
  n_plot_cols <- nrow(params)#((no_of_grads*n_study) + (n_study - 1))
  row_heights <- head(rep(c(1, -0.2, 1, 0.2), no_of_grads), -1)
  
  param_combos <- length(unique(params$threshold))
  n_total_cols <- n_plot_cols + floor(n_plot_cols / param_combos)
  col_widths <- head(c(1, rep(c(rep(1, param_combos), 0.3), floor(n_plot_cols / param_combos))), -1)
  
  layout <- c(
    area(3, 1))
  
  for (g in seq(7, no_of_grads + 8, by = 4)) {
    layout <- c(layout, area(g, 1))
  }
  
  for (g_txt in seq(1, no_of_grads + 8, by = 4) ) {
    layout <- c(layout, area(g_txt, 1))
  }
  
  
  rows <- layout$t - 1
  
  
  for (row in rows) {
    if (row %% 2 == 0) next
    for (col in 2:(n_total_cols + 1)) {
      if (((col - 2) %% (param_combos + 1)) == param_combos) next  
      layout <- c(layout, area(row, col))
    }
  }
  
  for (row in rows + 1) {
    if (row %% 2 == 0) next
    for (col in 2:(n_total_cols + 1)) {
      if (((col - 2) %% (param_combos + 1)) == param_combos) next  
      layout <- c(layout, area(row, col))
    }
  }
  
  
  
  
  patch_plots <- list()
  for (stud in names(connectome_est_list)){
    pp <- 
      std_grad_plots[[1]] + std_grad_plots[2:no_of_grads] +
      grad_text +
      plots[[stud]][1:(nrow(params)*no_of_grads)] +
      gradient_plots[[stud]][1:(nrow(params)*no_of_grads)] +
      plot_layout(design = layout, axis_titles = "collect", axes = "collect_y", guides = "collect",
                  heights = row_heights,
                  widths = col_widths
      ) & 
      theme(legend.position = "",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            axis.title.x =  element_blank(),
            plot.tag.position  = c(.9, .96))  
    
    patch_plots[[stud]] <- pp
  }
  
  
  return(patch_plots)
}



longitudinal_and_window_analysis_legacy <- function(long_df, 
                                                    b_size = 18,
                                                    run_windowing = FALSE,
                                                    processed_dir = "data/processed_and_cleaned",
                                                    mod_formula = formula(paste0("FC ~ age_bl + time + pathΔ + path_bl + sex + rsqa__MeanFD + (1 | sid)"))) {
  
  
  library(cowplot)
  library(tidyverse)
  source("src/util_vis.R")
  
  long_bf_ <- long_df
  
  old <- theme_set(theme_bw(base_size = b_size))
  theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA),
               legend.background = element_rect(fill = "transparent", colour = NA),
               legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  #######################
  # Linear longitudinal
  #######################
  
  
  bf_longitudinal <-  plot_gradient_relationships(long_bf_, 
                                                  gradient_data = grad_df %>% filter(study=="biofinder"), 
                                                  gradients = c(1, 3),
                                                  empty_row_height = -0.1,
                                                  gradient_colors = gradient_cols,
                                                  add_shade = TRUE, 
                                                  shade_alpha = 0.0075,
                                                  shade_size = 0.001,
                                                  rasterize = TRUE,
                                                  ggrastr_dpi = 300,
                                                  ggrastr_dev = "ragg",
                                                  r2_size = rel(3.8),
                                                  base_size_ = b_size,
                                                  list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                  covariates = c("time", "sex", "rsqa__MeanFD"),
                                                  filter_criteria = quo(),
                                                  show_networks = FALSE,
                                                  plt_title = NULL,
                                                  tag_prefix = "",
                                                  tag_sep = "",
                                                  layout_construction = "horizontal",
                                                  include_gradient_plots = TRUE,
                                                  right_term_side = FALSE,
                                                  cache_runs = FALSE,
                                                  longitudinal = TRUE,
                                                  sub_id = "sid",
                                                  longitudinal_formula = formula(paste0("FC ~ time  + age_bl + path_bl + pathΔ + sex + rsqa__MeanFD + (1 | sid)")))
  
  n_long <- bf_longitudinal$n %>% as.numeric()
  long_l_marg <- 0
  p_bf_long <- bf_longitudinal$plot + 
    plot_annotation(title = waiver(), #paste0("Longitudinal", " (N=", n_long, ")"), 
                    subtitle = waiver(), #bquote(FC[parcel] ~ "~" ~ age_bl + path_bl + Δpath + time + sex + motion +"(1|sub)"),
                    theme = theme(plot.background = element_rect(color = "black"),
                                  plot.subtitle = element_text(size = rel(0.8), 
                                                               hjust = 0, 
                                                               vjust = -0.05, 
                                                               margin = margin(l = 0.05, unit = "npc")), 
                                  plot.title = element_text(size = rel(1), 
                                                            hjust =0, margin = margin(l = 0.05, unit = "npc")))
    )  &
    theme(#text = element_text(size = text_size),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 16),
      plot.tag = element_blank())
  
  p_bf_long[[1]] <- p_bf_long[[1]] + theme(plot.title = element_text(vjust = -3))
  p_bf_long[[2]] <- p_bf_long[[2]] + theme(plot.title = element_text(vjust = -3))
  p_bf_long[[3]] <- p_bf_long[[3]] + labs(title = "Age BL") + theme(plot.title = element_text(vjust = -0.0))
  p_bf_long[[4]] <- p_bf_long[[4]] + labs(title = "Pathology BL") + theme(plot.title = element_text(vjust = -0.0))
  p_bf_long[[5]] <- p_bf_long[[5]] + labs(title = bquote(Δ* "Pathology" ~ (t[i]-t[0]) ) ) + theme(plot.title = element_text(vjust = -0.5))
  
  ###########################
  # Windowing
  ###########################
  
  long_bf_ %>% filter(fmri_bl) %>% 
    select(age_bl, path_bl) %>% 
    mutate(age_bl = cut(age_bl, 30),
           path_bl = cut(path_bl, 30)) %>% 
    group_by(age_bl, path_bl) %>% 
    summarise(n = n(), .groups = "drop")  -> x
  
  num_bins <- 30
  
  long_bf_ %>% 
    filter(fmri_bl) %>% 
    select(age_bl, path_bl) %>% 
    mutate(
      age_bin = cut(age_bl, breaks = num_bins, include.lowest = TRUE, labels = FALSE),
      path_bin = cut(path_bl, breaks = num_bins, include.lowest = TRUE, labels = FALSE)
    ) %>%
    mutate(
      age_mid = min(age_bl) + (age_bin - 0.5) * (max(age_bl) - min(age_bl)) / num_bins,
      path_mid = min(path_bl) + (path_bin - 0.5) * (max(path_bl) - min(path_bl)) / num_bins
    ) %>% 
    group_by(age_mid, path_mid) %>% 
    summarise(n = n(), .groups = "drop") -> x
  
  x %>%
    ggplot() + 
    geom_tile(aes(x=age_mid, y=path_mid, fill=n), color = "black") +  
    annotate("rect", xmin = 18.5, xmax = 88, ymin = 0, ymax = 0.35, fill = NA, color = "black") +
    annotate("rect", xmin = 19, xmax = 88.5, ymin = 0.015, ymax = 0.365, fill = NA, linetype = 2, color = "black") +
    annotate("rect", xmin = 20, xmax = 45, ymin = -0.025, ymax = 1, fill = NA, color = "black") +
    annotate("rect", xmin = 21, xmax = 46, ymin = -0.015, ymax = 1.01, fill = NA, linetype = 2, color = "black") +
    annotate("curve", x = 90, y = 0.20, xend = 90, yend = 0.365,
             arrow = arrow(length = unit(0.1, "inches"))) +
    annotate("curve", x = 35, y = 1.025, xend = 46, yend = 1.025,
             arrow = arrow(length = unit(0.1, "inches")), curvature = -0.3) +
    annotate("text", x = 40, y = 1.1, label = "1 Year", size = 6) +
    annotate("text", x = 95, y = 0.305, label = "0.015 Path", angle = -90, size = 6) +
    theme_gray(base_size = b_size) +
    scale_x_continuous("Age BL") + 
    scale_y_continuous("Pathology BL", breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) + 
    scale_fill_continuous("Count", breaks = c(2, 5, 8)) + 
    theme(#axis.text = element_text(size = 12),
      legend.position = "top",
      #legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-25,-10),
      legend.justification="right",
      legend.title = element_text(vjust = 0),
      legend.text.position = "top",
      legend.text = element_text(vjust = -1),
      plot.background = element_blank(),
      #legend.position.inside = c(0.95, 0.825),
      legend.background = element_blank(),
      #title = element_text(size = 12, face = "bold"),
      panel.border = element_rect(linewidth = 1, color = "black", fill = NA)) -> hist_map
  
  
  ######################
  # Age window analysis
  ######################
  
  long_bf_ %>% filter(fmri_bl) %>% 
    select(age_bl) -> age_df
  
  win_size <- 25
  frames <- data.frame(win_n = seq(min(age_df$age_bl) %>% round(), max(age_df$age_bl) %>% round() - win_size, by = 1)) %>%
    mutate(
      window_min = win_n,
      window_max = win_n + win_size,
    ) %>% 
    rowwise() %>%  
    mutate(sample_size = sum(age_df$age_bl >= window_min & age_df$age_bl <= window_max)) %>%
    mutate(mean_age = mean(age_df$age_bl[age_df$age_bl >= window_min & age_df$age_bl <= window_max])) %>%
    ungroup() %>% 
    mutate(window_range = paste0("[",window_min, ", ", window_max, "]"))
  
  
  if (run_windowing) {
    sliding_window_results <- data.frame()
    for (win_n in frames$win_n) {
      window_df <- long_bf_ %>% filter((age_bl >= win_n) , (age_bl <= win_size + win_n))
      ests_window <- nodal_lmm_ests(window_df, fc_measures_bf$affinity, roi_names = rois, model_formula = mod_formula)
      ests_window$window <- win_n
      sliding_window_results <- rbind(sliding_window_results, ests_window)
    }
    write_rds(sliding_window_results, file.path(processed_dir, "sliding_window_age_res.rds"))
  } else {
    sliding_window_results <- readRDS(file.path(processed_dir, "sliding_window_age_res.rds"))
  }
  
  grad_cor_age <- sliding_window_results %>% 
    filter(!(term %in% c("(Intercept)", "time", "rsqa__MeanFD", "sex"))) %>% 
    group_by(window, term, n) %>% 
    summarise(Gradient1 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient1), statistic), 
              Gradient3 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient3), statistic)) %>% 
    ungroup() %>% 
    inner_join(frames, join_by(window == win_n))
  
  ###########################
  # Pathology window analysis
  ###########################
  
  long_bf_ %>% filter(fmri_bl) %>% 
    select(path_bl) -> path_df
  #seq(min(path_df$path_bl) %>% round(), max(path_df$path_bl) %>% round() - win_size, by = 0.01)
  win_size <- 0.35
  steps <- seq(min(path_df$path_bl) %>% round(), (max(path_df$path_bl) %>% round() - win_size) + 0.01, by = 0.015)
  steps %>% length()
  frames_path <- data.frame(win_n = steps) %>%
    mutate(
      window_min = win_n,
      window_max = win_n + win_size,
    ) %>% 
    rowwise() %>%  # Ensures calculations are done row-by-row
    mutate(sample_size = sum(path_df$path_bl >= window_min & path_df$path_bl <= window_max)) %>%
    mutate(mean_path = mean(path_df$path_bl[path_df$path_bl >= window_min & path_df$path_bl <= window_max])) %>%
    ungroup() %>% 
    mutate(window_range = paste0("[",window_min, ", ", window_max, "]"))
  
  if (run_windowing) {
    sliding_window_results_path <- data.frame()
    for (win_n in frames_path$win_n) {
      window_df <- long_bf_ %>% filter((path_bl >= win_n) , (path_bl <= win_size + win_n)) 
      ests_window <- nodal_lmm_ests(window_df, fc_measures_bf$affinity, roi_names = rois, model_formula = mod_formula)
      ests_window$window <- win_n
      sliding_window_results_path <- rbind(sliding_window_results_path, ests_window)
    }
    write_rds(sliding_window_results_path, file.path(processed_dir, "sliding_window_path_res.rds"))
  } else {
    sliding_window_results_path <- readRDS(file.path(processed_dir, "sliding_window_path_res.rds"))
  }
  
  grad_cor_path <- sliding_window_results_path %>% 
    filter(!(term %in% c("(Intercept)", "time", "rsqa__MeanFD", "sex"))) %>% 
    group_by(window, term, n) %>% 
    summarise(Gradient1 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient1), statistic), 
              Gradient3 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient3), statistic)) %>% 
    inner_join(frames_path, join_by(window == win_n))
  
  ###########################
  # Plotting window analyses
  ###########################
  
  
  k = 5
  grad_cor_p_age <- grad_cor_age %>% 
    pivot_longer(starts_with("Grad"), names_to = "gradient", values_to = "grad_corr") %>% 
    filter((gradient == "Gradient3" & term == "age_bl") | (gradient == "Gradient1" & term == "path_bl") | (gradient == "Gradient1" & term == "pathΔ") ) %>% 
    mutate(term = case_when(
      term == "age_bl" ~ "Age~at~baseline",
      term == "path_bl" ~ "Pathology~at~baseline",
      term == "pathΔ" ~ "Δ*Pathology~(t[i]-t[0])",
    )) %>%
    ggplot(aes(
      window,
      grad_corr,
      color = gradient,
      fill = gradient
    )) +
    geom_point(aes(alpha = n), show.legend = FALSE) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = k, bs = "cs")) +
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5 )+
    facet_wrap(~term, nrow = 1, labeller = label_parsed) +
    #theme_bw(base_size = 14) +
    scale_x_continuous(labels = frames$window_range[seq(1, length(frames$win_n), by = 5)], 
                       breaks = frames$win_n[seq(1, length(frames$win_n), by = 5)],
                       guide = guide_axis(check.overlap = TRUE)) +
    labs(y = "Gradient corr (r)",
         x = "Baseline Age Window") +
    theme(
      legend.position = "",
      legend.title = element_blank()
    ) +
    geom_xsidecol(aes(y = sample_size, fill = NULL, color = NULL)) +
    scale_xsidey_continuous(breaks = c(0, 250))+
    #theme_ggside_void() +
    ggside(x.pos = "bottom") +
    theme(ggside.panel.scale = 0.2,
          axis.text.x = element_text(angle =-30, hjust =0, size = rel(0.8))) +
    scale_color_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_fill_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_alpha_continuous(range = c(0.3, 1))
  
  ####### PATHOLOGY
  
  
  grad_cor_p_path <- grad_cor_path %>% 
    pivot_longer(starts_with("Grad"), names_to = "gradient", values_to = "grad_corr") %>% 
    filter((gradient == "Gradient3" & term == "age_bl") | (gradient == "Gradient1" & term == "path_bl") | (gradient == "Gradient1" & term == "pathΔ") ) %>%
    mutate(term = case_when(
      term == "age_bl" ~ "Age~at~baseline",
      term == "path_bl" ~ "Pathology~at~baseline",
      term == "pathΔ" ~ "Δ*Pathology~(t[i]-t[0])",
    )) %>%
    ggplot(aes(
      window,
      grad_corr,
      color = gradient,
      fill = gradient
    )) +
    geom_point(aes(alpha = n), show.legend = FALSE) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs")) +
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5 )+
    facet_wrap(~term, nrow = 1, labeller = label_parsed) +
    #theme_bw(base_size = 14) +
    scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75) )+
    scale_x_continuous(labels = frames_path$window_range[seq(1, length(frames_path$win_n), by = 10)], 
                       breaks = frames_path$win_n[seq(1, length(frames_path$win_n), by = 10)],
                       guide = guide_axis(check.overlap = TRUE)) +
    labs(y = "Gradient corr (r)",
         x = "Baseline Pathology Window") +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.12, -0.45),
      legend.direction = "horizontal",
      legend.text = element_text(
        margin = margin(r = 20, l = 10,  unit = "pt")),
      legend.title = element_blank()
    ) +
    geom_xsidecol(aes(y = sample_size, fill = NULL, color = NULL)) +
    #theme_ggside_void() +
    ggside(x.pos = "bottom") +
    scale_xsidey_continuous(breaks = c(0, 250))+
    theme(ggside.panel.scale = 0.2,
          axis.text.x = element_text(angle =-30, hjust =0, size = rel(0.8))) +
    scale_color_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_fill_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_alpha_continuous(range = c(0.3, 1))
  
  window_plots <-  wrap_plots(list(grad_cor_p_age, grad_cor_p_path), nrow = 2)+ 
    plot_annotation(
      theme = theme(
        plot.background = element_rect(color = "black") 
      )
    ) 
  
  #############
  # Fonky way to get legend
  #############
  
  get_net_legend <- function(){
    #scale_factor <- 5
    x <- grad_df %>% filter(study == "biofinder") %>% 
      ggplot(aes(gradient1, gradient3, color = name)) +
      geom_point(alpha = 0.5) +
      labs(color = "Yeo Network") +
      guides(color = guide_legend(
        label.hjust=0,
        byrow = TRUE,
        nrow = 2, 
        reverse = TRUE,
        override.aes = list(size = rel(4))
      )) +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe) +
      theme(
        legend.position = "bottom",
        #legend.key.size = unit(0.0, "cm"),
        legend.key.spacing.x = unit(0, "cm"),
        legend.key.spacing.y = unit(-1.25, "cm"),
        legend.direction = "horizontal",
        #legend.title.position = "",
        legend.text.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.65), margin = margin(l = 5, r = 10, unit = "pt")),
        legend.background = element_blank()
      )
    
    leg <- ggpubr::get_legend(x)
    leg <- ggpubr::as_ggplot(leg)
    leg 
  }
  net_legend <- get_net_legend()
  
  ##########################
  # Put together all
  #########################
  
  p_bf_long_cp <- ggdraw() + 
    draw_plot(p_bf_long) + 
    draw_plot(net_legend, x = 0.15, y = 0.025, width = 0.2, height = 0.05) +
    draw_plot_label("A", size = 20) +
    draw_label(paste0("BioFINDER Longitudinal", " (N=", n_long, ")"), x = 0.065, y = 0.965, hjust = 0, size =  18) +
    draw_label("FCS ~ age_bl + path_bl + Δpath + \ntime + sex + motion + (1|sub)", x = 0.015, y = 0.89, hjust = 0, size =  12)
  
  window_plots_cp <- ggdraw() + 
    draw_plot(window_plots) + 
    draw_plot_label("C", size = 20)
  
  lmm_tab <- image_read("paper/figures/conceptual_plot/LMM_table.png")
  
  ggdraw() +
    draw_plot(hist_map, x = 0.01, y = 0.3, width = 0.98, height = 0.7) +
    draw_plot(p_bf_long[[5]] + labs(title = "ΔPath") + theme(plot.title = element_text(vjust = 0)),
              x = 0.39, y = 0.02,  width = 0.3, height = 0.3) +
    draw_plot(p_bf_long[[1]] + theme(plot.title = element_text(vjust = 0)),
              x = 0.69, y = 0.025, width = 0.29, height = 0.29) +
    draw_image(lmm_tab, x = 0.02, y = 0.00, width = 0.35, height = 0.25) +
    draw_plot_label("B", size = 20) + 
    annotate("curve",  x = 0.54, y = 0.09, xend = 0.83, yend = 0.09, linewidth = 1,
             arrow = arrow(length = unit(0.13, "inches"), ends = "both")) +
    annotate("segment",  x = 0.33, y = 0.425, xend = 0.225, yend = 0.22, linewidth = 1,
             arrow = arrow(length = unit(0.13, "inches"))) +
    annotate("label", label = "Parcel-wise LMM \n in age window", x = 0.25, y = 0.3, size = 6) +
    annotate("text", label = "Pearson Correlation", x = 0.7, y = 0.02, size = 6) +
    annotate("segment",  x = 0.38, y = 0.145, xend = 0.43, yend = 0.145, 
             arrow = arrow(length = unit(0.13, "inches")),
             linewidth = 1) +
    theme(plot.background = element_rect(color = "black", linewidth = 1)) -> slide_meth
  
  final_fig <- ggdraw() +
    draw_plot(p_bf_long_cp, x = 0, y = 0.51, width = 0.645, height = 0.49) +
    draw_plot(slide_meth, x = 0.655, y = 0.51, width = 0.345, height = 0.49) +
    draw_plot(window_plots_cp, x = 0, y = 0.0, width = 1, height = 0.5) +
    annotate("curve",  x = 0.89, y = 0.5105, xend = 0.706, yend = 0.4075,
             linewidth = 3,
             color = "white",
             curvature = 0.29) +
    annotate("curve",  x = 0.89, y = 0.5105, xend = 0.706, yend = 0.4075, linewidth = 1,
             arrow = arrow(length = unit(0.13, "inches")), 
             curvature = 0.29) 
  
  ###################
  # Supplementary fig
  ###################
  
  grad_cor_p_path_supp <- grad_cor_path %>% 
    pivot_longer(starts_with("Grad"), names_to = "gradient", values_to = "grad_corr") %>% 
    mutate(term = case_when(
      term == "age_bl" ~ "Age~at~baseline",
      term == "path_bl" ~ "Pathology~at~baseline",
      term == "pathΔ" ~ "Δ*Pathology~(t[i]-t[0])",
    )) %>%
    ggplot(aes(
      window,
      grad_corr,
      color = gradient,
      fill = gradient
    )) +
    geom_point(aes(alpha = n), show.legend = FALSE) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = k, bs = "cs")) +
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5 )+
    facet_wrap(~term, nrow = 1, labeller = label_parsed) +
    #theme_bw(base_size = 14) +
    scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75) )+
    scale_x_continuous(labels = frames_path$window_range[seq(1, length(frames_path$win_n), by = 10)], 
                       breaks = frames_path$win_n[seq(1, length(frames_path$win_n), by = 10)],
                       guide = guide_axis(check.overlap = TRUE)) +
    labs(y = "Gradient corr (r)",
         x = "Baseline Pathology Window") +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.12, -0.45),
      legend.direction = "horizontal",
      legend.text = element_text(
        margin = margin(r = 20, l = 10,  unit = "pt")),
      legend.title = element_blank()
    ) +
    geom_xsidecol(aes(y = sample_size/2, fill = NULL, color = NULL)) +
    #theme_ggside_void() +
    ggside(x.pos = "bottom") +
    scale_xsidey_continuous(breaks = c(0, 250))+
    theme(ggside.panel.scale = 0.2,
          axis.text.x = element_text(angle =-30, hjust =0, size = rel(0.8))) +
    scale_color_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_fill_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_alpha_continuous(range = c(0.3, 1))
  
  
  grad_cor_p_age_supp <- grad_cor_age %>% 
    pivot_longer(starts_with("Grad"), names_to = "gradient", values_to = "grad_corr") %>% 
    mutate(term = case_when(
      term == "age_bl" ~ "Age~at~baseline",
      term == "path_bl" ~ "Pathology~at~baseline",
      term == "pathΔ" ~ "Δ*Pathology~(t[i]-t[0])",
    )) %>%
    ggplot(aes(
      window,
      grad_corr,
      color = gradient,
      fill = gradient
    )) +
    geom_point(aes(alpha = n), show.legend = FALSE) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = k, bs = "cs")) +
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5 )+
    facet_wrap(~term, nrow = 1, labeller = label_parsed) +
    #theme_bw(base_size = 14) +
    scale_x_continuous(labels = frames$window_range[seq(1, length(frames$win_n), by = 5)], 
                       breaks = frames$win_n[seq(1, length(frames$win_n), by = 5)],
                       guide = guide_axis(check.overlap = TRUE)) +
    labs(y = "Gradient corr (r)",
         x = "Baseline Age Window") +
    theme(
      legend.position = "",
      legend.title = element_blank()
    ) +
    geom_xsidecol(aes(y = sample_size, fill = NULL, color = NULL)) +
    scale_xsidey_continuous(breaks = c(0, 250))+
    #theme_ggside_void() +
    ggside(x.pos = "bottom") +
    theme(ggside.panel.scale = 0.2,
          axis.text.x = element_text(angle =-30, hjust =0, size = rel(0.8))) +
    scale_color_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_fill_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_alpha_continuous(range = c(0.3, 1))
  
  window_plots_both_G <-  wrap_plots(list(grad_cor_p_age_supp, grad_cor_p_path_supp), nrow = 2)+ 
    plot_annotation(
      theme = theme(
        plot.background = element_rect(color = "black") 
      )
    ) 
  
  # Return figures:
  list(main_fig = final_fig, supp_fig = window_plots_both_G)
  
}