
plot_gams_v1 <- function(gam_predictions, grad_df, 
                         atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds"),
                         gradient_cols = data.frame(gradient1 = c("#3F596D", "#D38A4E"), 
                                                    gradient2 =  c("#4682B4", "#781286"),  
                                                    gradient3 =  c("#8A6081", "#738518")),
                         add_shade = FALSE,
                         shade_alpha = 0.1,
                         shade_size = 0.1,
                         grad_name = TRUE,
                         biofinder_data, b_size = 7,
                         rasterize_scatters = FALSE,
                         rasterize_brains = FALSE,
                         brain_raster_width_px = 1500,
                         brain_raster_dpi = 300,
                         scatter_raster_dev = "ragg",
                         scatter_raster_dpi = NULL,
                         point_size = 0.1,
                         point_alpha = 0.3,
                         axis_title_rel = 0.825,
                         axis_text_rel = 0.65,
                         gradient_title_rel = 0.7,
                         corr_label_size = 4.5,
                         pathology_label_size = 2.3,
                         pathology_legend_text_rel = 0.45,
                         pathology_legend_title_rel = 0.55,
                         legend_text_rel = 0.55,
                         legend_key_height = 0.55,
                         legend_key_width = 0.012,
                         net_legend_text_rel = 0.55,
                         net_legend_point_size = 2,
                         pathology_legend_width = 0.12,
                         pathology_legend_height = 0.06,
                         net_legend_width = 0.16,
                         net_legend_height = 0.07,
                         canvas_text_size = b_size * 0.75,
                         figure_label_size = b_size * 1.2,
                         struct_linewidth = 0.3,
                         plot_linewidth = 0.4,
                         spintest = TRUE,
                         perms = readRDS("data/atlas_data/permutations_1000_hungarian.rds"),
                         figure_pat = "paper/figures",
                         ...) {
  
  require(ggside)
  require(ggpmisc)
  require(sf)
  require(ggtext)
  if (rasterize_brains) {
    require(png)
    require(ragg)
  }
  if (rasterize_scatters) {
    require(ggrastr)
  }
  options("ggrastr.default.dpi" = scatter_raster_dpi)
  
  old <- theme_set(theme_bw(base_size = b_size))
  theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA),
               legend.background = element_rect(fill = "transparent", colour = NA),
               legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  make_raster_overlay <- function(plot, bbox,
                                  width_px = brain_raster_width_px,
                                  dpi = brain_raster_dpi,
                                  bg = "transparent") {
    ratio <- as.numeric((bbox["xmax"] - bbox["xmin"]) / (bbox["ymax"] - bbox["ymin"]))
    height_px <- max(1, round(width_px / ratio))
    tf <- tempfile(fileext = ".png")
    on.exit(unlink(tf), add = TRUE)
    
    ggsave(
      filename = tf,
      plot = plot,
      width = width_px / dpi,
      height = height_px / dpi,
      dpi = dpi,
      bg = bg,
      device = ragg::agg_png
    )
    
    annotation_raster(
      raster = as.raster(png::readPNG(tf)),
      xmin = bbox["xmin"],
      xmax = bbox["xmax"],
      ymin = bbox["ymin"],
      ymax = bbox["ymax"]
    )
  }
  
  make_brain_raster_layer <- function(atlas_sf,
                                      fill_var,
                                      fill_scale,
                                      shade_sf = NULL,
                                      shade_size = 0.01,
                                      shade_alpha = 0.01,
                                      linewidth = 0.1) {
    atlas_sf <- sf::st_as_sf(atlas_sf)
    bbox <- sf::st_bbox(atlas_sf)
    
    p_raster <- ggplot() +
      geom_sf(
        data = atlas_sf,
        aes(fill = .data[[fill_var]], geometry = geometry),
        linewidth = linewidth,
        show.legend = FALSE
      ) +
      (
        if (!is.null(shade_sf)) {
          geom_sf(data = shade_sf, size = shade_size, alpha = shade_alpha, show.legend = FALSE)
        } else {
          NULL
        }
      ) +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]),
        expand = FALSE
      ) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)
      ) +
      fill_scale
    
    make_raster_overlay(p_raster, bbox = bbox)
  }
  
  point_layer <- function(alpha = point_alpha, show.legend = FALSE) {
    layer <- geom_point(size = point_size, alpha = alpha, show.legend = show.legend)
    if (rasterize_scatters) {
      ggrastr::rasterise(layer, dev = scatter_raster_dev, dpi = scatter_raster_dpi)
    } else {
      layer
    }
  }
  
  compact_plot_theme <- theme(
    axis.title = element_text(size = rel(axis_title_rel)),
    axis.text = element_text(size = rel(axis_text_rel)),
    legend.text = element_text(size = rel(legend_text_rel)),
    legend.title = element_text(size = rel(legend_text_rel))
  )
  
  source_data <- list()
  
  legend_labs <-  c("Ab42/40", "Braak12", "Braak34", "Braak56")
  biof_path_scaled <- biofinder_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |> 
    select(pathology_ad, ab_ratio, 
           starts_with("braak"), starts_with("cho")) |> 
    mutate(across(where(is.numeric) & !pathology_ad, scale)) 
  ab_mean <- attr(biof_path_scaled$ab_ratio, c("scaled:center"))
  ab_sd <- attr(biof_path_scaled$ab_ratio, c("scaled:scale"))
  ab_pos <- (0.08 - ab_mean)/ab_sd
  
  biof_path_source <- biof_path_scaled  |> 
    pivot_longer(-pathology_ad, names_to = "scaled_pat_measures", values_to = "value") |>
    mutate(
      pathology_label = dplyr::recode(
        scaled_pat_measures,
        ab_ratio = legend_labs[1],
        braak12 = legend_labs[2],
        braak_12 = legend_labs[2],
        braak34 = legend_labs[3],
        braak_34 = legend_labs[3],
        braak56 = legend_labs[4],
        braak_56 = legend_labs[4],
        .default = scaled_pat_measures
      ),
      ab_positive_scaled_cutoff = ab_pos
    )
  
  source_data$pathology_plot <- biof_path_source
  
  gradient_source <- grad_df |>
    select(any_of(c("study", "region", "name", "gradient1", "gradient2", "gradient3"))) |>
    distinct()
  
  source_data$gradient_maps <- gradient_source
  
  pat_quartile_source <- gam_predictions$pat_derivs |>
    pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "fcs_slope") |>
    mutate(
      predictor = "pathology_ad",
      predictor_value = pathology_ad,
      predictor_quartile = cut(pathology_ad, 4)
    ) |>
    group_by(predictor, predictor_quartile, region) |>
    summarise(
      mean_predictor_value = mean(predictor_value, na.rm = TRUE),
      fcs_slope = mean(fcs_slope, na.rm = TRUE),
      .groups = "drop"
    ) |>
    inner_join(gradient_source, by = "region")
  
  age_quartile_source <- gam_predictions$age_derivs |>
    pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "fcs_slope") |>
    mutate(
      predictor = "age",
      predictor_value = age,
      predictor_quartile = cut(age, 4)
    ) |>
    group_by(predictor, predictor_quartile, region) |>
    summarise(
      mean_predictor_value = mean(predictor_value, na.rm = TRUE),
      fcs_slope = mean(fcs_slope, na.rm = TRUE),
      .groups = "drop"
    ) |>
    inner_join(gradient_source, by = "region")
  
  source_data$quartile_fcs_slopes <- bind_rows(
    pat_quartile_source,
    age_quartile_source
  )
  
  pat_pred_source <- gam_predictions$pat_pred |>
    pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "predicted_fc") |>
    inner_join(gradient_source |> select(region, any_of(c("study", "name", "gradient1", "gradient3"))), by = "region") |>
    mutate(
      predictor = "pathology_ad",
      predictor_value = pathology_ad,
      plotted_gradient = "gradient1",
      plotted_gradient_value = gradient1,
      grad_grp = cut(gradient1, breaks = quantile(gradient1, probs = seq(0, 1, length.out = 21)))
    )
  
  age_pred_source <- gam_predictions$age_pred |>
    pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "predicted_fc") |>
    inner_join(gradient_source |> select(region, any_of(c("study", "name", "gradient1", "gradient3"))), by = "region") |>
    mutate(
      predictor = "age",
      predictor_value = age,
      plotted_gradient = "gradient3",
      plotted_gradient_value = gradient3,
      grad_grp = cut(gradient3, breaks = quantile(gradient3, probs = seq(0, 1, length.out = 21)))
    )
  
  pred_source <- bind_rows(pat_pred_source, age_pred_source)
  
  pred_ventile_source <- pred_source |>
    group_by(predictor, predictor_value, plotted_gradient, grad_grp) |>
    summarise(
      mean_predicted_fc = mean(predicted_fc, na.rm = TRUE),
      mean_gradient_value = mean(plotted_gradient_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  pat_corr_source <- gam_predictions$pat_derivs |>
    pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "fcs_slope") |>
    inner_join(gradient_source |> select(region, any_of(c("gradient1", "gradient3"))), by = "region") |>
    group_by(pathology_ad) |>
    summarise(
      predictor = "pathology_ad",
      predictor_value = pathology_ad[1],
      plotted_gradient = "gradient1",
      gradient_correlation = cor(fcs_slope, gradient1),
      .groups = "drop"
    )
  
  age_corr_source <- gam_predictions$age_derivs |>
    pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "fcs_slope") |>
    inner_join(gradient_source |> select(region, any_of(c("gradient1", "gradient3"))), by = "region") |>
    group_by(age) |>
    summarise(
      predictor = "age",
      predictor_value = age[1],
      plotted_gradient = "gradient3",
      gradient_correlation = cor(fcs_slope, gradient3),
      .groups = "drop"
    )
  
  gradient_corr_source <- bind_rows(pat_corr_source, age_corr_source)
  
  source_data$predicted_fcs_and_gradient_correlations <- pred_source |>
    left_join(
      gradient_corr_source,
      by = c("predictor", "predictor_value", "plotted_gradient")
    )
  
  source_data$ventile_mean_predicted_fcs <- pred_ventile_source |>
    left_join(
      gradient_corr_source,
      by = c("predictor", "predictor_value", "plotted_gradient")
    )
  
  source_data$gradient_correlations <- gradient_corr_source |>
    select(
      predictor,
      predictor_value,
      plotted_gradient,
      gradient_correlation
    )
  
  biof_path_plot <- biof_path_source %>% 
    ggplot(aes(pathology_ad, value, color = scaled_pat_measures)) +
    geom_smooth(linewidth = plot_linewidth) +
    guides(color = guide_legend(
      ncol=1, byrow=TRUE,
      title.position="top", 
      title.hjust = 0)) +
    geom_hline(yintercept = ab_pos, linetype = 6, linewidth = struct_linewidth) +
    geom_text(inherit.aes = FALSE, 
              data = data.frame(x = 0.85, y = ab_pos),
              aes(x = x, y = y),
              label = "A*beta*'+'",
              parse = TRUE, nudge_y = 0.5, size = pathology_label_size) +
    labs(x = "Pathology score", y = "Scaled value")+
    ggsci::scale_color_nejm(name = "Pathology", labels = legend_labs) +
    compact_plot_theme +
    theme(legend.position = "right",
          legend.title.position = "top",
          legend.key.height = unit(0.35, "lines"),
          legend.key.width = unit(0.45, "lines"),
          legend.spacing.y = unit(0.02, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          #ggside.panel.scale = 0.2,
          legend.text = element_text(size = rel(pathology_legend_text_rel)),
          legend.title = element_text(size = rel(pathology_legend_title_rel)))
  # geom_xsidehistogram(aes(x = pathology_ad),
  #                     color = "black",
  #                     fill = "lightgray",
  #                     bins = 20,
  #                     data = biofinder_data %>% filter(fmri_bl, !is.na(age), !is.na(pathology_ad)),
  #                      inherit.aes = FALSE) +
  # scale_xsidey_continuous(breaks = c(0, 100, 200), position = "left") +
  # ggside(x.pos = "bottom")
  
  pathology_plot <- biof_path_plot +
    theme(legend.position = "")  
  
  pat_leg <- ggpubr::get_legend(biof_path_plot)
  pat_leg <- ggpubr::as_ggplot(pat_leg)
  
  
  gradient_plots <- list()
  i = 1
  grad_char <- c("gradient1", "gradient2", "gradient3")
  for (grad in grad_char) {
    atlas_plot_df <- grad_df %>% filter(study=="biofinder") %>% 
      #mutate(segregation = ifelse(segregation<0, 0, segregation)) %>% 
      right_join(atlas_geometry$atlas, by = "region") %>%
      sf::st_as_sf()
    
    bbox <- sf::st_bbox(atlas_plot_df)
    fill_limits <- range(atlas_plot_df[[grad]], na.rm = TRUE)
    gradient_fill_scale <- scale_fill_gradient2(
      low = gradient_cols[[grad]][1],
      mid = "white",
      high = gradient_cols[[grad]][2],
      limits = fill_limits
    )
    
    gradient_plots[[paste0(grad)]] <- ggplot() +
      (if (rasterize_brains)
        make_brain_raster_layer(
          atlas_sf = atlas_plot_df,
          fill_var = grad,
          fill_scale = gradient_fill_scale,
          shade_sf = if (add_shade) atlas_geometry$shade else NULL,
          shade_size = shade_size,
          shade_alpha = shade_alpha,
          linewidth = 0.1
        )
       else
         geom_sf(data = atlas_plot_df,
                 aes(fill = .data[[grad]], geometry = geometry),
                 linewidth= 0.1,
                 show.legend = FALSE))+
      # geom_sf(data = network_geometry %>% drop_na(),
      #         aes(
      #           #fill = region,
      #           color = name,
      #           geometry = geometry
      #         ), alpha = 0, linewidth = 0.5,
      #         show.legend = FALSE) +
      (if (!rasterize_brains && add_shade)
        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
       else NULL) +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]),
        expand = FALSE
      ) +
      theme_void(base_size = b_size)+
      (if (grad_name)
        labs(title = case_when(
          grad == "gradient1" ~ "Sens-Assoc",
          grad == "gradient2" ~ "Vis-Mot",
          grad == "gradient3" ~ "Rep-Exec",
          TRUE ~ grad
        ))  else NULL) +
      theme(legend.position = "",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(hjust = 0.5, size = rel(gradient_title_rel))
      ) +
      #guides(color = guide_legend(override.aes = list(size = 1))) +
      gradient_fill_scale
    i = i +1
  }
  
  
  brain_pat <- list()
  for (pat_grp in unique(cut(gam_predictions$pat_derivs$pathology_ad, 4))) {
    atlas_plot_df <- gam_predictions$pat_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(pathology_ad = cut(pathology_ad, 4)) %>% 
      filter(pathology_ad == pat_grp) %>% 
      group_by(region, pathology_ad) %>% 
      summarise(value=mean(value)) %>% 
      right_join(atlas_geometry$atlas, by = "region") %>%
      sf::st_as_sf()
    
    bbox <- sf::st_bbox(atlas_plot_df)
    fill_limits <- range(atlas_plot_df$value, na.rm = TRUE)
    brain_fill_scale <- scale_fill_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red"),
      limits = fill_limits
    )
    
    brain_pat[[pat_grp]] <- ggplot() +
      (if (rasterize_brains)
        make_brain_raster_layer(
          atlas_sf = atlas_plot_df,
          fill_var = "value",
          fill_scale = brain_fill_scale,
          shade_sf = if (add_shade) atlas_geometry$shade else NULL,
          shade_size = shade_size,
          shade_alpha = shade_alpha,
          linewidth = 0.1
        )
       else
         geom_sf(data = atlas_plot_df,
                 aes(fill = value, geometry = geometry),
                 linewidth= 0.1,
                 show.legend = FALSE))+
      (if (!rasterize_brains && add_shade)
        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
       else NULL) +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]),
        expand = FALSE
      ) +
      theme_void()+
      theme(legend.position = "",
            plot.margin = unit(c(0, 0, 0, 0), "npc"),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(color = "black", hjust = 0.5, size = rel(2)))+
      brain_fill_scale
  }
  
  scatter_pat <- list()
  for (pat_grp in unique(cut(gam_predictions$pat_derivs$pathology_ad, 4))) {
    gam_predictions$pat_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(pathology_ad = cut(pathology_ad, 4)) %>% 
      filter(pathology_ad == pat_grp) %>% 
      group_by(region, pathology_ad) %>% 
      summarise(value=mean(value)) %>% 
      inner_join(grad_df) -> plot_df
    
    x_range_df <- data.frame(x = seq(min(plot_df$gradient1), max(plot_df$gradient1), length.out = 100))
    colnames(x_range_df)[1] <- "gradient1"
    y_range_df <- data.frame(y = seq(min(plot_df$value), max(plot_df$value), length.out = 100))
    colnames(y_range_df)[1] <- "value"
    
    x_axis_colorbar <- gradient_cols[[1]]
    
    scatter_pat[[pat_grp]] <- plot_df %>% 
      ggplot(aes(x = value, y = gradient1,
                 color = name)) +
      point_layer(alpha = 0.1, show.legend = FALSE) +
      stat_poly_line(se = FALSE, color = "#323232", linewidth = struct_linewidth) +
      #xlim(-0.2, 0.05) +
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE, angle = 30), position = "bottom") +
      labs(
        #title = "Pathology AD {current_frame}",
        color = "Network",
        x = "FCS slopes averaged over pathology quartiles", 
        y = "") +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe())+
      scale_y_continuous(position = "right", limits = c(min(grad_df$gradient1), max(grad_df$gradient1)),
                         guide = guide_axis(check.overlap = TRUE)
                         ) +
      compact_plot_theme +
      theme(axis.title.y = element_blank()) +
      geom_xsidetile(data = y_range_df, aes(x = value, y = 0, fill = value), 
                     show.legend = FALSE, 
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red")
      ) +
      ggnewscale::new_scale_fill() +
      geom_ysidetile(data = x_range_df, aes(y = gradient1, x = 0, fill = gradient1), 
                     alpha = ifelse(pat_grp == unique(cut(gam_predictions$pat_derivs$pathology_ad, 4))[4], 1, 0),
                     show.legend = FALSE,
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = x_axis_colorbar[1],
        mid = "white",
        high = x_axis_colorbar[2] 
      ) +
      theme_ggside_void() +
      ggside(x.pos = "bottom", y.pos = "right") +
      theme(
        ggside.panel.scale = 0.05,
        axis.title.x = element_text(size = rel(axis_title_rel)),
        axis.text = element_text(size = rel(axis_text_rel))
      )
    
    if (spintest) {
      
      r <- cor(plot_df$gradient1, plot_df$value, method = "pearson")
      p_val <- perm_sphere_p(plot_df$gradient1, plot_df$value, perm.id = perms, corr.type='pearson')
      
      p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)

      label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)


      # p_lab <- scales::label_pvalue(
      #   accuracy = 0.001,
      #   prefix = c(
      #     "<i>p<sub>spin</sub></i> &lt; ",
      #     "<i>p<sub>spin</sub></i> = ",
      #     "<i>p<sub>spin</sub></i> &gt; "
      #   )
      # )(p_val)
      # 
      # label_corr <- paste0("<i>r</i> = ", round(r, 2), ", ", p_lab)
      
      scatter_pat[[pat_grp]] <- scatter_pat[[pat_grp]] + ggpp::annotate(geom = "text_npc", label = label_corr,
                                                                        size = corr_label_size,
                                                                        npcx = "left", npcy = "top",
                                                                        family = "sans",
                                                                        parse = TRUE)
      
        # labs(subtitle = label_corr) +
        # theme(
        #   plot.subtitle = element_markdown(
        #     size = corr_label_size,
        #     hjust = 0.95,
        #     margin = margin(b = 3)
        #   )
        # )
    } else {
      scatter_pat[[pat_grp]] <- scatter_pat[[pat_grp]] +       
        stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                       str_remove(after_stat(rr.confint.label), "95% CI "),
                                       sep = "*\" \"*")),
                     parse = TRUE, color = "#323232", label.x = "left", label.y = "top", size = corr_label_size) 
    }
    
  }
  
  
  brain_age <- list()
  for (age_grp in unique(cut(gam_predictions$age_derivs$age, 4))) {
    atlas_plot_df <- gam_predictions$age_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(age = cut(age, 4)) %>% 
      filter(age == age_grp) %>% 
      group_by(region, age) %>% 
      summarise(value=mean(value)) %>% 
      right_join(atlas_geometry$atlas, by = "region") %>%
      sf::st_as_sf()
    
    bbox <- sf::st_bbox(atlas_plot_df)
    fill_limits <- range(atlas_plot_df$value, na.rm = TRUE)
    brain_fill_scale <- scale_fill_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red"),
      limits = fill_limits
    )
    
    brain_age[[age_grp]] <- ggplot() +
      (if (rasterize_brains)
        make_brain_raster_layer(
          atlas_sf = atlas_plot_df,
          fill_var = "value",
          fill_scale = brain_fill_scale,
          shade_sf = if (add_shade) atlas_geometry$shade else NULL,
          shade_size = shade_size,
          shade_alpha = shade_alpha,
          linewidth = 0.1
        )
       else
         geom_sf(data = atlas_plot_df,
                 aes(fill = value, geometry = geometry),
                 linewidth= 0.1,
                 show.legend = FALSE))+
      (if (!rasterize_brains && add_shade)
        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
       else NULL) +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]),
        expand = FALSE
      ) +
      theme_void()+
      #labs(title = "Pathology AD {current_frame}")+
      theme(legend.position = "",
            plot.margin = unit(c(0, 0, 0, 0), "npc"),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(color = "black", hjust = 0.5, size = rel(2)))+
      brain_fill_scale
  }
  
  scatter_age <- list()
  for (age_grp in unique(cut(gam_predictions$age_derivs$age, 4))) {
    gam_predictions$age_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(age = cut(age, 4)) %>% 
      filter(age == age_grp) %>% 
      group_by(region, age) %>% 
      summarise(value=mean(value)) %>% 
      inner_join(grad_df) -> plot_df
    
    x_range_df <- data.frame(x = seq(min(plot_df$gradient3), max(plot_df$gradient3), length.out = 100))
    colnames(x_range_df)[1] <- "gradient3"
    y_range_df <- data.frame(y = seq(min(plot_df$value), max(plot_df$value), length.out = 100))
    colnames(y_range_df)[1] <- "value"
    
    x_axis_colorbar <- gradient_cols[[3]]
    
    scatter_age[[age_grp]] <-  plot_df %>% 
      ggplot(aes(x = value, y = gradient3,
                 color = name)) +
      point_layer(alpha = 0.1, show.legend = FALSE) +
      stat_poly_line(se = FALSE, color = "#323232", linewidth = struct_linewidth) +
      #ylim(min(grad_df$gradient3), max(grad_df$gradient3)) +
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE, angle = 30)#, position = "top"
      ) +
      scale_y_continuous(limits = c(min(grad_df$gradient3), max(grad_df$gradient3)),
                         guide = guide_axis(check.overlap = TRUE)
      ) +
      labs(
        color = "Network",
        x = "FCS slopes averaged over age quartiles", 
        y = "") +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
      compact_plot_theme +
      theme(axis.title.y = element_blank()) +
      geom_xsidetile(data = y_range_df, aes(x = value, y = 0, fill = value), 
                     show.legend = FALSE, 
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red")
      ) +
      ggnewscale::new_scale_fill() +
      geom_ysidetile(data = x_range_df, aes(y = gradient3, x = 0, fill = gradient3), 
                     alpha = ifelse(age_grp == unique(cut(gam_predictions$age_derivs$age, 4))[1], 1, 0),
                     show.legend = FALSE,
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = x_axis_colorbar[1],
        mid = "white",
        high = x_axis_colorbar[2] 
      ) +
      theme_ggside_void() +
      ggside(x.pos = "bottom", y.pos = "left") +
      theme(
        ggside.panel.scale = 0.05,
        axis.title.x = element_text(size = rel(axis_title_rel)),
        axis.text = element_text(size = rel(axis_text_rel))
      )
    
    
    if (spintest) {
      
      r <- cor(plot_df$gradient3, plot_df$value, method = "pearson")
      p_val <- perm_sphere_p(plot_df$gradient3, plot_df$value, perm.id = perms, corr.type='pearson')
      p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
      label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
      
      # p_lab <- scales::label_pvalue(
      #   accuracy = 0.001,
      #   prefix = c(
      #     "<i>p<sub>spin</sub></i> &lt; ",
      #     "<i>p<sub>spin</sub></i> = ",
      #     "<i>p<sub>spin</sub></i> &gt; "
      #   )
      # )(p_val)
      # 
      # label_corr <- paste0("<i>r</i> = ", round(r, 2), ", ", p_lab)

      
      scatter_age[[age_grp]] <- scatter_age[[age_grp]] + ggpp::annotate(geom = "text_npc", label = label_corr,
                                                                        size = corr_label_size,
                                                                        npcx = "left", npcy = "top",
                                                                        family = "sans",
                                                                        parse = TRUE)
      
      
        # labs(subtitle = label_corr) +
        # theme(
        #   plot.subtitle = element_markdown(
        #     size = corr_label_size,
        #     hjust = 0.95,
        #     margin = margin(b = 3)
        #   ))
        # 
        

    } else {
      scatter_age[[age_grp]] <- scatter_age[[age_grp]] + 
        stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                       str_remove(after_stat(rr.confint.label), "95% CI "),
                                       sep = "*\" \"*")),
                     parse = TRUE, color = "#323232", label.x = "left", label.y = "top", size = corr_label_size) 
    }
    
  }
  
  
  quantile_trajectories <- gam_predictions$pat_pred %>% pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "predicted_fc") %>% 
    inner_join(grad_df %>% select(region, gradient1)) %>% 
    mutate(grad_grp = cut(gradient1, breaks = quantile(gradient1, probs = seq(0, 1, length.out = 21)))) %>% 
    group_by(pathology_ad, grad_grp) %>% 
    summarise(mean_pred_fc = mean(predicted_fc),
              mean_grad_value = mean(gradient1)) %>% 
    ggplot(aes(pathology_ad, mean_pred_fc, group = mean_grad_value, color = mean_grad_value)) +
    geom_line(linewidth = plot_linewidth) +
    labs(x = "Pathology score", y = "Predicted FCS", color = "Mean SA\nin ventile") +
    #guides(colour = guide_colorbar(title.position="top", title.hjust = 0.0))+
    scale_color_gradient2(
      high = gradient_cols[[1]][2], mid = "white", low = gradient_cols[[1]][1]) +
    compact_plot_theme +
    theme(
      legend.key.height = unit(legend_key_height, "lines"),
      legend.key.width = unit(legend_key_width, "npc"),
      legend.box.spacing = unit(0.0025, "npc"),
      legend.title = element_blank(),
      legend.text = element_text(margin = margin(l = 0.8), size = rel(legend_text_rel)),
      legend.ticks = element_blank()
    ) 
  
  
  r2_pat <- gam_predictions$pat_derivs %>% pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>%
    inner_join(grad_df) %>%
    group_by(pathology_ad) %>%
    summarise(grad_R2 = cor(value, gradient1)) %>%
    ggplot(aes(pathology_ad, grad_R2)) +
    geom_line(linewidth = plot_linewidth) +
    labs(y = "Corr (r)", x = "Pathology score") +
    compact_plot_theme +
    theme(
      ggside.panel.scale = 0.2
    ) 
  # +
  #   geom_xsidehistogram(aes(x = pathology_ad), data = biofinder_data, inherit.aes = FALSE) +
  #   ggside(x.pos = "top")
  
  
  
  for(x_int in seq_range(gam_predictions$pat_pred$pathology_ad, 5)){
    quantile_trajectories <- quantile_trajectories +
      geom_vline(xintercept = x_int, linetype = "dashed",
                 linewidth = struct_linewidth)
    r2_pat <- r2_pat +
      geom_vline(xintercept = x_int, linetype = "dashed",
                 linewidth = struct_linewidth)
    pathology_plot <- pathology_plot +
      geom_vline(xintercept = x_int, linetype = "dashed",
                 linewidth = struct_linewidth)
  }
  
  age_range <- range(gam_predictions$age_pred$age)
  
  quantile_trajectories_term2 <- gam_predictions$age_pred %>% 
    pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "predicted_fc") %>% 
    inner_join(grad_df %>% select(region, gradient3)) %>% 
    mutate(grad_grp = cut(gradient3, breaks = quantile(gradient3, probs = seq(0, 1, length.out = 21)))) %>% 
    group_by(age, grad_grp) %>% 
    summarise(mean_pred_fc = mean(predicted_fc),
              mean_grad_value = mean(gradient3)) %>% 
    ggplot(aes(age, mean_pred_fc, group = mean_grad_value, color = mean_grad_value)) +
    geom_line(linewidth = plot_linewidth) +
    labs(y = "Predicted FCS", x = "Age", color = "Mean RE\nin ventile") +
    scale_color_gradient2(
      high = gradient_cols[[3]][2], mid = "white", low = gradient_cols[[3]][1]) +
    #guides(colour = guide_colorbar(title.position="top", title.hjust = 0.0))+
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    compact_plot_theme +
    theme(
      legend.key.height = unit(legend_key_height, "lines"),
      legend.key.width = unit(legend_key_width, "npc"),
      legend.box.spacing = unit(0.005, "npc"),
      legend.title = element_blank(),
      legend.text = element_text(margin = margin(l = 0.8), size = rel(legend_text_rel)),
      legend.ticks = element_blank()
    ) 
  
  r2_term2 <- gam_predictions$age_derivs %>% pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>%
    inner_join(grad_df ) %>%
    group_by(age) %>%
    summarise(grad_R2 = cor(value, gradient3)) %>%
    ggplot(aes(age, grad_R2)) +
    geom_line(linewidth = plot_linewidth) +
    compact_plot_theme +
    theme(ggside.panel.scale = 0.2) +
    labs(y = "Corr (r)", x = "Age") 
  # geom_xsidehistogram(aes(x = age), data = biofinder_data |> 
  #                       filter(age > age_range[1], age < age_range[2]), 
  #                     inherit.aes = FALSE) +
  # ggside(x.pos = "bottom")
  
  
  for(x_int in seq_range(gam_predictions$age_pred$age, 5)){
    quantile_trajectories_term2 <- quantile_trajectories_term2 +
      geom_vline(xintercept = x_int, linetype = "dashed",
                 linewidth = struct_linewidth)
    r2_term2 <- r2_term2 +
      geom_vline(xintercept = x_int, linetype = "dashed",
                 linewidth = struct_linewidth)
  }
  
  
  
  a <- wrap_plots(c(scatter_pat,
                    list(gradient_plots[[1]]),
                    list(plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer()),
                    brain_pat, 
                    list(plot_spacer())
  ), 
  nrow = 3) +
    plot_layout(
      guides = "collect",           
      axis_titles = "collect",      
      axes = "collect",
      heights = c(0.5, -0.15, 0.5)          
    )
  
  a <- ggdraw() +
    draw_plot(a) +
    draw_figure_label("A", size = figure_label_size)
  
  
  b <- 
    wrap_plots(list(pathology_plot + labs(tag = "B.1"), plot_spacer(),
                    r2_pat + labs(tag = "B.2"), r2_term2,
                    quantile_trajectories+ labs(tag = "B.3"), quantile_trajectories_term2),
               ncol = 2, byrow = TRUE) +
    plot_layout(
      axis_titles = "collect", axes = "collect"
      #guides = "collect" 
    )  & 
    theme(plot.tag.position  = c(0.23, 0.85),
          plot.tag = element_text(size = rel(0.6), hjust = 0, vjust = 0),
          axis.title = element_text(size = rel(axis_title_rel)),
          axis.text = element_text(size = rel(axis_text_rel)),
          legend.text = element_text(size = rel(legend_text_rel)),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA)
          #legend.box = "horizontal"
    ) 
  
  b <- ggdraw() +
    draw_plot(b) +
    draw_figure_label("B", size = figure_label_size)
  
  c <- wrap_plots(c(list(plot_spacer()),
                    brain_age, 
                    list(plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer()),
                    list(gradient_plots[[3]]),
                    scatter_age), nrow = 3) +
    plot_layout(
      guides = "collect",           
      axis_titles = "collect",      
      axes = "collect",
      heights = c(0.5, -0.2, 0.5)
    )
  
  c <- ggdraw() +
    draw_plot(c) +
    draw_figure_label("C", size = figure_label_size)
  
  get_net_legend <- function(ba_size = 11, point_size = net_legend_point_size){
    #scale_factor <- 5
    x <- grad_df %>% filter(study == "biofinder") |> 
      ggplot(aes(gradient1, gradient3, color = name)) +
      geom_point(alpha = 0.5) +
      theme_bw(base_size = ba_size) +
      guides(color = guide_legend(
        label.hjust=0,
        byrow = TRUE,
        nrow = 2, 
        override.aes = list(size = point_size)
      )) +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe) +
      theme(
        legend.position = "bottom",
        #legend.key.size = unit(0.0, "cm"),
        legend.key.spacing.x = unit(-0, "cm"),
        legend.key.spacing.y = unit(-0.85, "cm"),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.direction = "horizontal",
        #legend.title.position = "",
        legend.text.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = rel(net_legend_text_rel), margin = margin(l = -5, r = 0, unit = "pt")),
        legend.background = element_blank()
      )
    
    leg <- ggpubr::get_legend(x)
    leg <- ggpubr::as_ggplot(leg)
    leg 
  }
  
  net_legend <- get_net_legend(ba_size = b_size)
  
  age_hist <- biofinder_data |> filter(age > age_range[1], age < age_range[2]) |>  
    ggplot(aes(x = age)) +
    geom_histogram(bins = 20, color = "black", fill = "lightgray", linewidth = 0.05) +
    labs(y = "", x = "") +
    scale_y_continuous(breaks = c(0, 100), position = "right") +
    compact_plot_theme +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
    )
  
  path_hist <- biofinder_data  |>
    ggplot(aes(x = pathology_ad)) +
    geom_histogram(bins = 20, color = "black", fill = "lightgray", linewidth = 0.05) +
    labs(y = "", x = "") +
    scale_y_continuous(breaks = c(0, 200), position = "left") +
    compact_plot_theme +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
    )
  
  # 
  nonlin_p <-   ggdraw() +
    draw_plot(a, x = 0, y = 0.7, width = 1, height = 0.30) +
    draw_plot(b, x = 0, y = 0.30, width = 1, height = 0.40) +
    draw_plot(age_hist, x = 0.565, y = 0.5375, width = 0.412, height = 0.08) +
    draw_plot(path_hist, x = 0.01675, y = 0.645, width = 0.414, height = 0.08) +
    draw_plot(pat_leg, y = 0.6, x = 0.43, width = pathology_legend_width, height = pathology_legend_height) +
    draw_plot(c, x = 0, y = 0.0, width = 1, height = 0.30) +
    draw_plot(net_legend, x = 0.72, y = 0.63, width = net_legend_width, height = net_legend_height)+
    draw_text(text = c("Mean SA\nin ventile", 
                       "Mean RE\nin ventile"), x = c(0.475, 0.95), y = c(0.48, 0.48), size = canvas_text_size) +
    draw_line(
      x = c(0.12, 0.01, 0.01),
      y = c(0.68, 0.725, 0.81), 
      linetype = 2,
      linewidth = struct_linewidth
    ) +  
    draw_line(
      x = c(0.12 + 0.072, 0.01 + 0.19, 0.01 + 0.19),
      y = c(0.68, 0.725, 0.81), 
      linetype = 2,
      linewidth = struct_linewidth
    ) +
    draw_line(
      x = c(0.12 + 0.072*2, 0.01 + 0.19*2, 0.01 + 0.19*2),
      y = c(0.68, 0.725, 0.81), 
      linetype = 2,
      linewidth = struct_linewidth
    ) +
    draw_line(
      x = c(0.12 + 0.072*3, 0.01 + 0.19*3, 0.01 + 0.19*3),
      y = c(0.68, 0.725, 0.81), 
      linetype = 2,
      linewidth = struct_linewidth
    ) +
    draw_line(
      x = c(0.12 + 0.072*4, 0.01 + 0.19*4, 0.01 + 0.19*4),
      y = c(0.68, 0.725, 0.81), 
      linetype = 2,
      linewidth = struct_linewidth
    ) + # NEDRE DEL
    draw_line(
      x = c(0.515 + 0.072, 0.04 + 0.19, 0.04 + 0.19),
      y = c(0.368, 0.285, 0.19), 
      linetype = 2,
      linewidth = struct_linewidth
    ) +
    draw_line(
      x = c(0.515 + 0.072*2, 0.04 + 0.19*2, 0.04 + 0.19*2),
      y = c(0.368, 0.285, 0.19), 
      linetype = 2,
      linewidth = struct_linewidth
    ) +
    draw_line(
      x = c(0.515 + 0.072*3, 0.04 + 0.19*3, 0.04 + 0.19*3),
      y = c(0.368, 0.285, 0.19), 
      linetype = 2,
      linewidth = struct_linewidth
    ) +
    draw_line(
      x = c(0.515 + 0.072*4, 0.04 + 0.19*4, 0.04 + 0.19*4),
      y = c(0.368, 0.285, 0.19), 
      linetype = 2,
      linewidth = struct_linewidth
    )+
    draw_line(
      x = c(0.515 + 0.072*5, 0.04 + 0.19*5, 0.04 + 0.19*5),
      y = c(0.368, 0.285, 0.19), 
      linetype = 2,
      linewidth = struct_linewidth
    )
  list(
    plot = nonlin_p,
    source_data = source_data
  )
}

plot_gams_v1_legacy <- function(...) {
  args <- list(...)
  args$rasterize_scatters <- FALSE
  args$rasterize_brains <- FALSE
  if (is.null(args$axis_title_rel)) args$axis_title_rel <- 1
  if (is.null(args$axis_text_rel)) args$axis_text_rel <- 1
  if (is.null(args$gradient_title_rel)) args$gradient_title_rel <- 1
  if (is.null(args$corr_label_size)) args$corr_label_size <- 3.6
  if (is.null(args$pathology_label_size)) args$pathology_label_size <- 5
  if (is.null(args$pathology_legend_text_rel)) args$pathology_legend_text_rel <- 0.6
  if (is.null(args$pathology_legend_title_rel)) args$pathology_legend_title_rel <- 0.8
  if (is.null(args$legend_text_rel)) args$legend_text_rel <- 0.7
  if (is.null(args$net_legend_text_rel)) args$net_legend_text_rel <- 0.8
  if (is.null(args$net_legend_point_size)) args$net_legend_point_size <- 3
  if (is.null(args$pathology_legend_width)) args$pathology_legend_width <- 0.2
  if (is.null(args$pathology_legend_height)) args$pathology_legend_height <- 0.1
  if (is.null(args$net_legend_width)) args$net_legend_width <- 0.2
  if (is.null(args$net_legend_height)) args$net_legend_height <- 0.1
  if (is.null(args$canvas_text_size)) args$canvas_text_size <- 11
  if (is.null(args$figure_label_size)) args$figure_label_size <- 14
  out <- do.call(plot_gams_v1, args)
  if (is.list(out) && all(c("plot", "source_data") %in% names(out))) {
    out$plot
  } else {
    out
  }
}

plot_gams_legacy <- plot_gams_v1_legacy
