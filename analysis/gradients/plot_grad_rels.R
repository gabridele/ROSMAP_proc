
get_net_legend <- function(b_size){
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
      override.aes = list(size = rel(0.85))
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
      legend.text = element_text(size = rel(0.65), margin = margin(l = -4, r = 0, unit = "pt")),
      legend.background = element_blank()
    )

  leg <- ggpubr::get_legend(x)
  leg <- ggpubr::as_ggplot(leg)
  leg
}


plot_gradient_relationships <- function(subject_data,
                                        gradient_data,
                                        list_of_parcel_data,
                                        t_mat = NULL,
                                        rec_miss_contr = FALSE,
                                        atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds"),
                                        add_shade = FALSE,
                                        shade_alpha = 0.01,
                                        shade_size = 0.01,
                                        point_alpha = 0.2,
                                        point_size = 1,
                                        vect = FALSE,
                                        base_size_ = 11,
                                        grad_name = TRUE,
                                        gradients = c(1, 3),
                                        gradient_colors = NULL,
                                        empty_row_height = 0,
                                        padding = 0,
                                        r2_size = rel(4.6),
                                        r_spin_size = 0.75,
                                        spintest = TRUE,
                                        rectangle = FALSE,
                                        rasterize = FALSE,
                                        rasterize_scatters = FALSE,
                                        ggrastr_dev = "ragg",
                                        ggrastr_dpi = NULL,
                                        raster_scale_atlas = 1,
                                        raster_scale_shade = 1,
                                        l_width = 0.05,
                                        perms = readRDS("data/atlas_data/permutations_1000_hungarian.rds"),
                                        id_var = "image_file",
                                        mod_formula = formula(paste0("FC ~ age")),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        filter_criteria = quo(),
                                        plt_title = NULL,
                                        plt_title_parse = FALSE,
                                        plt_title_hjust = 0,
                                        plt_subtitle_hjust = 0,
                                        title_lmargin = 0,
                                        plt_title_size = rel(1),
                                        plt_subtitle_size = rel(0.9),
                                        brain_title_size = rel(1),
                                        axis_text_size = rel(1),
                                        axis_title_size = rel(1),
                                        brain_title_vjust = 0,
                                        brain_subtitle_vjust = 0,
                                        scatter_title_vjust = 0,
                                        plt_subtitle = FALSE,
                                        subtit_lookup = c(pathology_ad = "AD pathology", rsqa__MeanFD = "motion"),
                                        label_gradient_score = "Gradient Score",
                                        group_n_title = FALSE,
                                        brain_names = NULL,
                                        brain_names_parse = FALSE,
                                        tag_sep = "",
                                        tag_prefix = "",
                                        layout_construction = "horizontal",
                                        right_term_side = FALSE,
                                        include_gradient_plots = TRUE,
                                        gray_out = FALSE,
                                        side_color_bar = TRUE,
                                        plot_spacing = 0.3,
                                        show_networks = FALSE,
                                        network_geometry = NULL,
                                        plot_net_legend = FALSE,
                                        net_legend_x = 0.075,
                                        net_legend_y = 0.01,
                                        cache_runs = FALSE,
                                        longitudinal = FALSE,
                                        logistic_fit = FALSE,
                                        scale_fc = FALSE,
                                        sub_id = "sid",
                                        longitudinal_formula = formula(paste0("FC ~ age + (1 | ", sub_id, ")"))
){

  source("src/util.R")
  require(tidyverse)
  require(scales)
  require(patchwork)
  require(ggpmisc)
  require(sf)
  require(ggtext)
  require(ggrastr)

  options("ggrastr.default.dpi" = ggrastr_dpi)

  parse_plot_label <- function(label, parse = FALSE, label_name = "label") {
    if (!isTRUE(parse) || is.null(label) || inherits(label, "waiver")) {
      return(label)
    }

    if (length(label) != 1) {
      stop(label_name, " must have length 1 when parsed.", call. = FALSE)
    }

    tryCatch(
      parse(text = label)[[1]],
      error = function(e) {
        stop("Could not parse ", label_name, ": ", conditionMessage(e), call. = FALSE)
      }
    )
  }

  if (!is.logical(plt_title_parse) || length(plt_title_parse) != 1 || is.na(plt_title_parse)) {
    stop("plt_title_parse must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(rasterize_scatters) || length(rasterize_scatters) != 1 || is.na(rasterize_scatters)) {
    stop("rasterize_scatters must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(brain_names)) {
    if (!is.logical(brain_names_parse) || any(is.na(brain_names_parse))) {
      stop("brain_names_parse must be a logical vector without NA values.", call. = FALSE)
    }

    if (length(brain_names_parse) == 1) {
      brain_names_parse <- rep(brain_names_parse, length(brain_names))
    } else if (length(brain_names_parse) != length(brain_names)) {
      stop("brain_names_parse must have length 1 or the same length as brain_names.", call. = FALSE)
    }
  }


  old <- theme_set(theme_bw(base_size = base_size_))
  theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA),
               legend.background = element_rect(fill = "transparent", colour = NA),
               legend.box.background = element_rect(fill = "transparent", colour = NA))

  net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                          col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), #"#DCF8A4"
                          label = c(1:7))

  grad_char <- paste0("gradient", gradients)
  analysis_name <- names(list_of_parcel_data)

  if (is.null(gradient_colors)) {
    gradient_colors <- matrix(rep(c("#3A3A98", "#832424"), length(grad_char) ), ncol =length(grad_char))
    colnames(gradient_colors) <- grad_char
    gradient_colors <- as.data.frame(gradient_colors)
  }


  if (!longitudinal) {

    if (cache_runs) {
      prev_mod_formula <- function() {
        tryCatch(
          {
            list(suppressWarnings(readRDS("analysis_cache/model_formula.rds")),
                 readRDS("analysis_cache/analysis_name.rds"),
                 readRDS("analysis_cache/filter_crit.rds"))
          },
          error = function(cond) {
            message(paste("No cache or file exists"))
            message("Here's the original error message:")
            message(conditionMessage(cond))
            # Choose a return value in case of error
            NA
          })
      }
      prev_values <- prev_mod_formula()
      dir.create(file.path("analysis_cache"), showWarnings = FALSE)
      write_rds(mod_formula, "analysis_cache/model_formula.rds")
      write_rds(analysis_name, "analysis_cache/analysis_name.rds")
      write_rds(as_label(filter_criteria), "analysis_cache/filter_crit.rds")

      if (
        identical(mod_formula, prev_values[[1]]) & identical(analysis_name, prev_values[[2]]) & identical(as_label(filter_criteria), prev_values[[3]])
      ){
        list_of_fits <- read_rds("analysis_cache/list_of_fits.rds")
      } else {

        list_of_fits <- list()
        for (analysis in analysis_name) {
          print(paste0("Running linear models for ", analysis))
          list_of_fits[[analysis]] <-
            nodal_regression_fits(
              subject_data %>% filter(!!filter_criteria),
              list_of_parcel_data[[analysis]],
              vectorised = vect,
              roi_names = rois,
              id_var = id_var,
              logistic = logistic_fit,
              scale_fc = scale_fc,
              model_formula = mod_formula)
        }
        write_rds(list_of_fits, "analysis_cache/list_of_fits.rds")
      }

    } else {
      list_of_fits <- list()
      for (analysis in analysis_name) {
        print(paste0("Running linear models for ", analysis))


        if (!is.null(t_mat)) {
          list_of_fits[[analysis]] <-
            nodal_regression_fits_roiwise_pred(
              subject_data %>% filter(!!filter_criteria),
              list_of_parcel_data[[analysis]],
              tau_matrix = t_mat,
              vectorised = vect,
              roi_names = rois,
              id_var = id_var,
              model_formula = mod_formula)

        } else {
          list_of_fits[[analysis]] <-
            nodal_regression_fits(
              subject_data %>% filter(!!filter_criteria),
              list_of_parcel_data[[analysis]],
              vectorised = vect,
              roi_names = rois,
              id_var = id_var,
              model_formula = mod_formula)
        }

      }
    }

    list_of_ests <- list()
    print("Getting model estimates")
    for (analysis in names(list_of_fits)){
      list_of_ests[[analysis]] <- get_nodal_ests(list_of_fits[[analysis]], vectorised = vect,
                                                 recover_missing_contr = rec_miss_contr,
                                                 mc = TRUE) %>%
        select(term, region, statistic, n, model_formula)
    }

  }

  if (longitudinal) {
    list_of_ests <- list()
    print("fitting longitudinal")
    for (analysis in analysis_name){
      print(analysis)
      list_of_ests[[analysis]] <- nodal_lmm_ests(
        subject_data %>% filter(!!filter_criteria),
        list_of_parcel_data[[analysis]],
        roi_names = rois,
        id_var = id_var,
        subject_id = sub_id,
        model_formula = longitudinal_formula
      )
    }
  }


  ests <- list_of_ests[[1]] %>% rename_with(~ names(list_of_ests[1]), statistic)
  if (length(list_of_ests) > 1) {
    for(i in 2:length(list_of_ests)) {
      ests <- inner_join(ests, list_of_ests[[i]] %>% rename_with(~ names(list_of_ests[i]), statistic), by = c("term", "region"))
    }
  }

  if (group_n_title) {
    group_sizes <- get_factor_group_sizes(mod_formula, data = subject_data)
  }

  # shade_layer <- rasterise(
  #   geom_raster(
  #     data = shade_df,
  #     aes(x = x, y = y),
  #     inherit.aes = FALSE,
  #     interpolate = TRUE
  #   ),
  #   dev = ggrastr_dev,
  #   scale = raster_scale_shade
  # )

  # shade_layer <- rasterise(
  #   geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha),
  #   dev = ggrastr_dev,
  #   scale = raster_scale_shade
  # )

  library(png)
  library(ragg)

  make_raster_overlay <- function(plot, bbox,
                                  width_px = 3000,
                                  height_px = 3000,
                                  dpi = 300,
                                  bg = "transparent") {
    tf <- tempfile(fileext = ".png")

    ggsave(
      filename = tf,
      plot = plot,
      width = width_px / dpi,
      height = height_px / dpi,
      dpi = dpi,
      bg = bg,
      device = ragg::agg_png
    )

    img <- png::readPNG(tf)
    unlink(tf)

    annotation_raster(
      raster = as.raster(img),
      xmin = bbox["xmin"],
      xmax = bbox["xmax"],
      ymin = bbox["ymin"],
      ymax = bbox["ymax"]
    )
  }

  make_atlas_raster_layer <- function(atlas_sf,
                                      shade_sf = NULL,
                                      shade_size = 0.0005,
                                      shade_alpha = 0.005,
                                      width_px = 3000,
                                      dpi = 300,
                                      fill_limits = NULL,
                                      fill_scale = scale_fill_gradient2(
                                        low = scales::muted("blue"),
                                        mid = "white",
                                        high = scales::muted("red"),
                                        limits = fill_limits
                                      )) {
    bbox <- st_bbox(atlas_sf)

    ratio <- as.numeric((bbox["xmax"] - bbox["xmin"]) / (bbox["ymax"] - bbox["ymin"]))
    height_px <- max(1, round(width_px / ratio))

    p_raster <- ggplot() +
      geom_sf(
        data = atlas_sf,
        aes(fill = statistic, geometry = geometry),
        linewidth = 0.25,
        show.legend = FALSE
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


    if (!is.null(shade_sf)) {
      p_raster <- p_raster +
        geom_sf(
          data = shade_sf,
          size = shade_size,
          alpha = shade_alpha,
          show.legend = FALSE
        )
    }

    make_raster_overlay(
      plot = p_raster,
      bbox = bbox,
      width_px = width_px,
      height_px = height_px,
      dpi = dpi
    )
  }

  plot_brain_ests <- function(ests, tag = "a") {
    parcel_line_size = 0.1
    #bsize = 16

    terms_of_interest <- unique(ests$term)[!(unique(ests$term) %in% c(covariates, "(Intercept)"))]
    plots_of_terms <- list()
    tag_labs <- paste0(tag_prefix, tag_sep, tag, tag_sep, tolower(as.roman(1:length(terms_of_interest))))
    i = 1
    for (term_of_i in terms_of_interest) {

      if (group_n_title) {
        if(term_of_i %in% group_sizes$term){
          plot_sub <- paste0("(n = ", group_sizes |> filter(term == term_of_i) |> pull(n), ")")
        } else {
          plot_sub = ""
        }
      }


      ####### From here is new

      atlas_plot_df <- ests %>%
        filter(term == term_of_i) %>%
        right_join(atlas_geometry$atlas, by = "region")%>%
        sf::st_as_sf()

      fill_limits <- range(atlas_plot_df$statistic, na.rm = TRUE)

      bbox <- sf::st_bbox(atlas_geometry$atlas)

      raster_layer <- NULL
      if (rasterize) {
        raster_layer <- make_atlas_raster_layer(
          atlas_sf = atlas_plot_df,
          shade_sf = if (add_shade) atlas_geometry$shade else NULL,
          shade_size = shade_size,
          shade_alpha = shade_alpha,
          width_px = 1500,   # increase if needed
          dpi = 300,
          fill_limits = fill_limits
        )
      }

      p <- ggplot() +
        (
          if (rasterize) {
            raster_layer
          } else {
            geom_sf(
              data = atlas_plot_df,
              aes(fill = statistic, geometry = geometry),
              linewidth = l_width,
              show.legend = FALSE
            )
          }
        ) +
        (
          if (!rasterize && add_shade) {
            geom_sf(
              data = atlas_geometry$shade,
              size = shade_size,
              alpha = shade_alpha
            )
          } else {
            NULL
          }
        ) +
        # (
        #   if (rasterize) {
        #     geom_sf(
        #       data = atlas_plot_df,
        #       aes(geometry = geometry),
        #       fill = NA,
        #       colour = "black",
        #       linewidth = l_width,
        #       show.legend = FALSE
        #     )
        #   } else {
        #     NULL
        #   }
        # ) +
        coord_sf(
          xlim = c(bbox["xmin"], bbox["xmax"]),
          ylim = c(bbox["ymin"], bbox["ymax"]),
          expand = FALSE
        ) +
        theme_void(base_size = base_size_) +
        labs(
          fill = "t",
          title = stringr::str_to_title(stringr::str_replace(term_of_i, "_", " ")),
          subtitle = if (group_n_title) plot_sub else waiver()
        ) +
        theme(
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(color = "black", hjust = 0.5, vjust = brain_title_vjust, size = brain_title_size),
          plot.subtitle = element_text(hjust = 0.5, size = rel(0.7), vjust = brain_subtitle_vjust)
        ) +
        scale_fill_gradient2(
          low = scales::muted("blue"),
          mid = "white",
          high = scales::muted("red"),
          limits = fill_limits
        )

      # p <-  ests %>%
      #   filter(term == term_of_i) %>%
      #   right_join(atlas_geometry$atlas, by = "region") %>%
      #   ggplot() +
      #   (if (rasterize)
      #     rasterise(
      #       geom_sf(aes(
      #         fill = statistic,
      #         geometry = geometry), linewidth= l_width,
      #         show.legend = FALSE),
      #       dev = ggrastr_dev,
      #       scale = raster_scale_atlas
      #     )
      #     else
      #     geom_sf(aes(
      #       fill = statistic,
      #       geometry = geometry), linewidth= l_width,
      #       show.legend = FALSE)
      #   ) +
      #   (if (add_shade)
      #     (if (rasterize)
      #       shade_layer
      #       # rasterise(
      #       #   geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha),
      #       #   dev = ggrastr_dev,
      #       #   scale = raster_scale_shade
      #       # )
      #      else
      #        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
      #     )
      #    else NULL) +
      #   theme_void(base_size = base_size_)+
      #   labs(fill = 't', title = str_to_title(str_replace(term_of_i, "_", " ")),
      #        subtitle = if(group_n_title) plot_sub  else waiver()
      #        #tag = tag_labs[i]
      #        ) +
      #   theme(#legend.position = "",
      #     panel.background = element_rect(fill = "transparent", colour = NA),
      #     plot.background = element_rect(fill = "transparent", colour = NA),
      #     legend.background = element_rect(fill = "transparent", colour = NA),
      #     legend.box.background = element_rect(fill = "transparent", colour = NA),
      #     plot.title = element_text(color = "black", hjust = 0.5, vjust = brain_title_vjust),
      #     plot.subtitle = element_text(hjust = 0.5, size = rel(0.7), vjust = brain_subtitle_vjust),
      #   ) +
      #   scale_fill_gradient2(
      #     low = muted("blue"),
      #     mid = "white",
      #     high = muted("red")
      #   )

      if(!is.null(brain_names)) {
        if(!is.na(brain_names[i])) {
          b_name <- parse_plot_label(
            brain_names[i],
            parse = brain_names_parse[i],
            label_name = paste0("brain_names[", i, "]")
          )
          p <- p + ggtitle(b_name)
        }
      }

      i = i + 1

      if (show_networks) {
        p <- p +
          geom_sf(data = network_geometry %>% drop_na(),
                  aes(
                    #fill = region,
                    color = name,
                    geometry = geometry
                  ), alpha = 0, linewidth = 0.5,
                  show.legend = FALSE) +
          scale_color_manual(
            values = net_names %>% select(name, col) %>% deframe()
          )
      }
      plots_of_terms[[term_of_i]]<- p
    }

    return(plots_of_terms)
  }

  ## Gradient plots
  gradient_plots <- list()
  tag_labs <- paste0(tag_prefix, tag_sep, "a",tag_sep, tolower(as.roman(1:length(gradients))))
  i = 1
  for (stud in unique(gradient_data$study)) {
    for (grad in grad_char) {

      if (include_gradient_plots) {


        atlas_plot_df <- atlas_geometry$atlas %>%
          left_join(
            gradient_data %>% filter(study == stud),
            by = "region"
          ) %>%
          sf::st_as_sf()

        bbox <- sf::st_bbox(atlas_geometry$atlas)

        fill_limits <- range(atlas_plot_df[[grad]], na.rm = TRUE)

        raster_layer <- NULL
        if (rasterize) {
          raster_layer <- make_atlas_raster_layer(
            atlas_sf = atlas_plot_df %>%
              dplyr::rename(statistic = !!rlang::sym(grad)),
            shade_sf = if (add_shade) atlas_geometry$shade else NULL,
            shade_size = shade_size,
            shade_alpha = shade_alpha,
            width_px = 1500,
            dpi = 300,
            fill_limits = fill_limits,
            fill_scale = scale_fill_gradient2(
              low = gradient_colors[[grad]][1],
              mid = "white",
              high = gradient_colors[[grad]][2],
              limits = fill_limits
            )
          )
        }

        gradient_plots[[paste0(stud, "_", grad)]] <-
          ggplot() +
          (
            if (rasterize) {
              raster_layer
            } else {
              geom_sf(
                data = atlas_plot_df,
                aes(fill = .data[[grad]], geometry = geometry),
                linewidth = l_width,
                show.legend = FALSE
              )
            }
          ) +
          (
            if (!rasterize && add_shade) {
              geom_sf(
                data = atlas_geometry$shade,
                size = shade_size,
                alpha = shade_alpha
              )
            } else {
              NULL
            }
          ) +
          # (
          #   if (rasterize) {
          #     geom_sf(
          #       data = atlas_plot_df,
          #       aes(geometry = geometry),
          #       fill = NA,
          #       colour = "black",
          #       linewidth = l_width,
          #       show.legend = FALSE
          #     )
          #   } else {
          #     NULL
          #   }
          # ) +
          coord_sf(
            xlim = c(bbox["xmin"], bbox["xmax"]),
            ylim = c(bbox["ymin"], bbox["ymax"]),
            expand = FALSE
          ) +
          theme_void(base_size = base_size_) +
          labs(
            fill = "",
            title = stringr::str_to_title(stringr::str_replace(paste0("_", grad), "_", " "))
          ) +
          (
            if (grad_name) {
              labs(title = dplyr::case_when(
                grad == "gradient1" ~ "Sens-Assoc",
                grad == "gradient2" ~ "Vis-Mot",
                grad == "gradient3" ~ "Rep-Exec",
                TRUE ~ grad
              ))
            } else {
              NULL
            }
          ) +
          theme(
            legend.position = "",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(color = "black", hjust = 0.5, vjust = -3, size = brain_title_size)
          )



        # gradient_plots[[paste0(stud,"_", grad)]] <- gradient_data %>% filter(study==stud) %>%
        #   right_join(atlas_geometry$atlas, by = "region") %>%
        #   ggplot() +
        #   (if (rasterize)
        #     rasterise(
        #       geom_sf(aes(
        #         fill = .data[[grad]],
        #         geometry = geometry), linewidth= l_width,
        #         show.legend = FALSE),
        #       dev = ggrastr_dev,
        #       scale = raster_scale_atlas
        #     )
        #    else
        #      geom_sf(aes(
        #        fill = .data[[grad]],
        #        geometry = geometry), linewidth= l_width,
        #        show.legend = FALSE)
        #   ) +
        #   (if (add_shade)
        #     (if (rasterize)
        #       shade_layer
        #      # rasterise(
        #      #   geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha),
        #      #   dev = ggrastr_dev,
        #      #   scale = raster_scale_shade
        #      # )
        #      else
        #        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
        #     )
        #    else NULL) +
        #   theme_void(base_size = base_size_)+
        #   labs(fill = "", title = str_to_title(str_replace(paste0("_", grad), "_", " ")),
        #        #tag = tag_labs[i]
        #   ) +
        #   (if (grad_name)
        #     labs(title = case_when(
        #       grad == "gradient1" ~ "Sens-Assoc",
        #       grad == "gradient2" ~ "Vis-Mot",
        #       grad == "gradient3" ~ "Rep-Exec",
        #       TRUE ~ grad
        #     ))  else NULL) +
        #   theme(legend.position = "",
        #         panel.background = element_rect(fill = "transparent", colour = NA),
        #         plot.background = element_rect(fill = "transparent", colour = NA),
        #         legend.background = element_rect(fill = "transparent", colour = NA),
        #         legend.box.background = element_rect(fill = "transparent", colour = NA),
        #         plot.title = element_text(color = "black", hjust = 0.5, vjust = -3)
        #   ) +
        #   scale_fill_gradient2(
        #     low = gradient_colors[[grad]][1],
        #     mid = "white",
        #     high = gradient_colors[[grad]][2]
        #   )

        i = i +1
      } else {
        gradient_plots[[paste0(stud,"_", grad)]] <- gradient_data %>% filter(study==stud) %>%
          inner_join(atlas_geometry$atlas, by = "region") %>%
          ggplot() +
          (
            if (rasterize) {
              rasterise(
                geom_sf(aes(
                  fill = .data[[grad]],
                  geometry = geometry), alpha = 0, linewidth= l_width, color = NA,
                  show.legend = FALSE)
              )
            } else {
              geom_sf(aes(
                fill = .data[[grad]],
                geometry = geometry), alpha = 0, linewidth= l_width, color = NA,
                show.legend = FALSE)
            }
          ) +
          theme_void(base_size = base_size_)+
          labs(fill = "", title = str_to_title(str_replace(paste0("_", grad), "_", " ")),
               #tag = tag_labs[i]
          ) +
          theme(legend.position = "",
                panel.background = element_rect(fill = "transparent", colour = NA),
                plot.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                plot.title = element_text(color = NA, hjust = 0.5, size = brain_title_size)
          ) +
          scale_fill_gradient2(
            low = gradient_colors[[grad]][1],
            mid = "white",
            high = gradient_colors[[grad]][2]
          )
        i = i +1
      }


    }
  }



  vars_of_interest <- analysis_name
  n_analysis <- length(analysis_name)
  terms_of_interest <- unique(ests$term)[!(unique(ests$term) %in% c(covariates, "(Intercept)"))]
  n_terms <- length(terms_of_interest)
  n_plts_row <- length(terms_of_interest)*length(vars_of_interest)
  tag_labs <- c(length(analysis_name))

  n_gradients <- length(grad_char)

  labels <- letters[2:(n_analysis+1)]
  tag_labs  <- c()

  for (g in 0:(n_gradients - 1)) {
    num_seq <- (g * n_terms + 1):((g + 1) * n_terms)
    for (label in labels) {
      temp <- paste0(tag_prefix, tag_sep, label, tag_sep, num_seq)
      tag_labs  <- c(tag_labs , temp)
    }
  }

  plots <- list()
  count <- 1
  for (g in grad_char) {
    for (t in vars_of_interest) {
      for (term_ in terms_of_interest){
        plot_data  <- gradient_data %>% inner_join(ests %>% filter(term == term_), by = "region")

        i <- t
        j <- g

        # Swap axes if layout is "vertical"
        if (layout_construction == "vertical") {
          i <- g
          j <- t
        }

        lab_grad <- label_gradient_score
        lab_term <- str_to_sentence(paste(str_replace_all(t, "_", " "), "t-value"))

        if (layout_construction == "horizontal") {
          x_lab <- lab_term
          y_lab <- lab_grad
        } else {
          y_lab <- lab_term
          x_lab <- lab_grad
        }



        x_min <- min(plot_data[i])
        x_max <- max(plot_data[i])

        y_min <- min(plot_data[j])
        y_max <- max(plot_data[j])



        p <- plot_data %>%
          ggplot(aes(x = .data[[i]], y =.data[[j]],
                     color = if (!gray_out) name else "gray")) +
          (if (rasterize_scatters)
            ggrastr::rasterise(
              geom_point(alpha = ifelse(gray_out, 0.05, point_alpha),
                         size = point_size) ,
              dev = ggrastr_dev
            )
           else
             geom_point(alpha = ifelse(gray_out, 0.05, point_alpha),
                        size = point_size)
          )+
          stat_poly_line(se = FALSE,
                         linewidth = 0.5,
                         color = ifelse(gray_out, "#808080", "#323232")) +
          labs(#title = term_,
            x = x_lab,
            y = y_lab,
            #tag = tag_labs[count],
            color = "Network") +
          xlim(x_min, x_max)+
          #ylim(y_min, y_max)+
          theme_bw(base_size = base_size_) +
          theme(
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            axis.text = element_text(size = axis_text_size),
            axis.title = element_text(size = axis_title_size)
          ) +
          scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
          scale_y_continuous(limits = c(y_min, y_max), position = ifelse(right_term_side, "right", "left"))

        if (spintest) {
          r <- cor(plot_data[[i]], plot_data[[j]], method = "pearson")
          p_val <- perm_sphere_p(plot_data[[i]], plot_data[[j]], perm.id = perms, corr.type='pearson')


          # p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("<i>p<sub>spin</sub></i> &lt; ", "<i>p<sub>spin</sub></i> = ", "<i>p<sub>spin</sub></i> &gt; "))(p_val)
          # label_corr <- paste0("<i>r</i> = ", r |> round(2), " ,&#8203;&nbsp;" , p_lab)
          #
          #
          # p <- p + labs(subtitle = label_corr) +
          #   theme(plot.subtitle = element_markdown(size = rel(r_spin_size), hjust = 0.95, margin = margin(b = 3),
          #                                          color = ifelse(gray_out, scales::alpha("black", 0.2), "black")))


          p_lab <- scales::label_pvalue(
            accuracy = 0.001,
            prefix = c(
              "<i>p<sub>spin</sub></i> &lt; ",
              "<i>p<sub>spin</sub></i> = ",
              "<i>p<sub>spin</sub></i> &gt; "
            )
          )(p_val)

          label_corr <- paste0("<i>r</i> = ", round(r, 2), ", ", p_lab)

          p <- p + labs(subtitle = label_corr) +
            theme(
              plot.subtitle = element_markdown(
                size = rel(r_spin_size),
                hjust = 0.95,
                margin = margin(b = scatter_title_vjust),
                color = ifelse(gray_out, scales::alpha("black", 0.2), "black")
              )
            )

          #   subtitle = label = label_corr, alpha = ifelse(gray_out, 0.2, 1)
          # p <- p + ggpp::annotate(geom = "text_npc", label = label_corr,
          #                         alpha = ifelse(gray_out, 0.2, 1),
          #                           size = r2_size,
          #                           npcx = "left", npcy = "top",
          #                           family = "sans",
          #                           parse = TRUE)
        } else {
          p <- p + stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                                  str_remove(after_stat(rr.confint.label), "95% CI "),
                                                  sep = "*\" \"*")),
                                parse = TRUE, color = "#323232", label.x = "left", label.y = "top", size = r2_size)
        }

        if (gray_out) {
          p <- p + ggpp::annotate("text_npc",
                                  npcx = "middle",
                                  npcy = "middle",
                                  label = "?",
                                  size = 28
          )
        }


        if (side_color_bar) {
          require(ggside)
          x_range_df <- data.frame(x = seq(x_min, x_max, length.out = 100))
          colnames(x_range_df)[1] <- i

          y_range_df <- data.frame(y = seq(y_min, y_max, length.out = 100))
          colnames(y_range_df)[1] <- j

          if(layout_construction == "horizontal") {
            alpha_switch_xside <- ifelse(g == tail(grad_char, 1), 1, 0)
            alpha_switch_yside <- ifelse(term_ == terms_of_interest[1], 1, 0)
            x_axis_colorbar <- c(muted("blue"), muted("red"))
            y_axis_colorbar <- gradient_colors[[g]]
          }
          if(layout_construction == "vertical") {
            alpha_switch_xside <- ifelse(term_ == tail(terms_of_interest, 1), 1, 0)
            alpha_switch_yside <- ifelse(g == grad_char[1], 1, 0)
            x_axis_colorbar <- gradient_colors[[g]]
            y_axis_colorbar <- c(muted("blue"), muted("red"))
          }

          if(right_term_side & layout_construction == "vertical") alpha_switch_yside <- ifelse(g == tail(grad_char, 1), 1, 0)



          p <- p +
            ggrastr::rasterise(
              geom_xsidetile(data = x_range_df, aes(x = .data[[i]], y = 0, fill = .data[[i]]),
                             alpha = alpha_switch_xside,
                             show.legend = FALSE,
                             inherit.aes = FALSE),
              dev = "ragg",
              dpi = 300
            ) +
            scale_fill_gradient2(
              low = x_axis_colorbar[1],
              mid = "white",
              high = x_axis_colorbar[2] 
            ) +
            ggnewscale::new_scale_fill() +
            ggrastr::rasterise(
              geom_ysidetile(data = y_range_df, aes(y = .data[[j]], x = 0, fill = .data[[j]]),
                             alpha = alpha_switch_yside,
                             show.legend = FALSE,
                             inherit.aes = FALSE),
              dev = "ragg",
              dpi = 300
            ) +
            scale_fill_gradient2(
              low = y_axis_colorbar[1],
              mid = "white",
              high = y_axis_colorbar[2]
            ) +
            theme_ggside_void() +
            ggside(x.pos = "bottom", y.pos = ifelse(right_term_side, "right", "left")) +
            theme(ggside.panel.scale = 0.04)
        }

        plots[[count]] <- p
        count <- count + 1
      }
    }
  }

  n_terms = length(unique(ests$term)[!(unique(ests$term) %in% c(covariates, "(Intercept)"))])
  n_analysis <- length(list_of_parcel_data)

  list_of_brain_plots_ests <- list()
  letter_tag <- letters[2:(n_analysis+1)]
  i = 1
  for (analysis in names(list_of_ests)){
    list_of_brain_plots_ests[[analysis]] <- plot_brain_ests(list_of_ests[[analysis]], tag = letter_tag[i])
    i = i + 1
  }


  if (layout_construction == "horizontal") {
    n_plot_cols <- length(1:((n_terms*n_analysis) + (n_analysis - 1) + 1))
    col_widths <- c(1, rep(c(rep(1, n_terms), plot_spacing), n_analysis))[-(n_plot_cols+1)]
    empty_cols <- seq(1, n_plot_cols, by = n_terms+1)[-1]

    n_plot_rows <- n_gradients + 1 + (n_gradients)
    row_heights <- rep(1, n_plot_rows)
    empty_rows <- seq(2, n_plot_rows, by = 2)

    row_heights[empty_rows] <- empty_row_height

    layout <- c(
      area(3, 1)  # First gradient
    )

    if (n_gradients>1) {
      for (g in seq(5, n_plot_rows, by =2)){
        layout <- c(layout, area(g, 1))
      }
    }

    for (col in 2:n_plot_cols){
      if (col %in% empty_cols) {
      } else {
        layout <- c(layout, area(1, col))
      }
    }

    for(i in seq(3, (n_gradients*2-1)+2, by = 2)){
      for (j in 2:((n_terms*n_analysis) + (n_analysis - 1) + 1)) {
        if ((j %in% empty_cols)) {
        } else {
          layout <- c(layout, area(i, j))
        }
      }
    }

    for(empt in empty_cols){
      layout <- c(layout, area(1, empt , b = n_gradients + 1))
    }

    if (padding != 0) {
      for(pad in 1:padding){
        layout <- c(layout, area(1, n_plot_cols + pad, b = n_gradients + 1))
      }
    }


  }


  if (layout_construction == "vertical") {
    n_plot_rows <- (n_terms * n_analysis) + (n_analysis - 1) + 1
    row_heights <- c(1, rep(c(rep(1, n_terms), plot_spacing), n_analysis))[-(n_plot_rows+1)]
    col_widths <- NULL
    empty_rows <- seq(1, n_plot_rows, by = n_terms + 1)[-1]

    if (right_term_side) {
      layout <- c(
        area(1, 1)
      )

      for (g in 2:(length(gradients))) {
        layout <- c(layout, area(1, g))
      }

      for (row in 2:n_plot_rows) {
        if (!(row %in% empty_rows)) {
          layout <- c(layout, area(row, (length(gradients) + 1)))
        }
      }

      for (j in 1:(length(gradients))) {
        for (i in 2:n_plot_rows) {
          if (!(i %in% empty_rows)) {
            layout <- c(layout, area(i, j))
          }
        }
      }

    } else {
      layout <- c(
        area(1, 2)
      )

      for (g in 3:(length(gradients) + 1)) {
        layout <- c(layout, area(1, g))
      }

      for (row in 2:n_plot_rows) {
        if (!(row %in% empty_rows)) {
          layout <- c(layout, area(row, 1))
        }
      }

      for (j in 2:(length(gradients) + 1)) {
        for (i in 2:n_plot_rows) {
          if (!(i %in% empty_rows)) {
            layout <- c(layout, area(i, j))
          }
        }
      }
    }

  }


  filt_char = as_label(filter_criteria )


  brain_plots <- unlist(list_of_brain_plots_ests, recursive = FALSE)

  plots_to_include <- c(gradient_plots, brain_plots, plots)

  if (plt_subtitle) {
    #f <- mod_formula
    rhs <- ests$model_formula[1]
    rhs <- str_remove(rhs, "~")
    #rhs <- paste0(rhs)[2]
    rhs <- gsub("\\*", "×", rhs)
    rhs_expr <- parse(text = rhs)[[1]]

    if (!is.null(subtit_lookup)){
      for (var_name in names(subtit_lookup)) {
        rhs_expr <- str_replace_all(
          rhs_expr,
          regex(paste0("(?<![A-Za-z0-9_])", str_replace_all(var_name, "([\\W])", "\\\\\\1"), "(?![A-Za-z0-9_])")),
          subtit_lookup[[var_name]]
        )
      }
    }

    #subtit_expr <- parse(text = expr_str)[[1]]
    subtit_expr <- bquote(italic(FCS[parcel] ~ '~' ~ .(rhs_expr)))
  } else {
    subtit_expr <- waiver()
  }

  if (!is.null(plt_title)) {
    n = ests %>% pull(n) %>% unique()
    if (isTRUE(plt_title_parse)) {
      plt_title <- paste0(plt_title, "~'(N = ", n, ")'")
    } else {
      plt_title <- paste0(plt_title, " (N = ", n, ")")
    }
    plt_title <- parse_plot_label(plt_title, parse = plt_title_parse, label_name = "plt_title")
  } else {
    plt_title <- waiver()
  }

  p <- Reduce(`+`, plots_to_include) +
    plot_annotation(title = plt_title, subtitle = subtit_expr,
                    #tag_levels = c('A', '1'),
                    theme = theme(
                      plot.title = element_text(size = plt_title_size,
                        hjust = plt_title_hjust,
                        margin = margin(l = title_lmargin)),
                      plot.subtitle = element_text(size = plt_subtitle_size,
                                                   hjust = plt_subtitle_hjust,
                                                   margin = margin(l = title_lmargin),
                                                   vjust = -0.05,
                                                   family = "mono",
                                                   face = "italic")
                    )) +
    plot_layout(design = layout, axis_titles = "collect", axes = "collect", guides = "collect",
                widths = col_widths,
                heights = row_heights
    ) &
    theme(legend.position = "",
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.tag.position  = if (right_term_side) c(0.1, 0.96) else c(0.9, 0.96),
          #plot.tag = element_text(size = 8, hjust = 0, vjust = 0)
    )

  if (rectangle) {
    p <- p + plot_annotation(theme = theme(plot.background = element_rect(color = "black", fill = NA, linewidth = 0.5)))
  }

  if (plot_net_legend) {
    net_legend <- get_net_legend(b_size = base_size_)
    p <- ggdraw() +
      draw_plot(p) +
      draw_plot(net_legend, x = net_legend_x, y = net_legend_y, width = 0.1, height = 0.05)
  }

  list(plot = p, n = ests %>% pull(n) %>% unique(), model_formula = ests %>% pull(model_formula) %>% unique(), tmaps = ests)

}
