

rois_vec <- if (is.data.frame(rois)) {
if ("region" %in% names(rois)) {
  as.character(rois[["region"]])
} else {
  as.character(rois[[1]])
}
} else {
as.character(rois)
}

diffusion_map_embedding <- function(L, t, N, alpha = 0.5) {
  
  
  D <- diag(rowSums(L))
  
  D_alpha_inv <- diag(diag(D)^(-alpha))
  
  L_star <- D_alpha_inv %*% L %*% D_alpha_inv
  
  D_star <- diag(rowSums(L_star))
  
  
  D_star_inv <- diag(1 / diag(D_star))
  P <- D_star_inv %*% L_star
  
  
  eig <- eigen(P)
  lambda <- eig$values
  psi <- eig$vectors
  
  total_variance <- sum(lambda)
  variance_explained <- lambda / total_variance
  
  for (i in 1:N) {
    if (t == 0) {
      lambda[i] <- lambda[i] / (1 - lambda[i])
    } else {
      lambda[i] <- lambda[i]^t
    }
  }
  
  psi <- sweep(psi, 2, psi[, 1], "/")  
  phi <- sweep(psi, 2, lambda, "*")
  
  
  return(list(dmaps = phi[, 2:(N + 1)], var_exp = variance_explained))
}


reorder_gradients <- function(original, derived) {
  
  num_components <- ncol(original)
  correlation_matrix <- matrix(NA, nrow = num_components, ncol = num_components)
  
  for (i in 1:num_components) {
    for (j in 1:num_components) {
      correlation_matrix[i, j] <- cor(derived[, i], original[, j])
    }
  }
  
  
  max_cor_indices <- apply(correlation_matrix, 2, function(x) which.max(abs(x))) 
  
  # Reorder the derived components based on max correlations
  reordered_derived <- derived[, max_cor_indices]
  
  # Return the reordered and sign-aligned derived components
  return(reordered_derived)
}

align_gradients <- function(original, derived){
  
  num_components <- ncol(derived)
  for (i in 1:num_components) {
    # If the correlation is negative, flip the sign of the derived component
    if (cor(original[, i], derived[, i]) < 0) {
      derived[, i] <- -derived[, i]
    }
  }
  return(derived)
}

GRADIENT_LIMITS <- c(-2, 2)
GRADIENT_BREAKS <- seq(-2, 2, by = 1)

gradient_color_scale <- function() {
  scale_fill_viridis_c(
    option = "viridis",
    direction = 1,
    limits = GRADIENT_LIMITS,
    breaks = GRADIENT_BREAKS,
    oob = scales::squish,
    na.value = "grey85",
    name = "Gradient value",
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(35, "mm"),
      barwidth = grid::unit(5, "mm"),
      ticks = TRUE
    )
  )
}

#atlas_geometry_og <- readRDS("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/schaef_ggseg2.rds")
get_gradients <- function(connectome_ests, 
                          n_gradients = 1:3,
                          on_affinity = TRUE, 
                          threshold = 0.5,
                          atlas_geometry = readRDS("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/schaef400_ggseg2.rds"),
                          method = c("diffusion", "pca"),
                          similarity_method = "cosine",
                          reorder_components = FALSE,
                          align_components = TRUE, 
                          visualize = TRUE,
                          side_density = TRUE,
                          reference_gradients = marg_gradients) {
  
  require(scales)
  
  zero_out_mat <- function(mat, thresh) {
    library(matrixStats)
    thresholds <- colQuantiles(abs(mat), probs = thresh)
    mat[abs(mat) <= thresholds[col(mat)]] <- 0
    return(mat)
  }
  
  method = match.arg(method)
  grad_char <- paste0("gradient", n_gradients)
  
  reference_gradients <- reference_gradients[, n_gradients]
  
  grads_study <- list()
  var_exp_df <- matrix(nrow = length(connectome_ests), ncol = length(n_gradients))
  rownames(var_exp_df) <- names(connectome_ests)
  colnames(var_exp_df) <- grad_char
  for (ave_conn_name in names(connectome_ests)) {
    
    ave_conn <- connectome_ests[[ave_conn_name]]
    
    if (on_affinity) {
      library(proxy)
      
      L <- zero_out_mat(ave_conn, threshold)
      L[L<0] <- 0
      L <- proxy::simil(L, method = similarity_method, by_rows = FALSE)
      L <- as.matrix(L)
      diag(L) <- 1
    } else {
      L <- zero_out_mat(ave_conn, threshold)
      L[L<0] <- 0
    }
    
    
    if (method == "diffusion") {
      dmap <- diffusion_map_embedding(L, t = 0, N = 10)
      var_exp <- dmap$var_exp
      grads <- dmap$dmaps[, n_gradients]
    }
    
    if (method == "pca") {
      
      pca_result <- prcomp(L, scale. = TRUE)
      eigenvalues_pca <- (pca_result$sdev)^2
      total_variance_pca <- sum(eigenvalues_pca)
      var_exp <- eigenvalues_pca / total_variance_pca
      
      grads <- pca_result$x[, n_gradients]
      
    }
    
    if (reorder_components) {
      
      grads <- reorder_gradients(reference_gradients, grads)
      
    }
    
    if (align_components) {
      
      grads <- align_gradients(reference_gradients, grads)
      
    }
    colnames(grads) <- colnames(reference_gradients)
    grads <- grads %>% as_tibble()
    grads_study[[ave_conn_name]] <- grads
    
    var_exp_df[ave_conn_name, ] <- var_exp[n_gradients]
    
  }
  
  if (visualize) {
    library(patchwork)
    library(ggside)
    require(tidyverse)
    require(scales)
    require(patchwork)
    require(ggpmisc)
    require(sf)
    
    net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                            col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), #"#DCF8A4"
                            label = c(1:7))

    
    
    std_grad_plots <- list()
    
    for (grad in grad_char) {
      std_grad_plots[[grad]] <- reference_gradients %>%
  mutate(region = rois_vec) %>%
  inner_join(atlas_geometry$atlas, by = "region") %>%
  ggplot() +
  geom_sf(
    aes(
      fill = .data[[grad]],
      geometry = geometry
    ),
    linewidth = 0.2,
    show.legend = TRUE
  ) +
  geom_sf(
    data = atlas_geometry$shade,
    linewidth = 0.01,
    alpha = 0.01
  ) +
  gradient_color_scale() +
  labs(
    title = str_to_title(
      str_replace(
        paste0("Margulies_", grad),
        "_",
        " "
      )
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    plot.title = element_text(
      color = "black",
      hjust = 0.5
    )
  )
        # scale_fill_gradient2(
        #   low = muted("blue"),
        #   mid = "white",
        #   high = muted("red") 
        # ) 
      
    }
    
    gradient_plots <- list()
    
    for (ave_conn_name in names(grads_study)) {
      print(paste0("Visualizing gradients for ", ave_conn_name))
      for (grad in grad_char) {
        print(paste0("Visualizing ", grad))
        gradient_plots[[paste0(ave_conn_name,"_", grad)]] <- grads_study[[ave_conn_name]] %>% 
          mutate(region = rois_vec) %>% 
          inner_join(atlas_geometry$atlas, by = "region") %>%
          ggplot() +
          geom_sf(aes(
            fill = .data[[grad]],
            geometry = geometry), linewidth= 0.2,
            show.legend = TRUE)+
          geom_sf(data = atlas_geometry$shade, size = 0.01, alpha = 0.01) +
          theme_void()+
          labs(fill = "", title = str_to_title(str_replace(paste0(ave_conn_name, "_", grad), "_", " ")),
               subtitle = paste0(round(var_exp_df[ave_conn_name, grad]*100), "% explained variance")
          ) +
          theme(legend.position = "",
                panel.background = element_rect(fill = "transparent", colour = NA),
                plot.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                plot.title = element_text(color = "black", hjust = 0.5),
                plot.subtitle = element_text(color = "black", hjust = 0.5)
          ) +
          scale_fill_viridis_c(
          option = "viridis",
          direction = 1,
          begin = 0,
          end = 1,
          name = "Gradient value",
          na.value = "grey85"
        )
          # scale_fill_gradient2(
          #   low = muted("blue"),
          #   mid = "white",
          #   high = muted("red") 
          # ) 
        
      }
    }
    
    plots <- list()
    #count <- 1
    for (std_grad in grad_char) {
      
      for (ave_conn_name in names(grads_study)) {
        
        plot_data <- reference_gradients %>% mutate(region = rois_vec, study = "margulies") %>% 
          pivot_longer(starts_with("gradient"), names_to = "gradient", values_to = "value") %>% 
          filter(gradient == std_grad) %>% 
          pivot_wider(names_from = study, values_from = "value") %>% 
          inner_join(grads_study[[ave_conn_name]] %>% mutate(region = rois_vec, study = ave_conn_name), by = "region") %>% 
          inner_join(data.frame(region = rois_vec, label = yeo_msk) %>% inner_join(net_names, by = "label") %>% select(region, name), by = "region") 
        
        for (grad in grad_char) {
          #plot_data  <- std_grad_long %>% filter(gradient == grad) %>% pivot_wider(names_from = "study", values_from = "value")
          
          p <- plot_data %>% 
            ggplot(aes(x = .data[[grad]], y = margulies,
                       color = name)) +
            geom_point(alpha = 0.2) +
            stat_poly_eq(color = "#323232", label.x = "left", label.y = "top", size = 5) +
            stat_poly_line(se = FALSE, color = "#323232") +
            labs(
              #title = str_to_title(grad),
              y = "Margulies",
              x = str_to_title(ave_conn_name),
              #tag = tag_labs[count],
              color = "Network") +
            theme_bw() +
            theme(
              legend.position = "",
              panel.background = element_rect(fill = "transparent", colour = NA),
              plot.background = element_rect(fill = "transparent", colour = NA),
              legend.background = element_rect(fill = "transparent", color = NA),
              legend.box.background = element_rect(fill = "transparent", colour = NA)) +
            scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
            scale_y_continuous(limits = c(min(plot_data["margulies"]), max(plot_data["margulies"])))+
            scale_x_continuous(limits = c(min(plot_data[grad]), max(plot_data[grad])))
          
          
          if (side_density){
            if(std_grad ==grad_char[1]) {
              if(grad == tail(grad_char, 1)){
                p <- p +
                  geom_xsideboxplot(aes(y = as.numeric(factor(name))), orientation = "y", outlier.shape = NA ) +
                  geom_ysideboxplot(aes(x = as.numeric(factor(name))), orientation = "x", outlier.shape = NA ) +
                  #geom_xsidedensity(aes(y = after_stat(density)), show.legend = FALSE) +
                  #geom_ysidedensity(aes(x = after_stat(density)), show.legend = FALSE) +
                  theme_ggside_void()
                
              } else {
                p <- p +
                  geom_xsideboxplot(aes(y = as.numeric(factor(name))), orientation = "y", outlier.shape = NA ) +
                  geom_ysidedensity(color = NA, show.legend = TRUE) +
                  theme_ggside_void()
              }
            } else if(grad == tail(grad_char, 1) & std_grad != grad_char[1]){
              p <- p + 
                geom_ysideboxplot(aes(x = as.numeric(factor(name))), orientation = "x", outlier.shape = NA ) +
                geom_xsidedensity(color = NA, show.legend = TRUE) +
                theme_ggside_void()
            } else {
              p <- p +
                geom_xsidedensity(color = NA, show.legend = TRUE) +
                geom_ysidedensity(color = NA, show.legend = TRUE) +
                theme_ggside_void()
            }
          }
          plot_name <- paste(
  "scatter",
  ave_conn_name,
  grad,
  "vs_margulies",
  std_grad,
  sep = "_"
)

plots[[plot_name]] <- p
          
        }
      }
    }
    
    
    # Store all plot collections
plot_output <- list(
  margulies_brains = std_grad_plots,
  study_brains = gradient_plots,
  scatterplots = plots
)

# Print every plot separately
for (plot_group in names(plot_output)) {

  for (plot_name in names(plot_output[[plot_group]])) {

    print(
      plot_output[[plot_group]][[plot_name]]
    )
  }
}
    
  }
  
  gradient_data <- reference_gradients %>% mutate(study = "margulies", region = rois_vec)%>% 
    mutate(method = NA,
           affinity = NA,
           sim_method = NA,
           threshold = NA)
  print("Calculating gradients over various parameters")
  for(ave_conn_name in names(connectome_ests)){
    gradient_data <- rbind(gradient_data, grads_study[[ave_conn_name]] %>% mutate(study = ave_conn_name, 
                                                                                  region = rois_vec,
                                                                                  method = method,
                                                                                  affinity = on_affinity,
                                                                                  sim_method = ifelse(on_affinity, similarity_method, NA),
                                                                                  threshold = threshold))
  }
  
  gradient_data <- gradient_data %>% 
    inner_join(data.frame(region = rois_vec, label = yeo_msk) %>% 
                 inner_join(net_names, by = "label") %>% select(region, name), by = "region") 
  
  var_exp_df <-  var_exp_df %>% as_tibble(rownames = NA) %>% rownames_to_column("study") %>% 
    mutate(method = method,
           affinity = on_affinity,
           sim_method = ifelse(on_affinity, similarity_method, NA),
           threshold = threshold)
  
  return(list(gradients = gradient_data, varexp = var_exp_df))
}




# plot_grads_over_params <- function(connectome_list,
#                                    atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds"),
#                                    param_grid = NULL, 
#                                    n_gradients = 1:3) {
#   
#   require(sf)
#   
#   net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
#                           col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), #"#DCF8A4"
#                           label = c(1:7))
# 
#   
#   connectome_est_list <- connectome_list
#   
#   if (is.null(param_grid)){
#     params <- expand_grid(method = c("pca", "diffusion"), affinity = c(FALSE, TRUE), sim_method = "cosine", threshold = c(0.0, 0.25, 0.5, 0.75)) %>% 
#       filter(!(method == "diffusion" & !affinity)) %>% 
#       filter(!(method == "pca" & affinity)) %>% 
#       mutate(sim_method = ifelse(!affinity, NA, sim_method))
#   } else {
#     params <- param_grid
#   }
#   
#   gradient_data <- c()
#   varexp_df <- c()
#   for (i in 1:nrow(params)) {
#     param_i <- params[i, ]
#     grad_list <- get_gradients(connectome_ests = connectome_est_list,
#                                n_gradients = n_gradients,
#                                threshold = param_i$threshold,
#                                similarity_method = param_i$sim_method,
#                                on_affinity = param_i$affinity,
#                                method = param_i$method,
#                                visualize = FALSE,
#                                side_density = FALSE)
#     
#     gradient_data <- rbind(gradient_data, grad_list$gradients)
#     varexp_df <- rbind(varexp_df, grad_list$varexp)
#   }
#   
#   gradient_data <- gradient_data %>% distinct()
#   grad_char <- paste0("gradient", n_gradients)
#   
#   std_grad_plots <- list()
#   
#   for (grad in grad_char) {
#     std_grad_plots[[grad]] <- gradient_data %>% filter(study=="margulies") %>% 
#       mutate(region = rois_vec) %>% 
#       inner_join(atlas_geometry$atlas, by = "region") %>%
#       ggplot() +
#       geom_sf(aes(
#         fill = .data[[grad]],
#         geometry = geometry), linewidth= 0.2,
#         show.legend = FALSE)+
#       geom_sf(data = atlas_geometry$shade, size = 0.01, alpha = 0.01) +
#       theme_void()+
#       labs(fill = "", title = str_to_title(str_replace(paste0("Margulies_", grad), "_", " "))) +
#       theme(legend.position = "",
#             panel.background = element_rect(fill = "transparent", colour = NA),
#             plot.background = element_rect(fill = "transparent", colour = NA),
#             legend.background = element_rect(fill = "transparent", colour = NA),
#             legend.box.background = element_rect(fill = "transparent", colour = NA),
#             plot.title = element_text(color = "black", hjust = 0.5)
#       ) +
#       scale_fill_gradient2(
#         low = muted("blue"),
#         mid = "white",
#         high = muted("red") 
#       ) 
#     
#   }
#   
#   grad_text <- list()
#   
#   for(g_txt in paste0("G", n_gradients)) {
#     grad_text[[g_txt]] <- 
#       ggplot() +
#       theme_void() +  
#       annotate("text", x = 0.5, y = 0.5, label = g_txt, size = rel(20), 
#                fontface = "bold",
#                color = "#323232",
#                hjust = 0.5, vjust = 0.5)
#   }
#   
#   
#   
#   gradient_plots <- list()
#   for (stud in names(connectome_est_list)){
#     plt_idx <- 0
#     for (grad in grad_char) {
#       for (params_row in 1:nrow(params) ){
#         
#         pars <- params[params_row, ]
#         
#         plot_pars <- pars %>% mutate(affinity = ifelse(is.na(sim_method), "No", sim_method), 
#                                      method = ifelse(method=="pca", "PCA", "Diffusion")) %>% select(-sim_method) %>% 
#           rename(thresh = threshold)
#         par_char <- paste(names(plot_pars), plot_pars[1, ], sep = ": ") %>% str_to_sentence() 
#         
#         
#         variance_explained <- varexp_df %>% semi_join(pars, by = join_by(method, affinity, sim_method, threshold)) %>% 
#           filter(study == stud)
#         
#         
#         plt_idx <- plt_idx + 1
#         p <- gradient_data %>% 
#           filter(study == stud) %>% 
#           semi_join(pars, by = join_by(method, affinity, sim_method, threshold)) %>% 
#           inner_join(atlas_geometry$atlas, by = "region") %>%
#           ggplot() +
#           geom_sf(aes(
#             fill = .data[[grad]],
#             geometry = geometry), linewidth= 0.2,
#             show.legend = FALSE)+
#           geom_sf(data = atlas_geometry$shade, size = 0.01, alpha = 0.01) +
#           theme_void()+
#           labs(fill = "", title = paste0(plot_pars[, "method"],", ", par_char[3]),
#                subtitle = ifelse(grad == grad_char[1], par_char[2], NA),
#                caption = paste0(round(variance_explained %>% pull(grad)*100), "% explained variance")
#           ) +
#           theme(legend.position = "",
#                 panel.background = element_rect(fill = "transparent", colour = NA),
#                 plot.background = element_rect(fill = "transparent", colour = NA),
#                 legend.background = element_rect(fill = "transparent", colour = NA),
#                 legend.box.background = element_rect(fill = "transparent", colour = NA),
#                 plot.title = element_text(color = "black"),
#                 plot.subtitle = element_text(color = "black")
#           ) +
#           scale_fill_gradient2(
#             low = muted("blue"),
#             mid = "white",
#             high = muted("red") 
#           ) 
#         
#         if (grad != grad_char[1]) p <- p + theme(plot.title = element_blank(),
#                                                  plot.subtitle = element_blank())
#         gradient_plots[[stud]][[plt_idx]] <- p
#         
#       }
#     }
#   }
#   
#   
#   plots <- list()
#   
#   ref_grads <- gradient_data %>% filter(study == "margulies")
#   
#   for (stud in names(connectome_est_list)){
#     plt_idx <- 0
#     plot_data_stud <- gradient_data %>% filter(study %in% c(stud))
#     for (grad in grad_char) {
#       
#       plot_data_grad <- plot_data_stud %>% pivot_longer(starts_with("gradient"), names_to = "gradient", values_to = "value") %>% 
#         filter(gradient == grad) %>% 
#         pivot_wider(names_from = "study", values_from = "value")
#       
#       for (params_row in 1:nrow(params) ){
#         pars <- params[params_row, ]
#         
#         plot_data <- plot_data_grad %>% semi_join(pars) %>% 
#           mutate(margulies = ref_grads %>% pull(grad))
#         
#         p <- plot_data %>% 
#           ggplot(aes(x = .data[[stud]], y = margulies,
#                      color = name)) +
#           geom_point(alpha = 0.2) +
#           stat_poly_eq(color = "#323232", label.x = "left", label.y = "top", size = 5) +
#           stat_poly_line(se = FALSE, color = "#323232") +
#           labs(
#             #title = str_to_title(grad),
#             y = "Margulies",
#             x = str_to_title(stud),
#             #tag = tag_labs[count],
#             color = "Network") +
#           theme_bw() +
#           theme(
#             legend.position = "",
#             panel.background = element_rect(fill = "transparent", colour = NA),
#             plot.background = element_rect(fill = "transparent", colour = NA),
#             legend.background = element_rect(fill = "transparent", color = NA),
#             legend.box.background = element_rect(fill = "transparent", colour = NA)) +
#           scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
#           scale_y_continuous(limits = c(min(plot_data["margulies"]), max(plot_data["margulies"])))+
#           scale_x_continuous(limits = c(min(plot_data[stud]), max(plot_data[stud])))
#         
#         plt_idx <- plt_idx + 1
#         plots[[stud]][[plt_idx]] <- p
#       }
#     }
#   }
#   
#   
#   no_of_grads = length(n_gradients)
#   n_study <- length(connectome_est_list)
#   n_plot_cols <- nrow(params)#((no_of_grads*n_study) + (n_study - 1))
#   row_heights <- head(rep(c(1, 1, 0.1), no_of_grads), -1)
#   
#   param_combos <- length(unique(params$threshold))
#   n_total_cols <- n_plot_cols + floor(n_plot_cols / param_combos)
#   col_widths <- head(c(1, rep(c(rep(1, param_combos), 0.3), floor(n_plot_cols / param_combos))), -1)
#   
#   layout <- c(
#     area(2, 1))
#   
#   for (g in seq(5, no_of_grads + 5, by = 3)) {
#     layout <- c(layout, area(g, 1))
#   }
#   
#   for (g_txt in seq(1, no_of_grads + 5, by = 3)) {
#     layout <- c(layout, area(g_txt, 1))
#   }
#   
#   
#   rows <- layout$t - 1
#   
#   
#   for (row in rows) {
#     if (row %% 3 == 0) next
#     for (col in 2:(n_total_cols + 1)) {
#       if (((col - 2) %% (param_combos + 1)) == param_combos) next  
#       layout <- c(layout, area(row, col))
#     }
#   }
#   
#   for (row in rows + 1) {
#     if (row %% 3 == 0) next
#     for (col in 2:(n_total_cols + 1)) {
#       if (((col - 2) %% (param_combos + 1)) == param_combos) next  
#       layout <- c(layout, area(row, col))
#     }
#   }
#   
#   
#   
#   
#   
#   patch_plots <- list()
#   for (stud in names(connectome_est_list)){
#     pp <- 
#       std_grad_plots[[1]] + std_grad_plots[2:no_of_grads] +
#       grad_text +
#       gradient_plots[[stud]][1:(nrow(params)*no_of_grads)] +
#       plots[[stud]][1:(nrow(params)*no_of_grads)] +
#       plot_layout(design = layout, axis_titles = "collect", axes = "collect_y", guides = "collect",
#                   heights = row_heights,
#                   widths = col_widths
#       ) & 
#       theme(legend.position = "",
#             panel.background = element_rect(fill = "transparent", colour = NA),
#             plot.background = element_rect(fill = "transparent", colour = NA),
#             legend.background = element_rect(fill = "transparent", color = NA),
#             legend.box.background = element_rect(fill = "transparent", colour = NA),
#             plot.tag.position  = c(.9, .96))  
#     
#     patch_plots[[stud]] <- pp
#   }
#   return(patch_plots)
# }