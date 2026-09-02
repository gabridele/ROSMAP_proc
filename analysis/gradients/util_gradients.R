

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

GRADIENT_LIMITS <- c(-15, 20)
GRADIENT_BREAKS <- seq(-20, 20, by = 10)

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

get_gradients <- function(connectome_ests, 
                          n_gradients = 1:3,
                          on_affinity = TRUE, 
                          threshold = 0.5,
                          atlas_geometry = readRDS("analysis/gradients/atlas_data/schaef400_ggseg2.rds"),
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
    legend.position = "bottom",
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
          theme(legend.position = "bottom",
                panel.background = element_rect(fill = "transparent", colour = NA),
                plot.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                plot.title = element_text(color = "black", hjust = 0.5),
                plot.subtitle = element_text(color = "black", hjust = 0.5)
          ) +
          
        #   scale_fill_viridis_c(
        #   option = "magma",
        #   direction = 1,
        #   begin = 0,
        #   end = 1,
        #   name = "Gradient value",
        #   na.value = "grey85"
        # )
          scale_fill_gradient2(
            low = "#828F9D",
            mid = "white",
            high = "#CB915B",
            limits = GRADIENT_LIMITS,
            breaks = GRADIENT_BREAKS,
            oob = scales::squish,
            na.value = "grey85",
            name = "Gradient value"
          ) 
        
      }
    }
    
    plots <- list()

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
    ggsave(
      filename = paste0("analysis/gradients/plots/", plot_name, ".pdf"),
      plot = plot_output[[plot_group]][[plot_name]],
      width = 8,
      height = 6,
      dpi = 300
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