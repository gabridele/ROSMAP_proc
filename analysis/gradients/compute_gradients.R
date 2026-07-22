library(tidyverse)
library(broom)
library(SCORPIUS)
library(conflicted)
library(magick)
library(sf)
library(readxl)
library(ggpmisc)
library(cowplot)
library(pbapply)
library(RcppCNPy)
library(reticulate)
conflicts_prefer(ggpp::annotate)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::lag)

rois <- read.csv2(file.path("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/Schaefer2018_400Parcels_order.txt"), header = TRUE)
rois <- rois %>% rename(region = V2) 
# make character
rois$region <- as.character(rois$region)


# Source file for calculations related to deriving gradients
source("~/github_dir/ROSMAP_proc/analysis/gradients/util_gradients.R")
source("~/github_dir/ROSMAP_proc/analysis/gradients/util_vis.R")

##############
# Atlas setup
##############

# atlas_dir <- "data/atlas_data"

# schaef1k <- readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds"))
# schaef1k <- readRDS(file.path(atlas_dir, "schaef_ggseg2.rds"))
# yeo7 <- readRDS(file.path(atlas_dir, "Yeo2011_7_geometry.rds"))

# Yeo 7 network colors and names
net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                        col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), 
                        label = c(1:7))

# How they map to the Schaefer 1000 parcels
# GOTTA CREATE MYSELF
# ignore first row aka header
yeo_msk <- as.numeric(read_lines(file.path("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/org_mask_yeo_400.txt"), skip = 1))
# Schaefer 1000 parcel names
# read and delete superfluous cols


# Everything put together
roi_data <- data.frame(region = rois, 
                       yeo_label = yeo_msk) |> inner_join(net_names, join_by(yeo_label == label))

# Colors for the different ends of the gradients
gradient_cols <- data.frame(gradient1 = c("#3F596D", "#D38A4E"), 
                            gradient2 = c("#4682B4", "#781286"),  
                            gradient3 = c("#8A6081", "#738518"))


# if (create_brain_permutations) {
#   # Create the brain permutations to use for spintest, long compute time
  
#   # File containing an array of functions to facillitate calculations
#   source("src/util.R")
  
#   # Coords downloaded from https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/Centroid_coordinates
#   coords <- read_csv(file.path(atlas_dir, "Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv"))
  
#   perms <- rotate.parcellation(coord.l = coords[1:500, 3:5] |> as.matrix(), 
#                                coord.r = coords[501:1000, 3:5] |> as.matrix(), 
#                                nrot=1000, 
#                                method='hungarian')
#   write_rds(perms, file.path(atlas_dir, "permutations_1000_hungarian.rds"))
# } 


np <- import("numpy")

gradient_dir <- "gradients"
  

marg_gradients <- read_csv("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/marg/volumetric_marg_customgroup_400.csv") |> select(gradient1, gradient2, gradient3)
  
  ###################################################
  # Create gradients from average healthy connectomes
  ###################################################
  
#   healthy_young_connectomes <- con_cube_bf[, , biofinder_df |> filter(fmri_bl, age < 61, diagnosis == "Normal", abnorm_ab==0, !apoe4) |> pull(image_file)]
#   average_connectome <- apply(healthy_young_connectomes, c(1, 2), mean)
#   #write_rds(average_connectome, file.path(atlas_dir, "average_connectome_normalyoung.rds"))
#   average_connectome = read_rds(file.path(atlas_dir, "average_connectome_normalyoung.rds"))
#   rm(healthy_young_connectomes)
  
  # This takes some time, creates supplementary figure 9
  #gradient_comparison <- plot_grads_over_params(connectome_list = list(biofinder=average_connectome))
  
  #comp_plot <- ggdraw() +
#     draw_plot(gradient_comparison[["biofinder"]], x = 0, y = 0, width = 1, height = 0.95) +
#     draw_label(label = "PCA", x = 0.34, y = 0.985, size = 28) +
#     draw_label(label = "DME", x = 0.785, y = 0.985, size = 28) +
#     draw_label(label = "Threshold:", x = 0.34, y = 0.955, size = 18) +
#     draw_label(label = "Threshold:", x = 0.785, y = 0.955, size = 18) +
#     draw_line(x = c(0.15, 0.545), y = c(0.9725, 0.9725)) +
#     draw_line(x = c(0.59, 0.985), y = c(0.9725, 0.9725)) 
  
#   comp_plot <- comp_plot + theme(plot.background = element_rect(color = "black"))
#   # This is max millimiter (of journals) in inches
#   img_width = 180 / 25.4
  
#   # Scale image to balance figure elements and fontsize, but need "shrink" image after
#   # so the fontsize chosen should be large enough to shrink it by the scaling factor
#   scaling_factor <-  3
#   magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")
#   p_name <- "gradient_param_comparison_bf.png"
  
#   ggsave(file.path(figure_path, p_name), comp_plot , #patch_plots[["biofinder"]], 
#          width = img_width*scaling_factor, height = img_width*scaling_factor*0.675, bg ="white")
#   img <- magick::image_read(file.path(figure_path, p_name))
#   img_resized <- magick::image_resize(img, magick_geom_scaling)
#   magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)

# load my connectome
xcpd_avg_fc <- np$load("/Users/ga0034de/xcpd_456parcels/xcpd_avg_fc_matrix.npy", allow_pickle = TRUE)
custom_avg_fc <- np$load("/Users/ga0034de/custom_atlas_fc_matrices_442/customatlas400_avg_fc_matrix.npy", allow_pickle = TRUE)
individ_avg_fc <- np$load("/Users/ga0034de/individatlas_fc_matrices_442/individuatlas400_avg_fc_matrix.npy", allow_pickle = TRUE)
xcpd_avg_fc600 <- np$load("/Users/ga0034de/xcpd_600parcels/600avg_fc_matrix.npy", allow_pickle = TRUE)

#site_wise
bnk400xcpd <- np$load("/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/bnk/avg_bnk_fc_matrix.npy", allow_pickle = TRUE)
uc400xcpd <- np$load("/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/uc/avg_uc_fc_matrix.npy", allow_pickle = TRUE)
mg400xcpd <- np$load("/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/mg/avg_mg_fc_matrix.npy", allow_pickle = TRUE)
rirc400xcpd <- np$load("/Users/ga0034de/Desktop/timeseries_2306_xcpd400/site_wise/FC/rirc/avg_rirc_fc_matrix.npy", allow_pickle = TRUE)

bnk600xcpd <- np$load("/Users/ga0034de/Desktop/ts_656/site_wise/FC/bnk/avg_bnk_fc_matrix.npy", allow_pickle = TRUE)
uc600xcpd <- np$load("/Users/ga0034de/Desktop/ts_656/site_wise/FC/uc/avg_uc_fc_matrix.npy", allow_pickle = TRUE)
mg600xcpd <- np$load("/Users/ga0034de/Desktop/ts_656/site_wise/FC/mg/avg_mg_fc_matrix.npy", allow_pickle = TRUE)
rirc600xcpd <- np$load("/Users/ga0034de/Desktop/ts_656/site_wise/FC/rirc/avg_rirc_fc_matrix.npy", allow_pickle = TRUE)

bnk400custom <- np$load("/Users/ga0034de/Desktop/extracted_ts_customgroup/site_wise/FC/bnk/avg_fc_matrix.npy", allow_pickle = TRUE)
uc400custom <- np$load("/Users/ga0034de/Desktop/extracted_ts_customgroup/site_wise/FC/uc/avg_fc_matrix.npy", allow_pickle = TRUE)
mg400custom <- np$load("/Users/ga0034de/Desktop/extracted_ts_customgroup/site_wise/FC/mg/avg_fc_matrix.npy", allow_pickle = TRUE)
rirc400custom <- np$load("/Users/ga0034de/Desktop/extracted_ts_customgroup/site_wise/FC/rirc/avg_fc_matrix.npy", allow_pickle = TRUE)

bnk400invid <- np$load("/Users/ga0034de/Desktop/output_ts_individualatlas_SEND/site_wise/FC/bnk/avg_fc_matrix.npy", allow_pickle = TRUE)
uc400invid <- np$load("/Users/ga0034de/Desktop/output_ts_individualatlas_SEND/site_wise/FC/uc/avg_fc_matrix.npy", allow_pickle = TRUE)
mg400invid <- np$load("/Users/ga0034de/Desktop/output_ts_individualatlas_SEND/site_wise/FC/mg/avg_fc_matrix.npy", allow_pickle = TRUE)
rirc400invid <- np$load("/Users/ga0034de/Desktop/output_ts_individualatlas_SEND/site_wise/FC/rirc/avg_fc_matrix.npy", allow_pickle = TRUE)

grad_list_individ <- get_gradients(connectome_ests = list(bnk400invid = bnk400invid, uc400invid = uc400invid, mg400invid = mg400invid, rirc400invid = rirc400invid),
                                reference_gradients = marg_gradients,
                                n_gradients = c(1,2,3),
                                threshold = 0.0,
                                similarity_method = "cosine",
                                on_affinity = FALSE,
                                method = "pca",
                                visualize = TRUE)
  
  
  params <- expand_grid(method = c("pca", "diffusion"), affinity = c(FALSE, TRUE), sim_method = "cosine", threshold = c(0.0, 0.25, 0.5, 0.75)) |> 
    filter(!(method == "diffusion" & !affinity)) |> 
    filter(!(method == "pca" & affinity)) |> 
    mutate(sim_method = ifelse(!affinity, NA, sim_method))
  
  print("Calculating gradients over various parameters")
  gradient_data <- c()
  varexp_df <- c()
  for (i in 1:nrow(params)) {
    param_i <- params[i, ]
    grad_list <- get_gradients(connectome_ests = list(xcpd600 = xcpd_avg_fc600),
                               n_gradients = c(1,2,3),
                               threshold = param_i$threshold,
                               similarity_method = param_i$sim_method,
                               on_affinity = param_i$affinity,
                               method = param_i$method,
                               visualize = TRUE,
                               side_density = FALSE)
    
    gradient_data <- rbind(gradient_data, grad_list$gradients)
    varexp_df <- rbind(varexp_df, grad_list$varexp)
  }
  
  gradient_data <- gradient_data |> distinct()
