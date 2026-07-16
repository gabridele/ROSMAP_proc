

library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(sf)

mask400 <- read.csv("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/Schaefer2018_400Parcels_7Networks_order.txt", header = FALSE, sep = "\t") %>% select(c(2))

write.csv(mask400, "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/Schaefer2018_400Parcels_order.txt", row.names = FALSE, quote = FALSE)

mask400 <- mask400 %>%
  mutate(network = case_when(
    str_detect(V2, "Vis") ~ 1,
    str_detect(V2, "SomMot") ~ 2,
    str_detect(V2, "DorsAttn") ~ 3,
    str_detect(V2, "SalVentAttn") ~ 4,
    str_detect(V2, "Limbic") ~ 5,
    str_detect(V2, "Cont") ~ 6,
    str_detect(V2, "Default") ~ 7
  ))

mask400 <- mask400 %>% select(network)



# save
write.csv(mask400, "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/org_mask_yeo_400.txt", row.names = FALSE, quote = FALSE)

mask600 <- read.csv("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/Schaefer2018_600Parcels_7Networks_order.txt", header = FALSE, sep = "\t") %>% select(c(2))

write.csv(mask600, "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/Schaefer2018_600Parcels_order.txt", row.names = FALSE, quote = FALSE)

mask600 <- mask600 %>%
  mutate(network = case_when(
    str_detect(V2, "Vis") ~ 1,
    str_detect(V2, "SomMot") ~ 2,
    str_detect(V2, "DorsAttn") ~ 3,
    str_detect(V2, "SalVentAttn") ~ 4,
    str_detect(V2, "Limbic") ~ 5,
    str_detect(V2, "Cont") ~ 6,
    str_detect(V2, "Default") ~ 7
  ))

mask600 <- mask600 %>% select(network)

# save
write.csv(mask600, "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/org_mask_yeo_600.txt", row.names = FALSE, quote = FALSE)

# get ggseg2 atlas
PARCELS <- 600
atlas_fun <- get(
  paste0("schaefer7_", PARCELS),
  envir = asNamespace("ggsegSchaefer")
)

atlas <- atlas_fun()

atlas_full <- ggseg.formats::atlas_sf(atlas)

names(atlas_full)
head(atlas_full)
atlas_full %>%
  filter(is.na(region)) %>%
  select(label, hemi, view, colour, geometry)

atlas_full %>%
  filter(is.na(region) | region == "") %>%
  st_drop_geometry()

atlas_df <- atlas_full %>%
  filter(!is.na(region)) %>%
  transmute(
    region,
    hemisphere = hemi,
    view,
    color = colour,
    geometry
  )

shade_df <- atlas_full %>%
  filter(is.na(region)) %>%
  transmute(
    hemisphere = hemi,
    view,
    geometry
  )

silhouette_df <- atlas_full %>%
  group_by(hemi, view) %>%
  summarise(
    geometry = st_union(geometry),
    .groups = "drop"
  ) %>%
  mutate(
    group_id = paste(hemi, view, sep = "_"),
    geom_type = as.character(st_geometry_type(geometry)),
    geom_length = vapply(
      geometry,
      function(x) nrow(st_coordinates(x)),
      integer(1)
    ),
    hemisphere = hemi,
    int_view = match(view, unique(view)),
    hemi_n = match(hemisphere, unique(hemisphere)),
    view_n = match(view, unique(view)),
    idx = row_number(),
    row = hemi_n,
    col = view_n,
    shift_x = 0,
    shift_y = 0
  ) %>%
  select(
    group_id,
    geometry,
    geom_type,
    geom_length,
    hemisphere,
    view,
    int_view,
    hemi_n,
    view_n,
    idx,
    row,
    col,
    shift_x,
    shift_y
  )
atlas_geometry <- list(
  atlas = atlas_df,
  silhouette = silhouette_df,
  shade = shade_df
)
names(atlas_geometry)

str(atlas_geometry, max.level = 1)

n_parcels <- atlas_df %>%
  filter(!is.na(region), region != "") %>%
  distinct(region) %>%
  nrow()
print(paste0("Number of parcels in atlas: ", n_parcels))
#save to rds
saveRDS(atlas_geometry, paste0("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/schaef", PARCELS, "_ggseg2.rds"))

