## 03_TopDown_PlantDiv: Fixed with latest-visit filtering and geometric clipping
library(dplyr)
library(tidyr)
library(sf)
library(purrr)
library(ggplot2)

#######################
# function for dealing with overlapping crowns below
clip_overlapping_crowns <- function(df) {
  if (nrow(df) == 0) return(NULL)
  
  visible_geoms <- list()
  covered_area <- st_sfc(crs = st_crs(df))
  
  for (i in 1:nrow(df)) {
    current_tree <- df[i, ]
    geom <- st_make_valid(current_tree$geometry)
    geom <- st_intersection(geom, plot_bbox)
    
    if (length(covered_area) == 0) {
      visible_part <- geom
    } else {
      mask <- st_make_valid(st_union(covered_area))
      visible_part <- st_difference(geom, mask)
    }
    
    if (length(visible_part) > 0 && !st_is_empty(visible_part)) {
      current_tree$geometry <- visible_part
      current_tree$visible_area <- as.numeric(st_area(visible_part))
      
      if (current_tree$visible_area > 0.01) {
        visible_geoms[[length(visible_geoms) + 1]] <- current_tree
        covered_area <- st_make_valid(st_union(covered_area, visible_part))
      }
    }
  }
  
  if (length(visible_geoms) == 0) return(NULL)
  return(bind_rows(visible_geoms))
}

#######################
# Data clean and set up mapping 
# Took all the str since they don't do it every year... so took just the latest year/ most recent
latest_structure <- structure_all$vst_apparentindividual %>%
  group_by(individualID) %>% filter(date == max(date, na.rm = TRUE)) %>% slice(1) %>% ungroup()

# get the mapped locations of the trees
latest_mapping <- structure_all$vst_mappingandtagging %>% distinct(individualID, .keep_all = TRUE)

# get relative location in the plots (NEON does it by azimuth)
mapped_trees <- latest_mapping %>%
  select(individualID, plotID, stemDistance, stemAzimuth, taxonID) %>%
  inner_join(latest_structure %>% select(individualID, height, maxCrownDiameter), by = "individualID") %>%
  mutate(siteID = substr(plotID, 1, 4), crown_radius = maxCrownDiameter / 2,
         theta = stemAzimuth * pi / 180, x = stemDistance * sin(theta), y = stemDistance * cos(theta)) %>%
  filter(!is.na(x), !is.na(y), !is.na(crown_radius))

# create poly around location and crown radius
canopy_polys <- mapped_trees %>% st_as_sf(coords = c("x", "y")) %>% st_buffer(dist = .$crown_radius)
plot_bbox <- st_polygon(list(matrix(c(-10,-10, 10,-10, 10,10, -10,10, -10,-10), ncol=2, byrow=TRUE))) %>% 
  st_sfc() %>% st_set_crs(st_crs(canopy_polys))

# arrange by height and then deal with the overlapping crowns
canopy_visible_final <- canopy_polys %>% arrange(plotID, desc(height)) %>%
  group_by(plotID) %>% group_split() %>%
  map_dfr(clip_overlapping_crowns) 

# Calculate how much of the plot area is NOT covered by thee trees 
gap_summary <- canopy_visible_final %>% st_drop_geometry() %>%
  group_by(plotID) %>% summarize(total_visible_canopy = sum(visible_area, na.rm = TRUE)) %>%
  mutate(gap_fraction = pmax(0, (400 - total_visible_canopy) / 400))

# merge all the layers (vis trees, shrub, herb) and scale by gap calculated 
topdown_master_table <- bind_rows(
  canopy_visible_final %>% st_drop_geometry() %>% select(siteID, plotID, taxonID, area_m2 = visible_area),
  latest_shrubgroups %>% select(siteID, plotID, taxonID, area_m2),
  ground_df %>% left_join(gap_summary, by = "plotID") %>% 
    mutate(area_m2 = (mean_pct / 100) * (400 * replace_na(gap_fraction, 1)))
) %>% filter(area_m2 > 0)

##############################
# Calculate div whoop whoop 

# join NEon names to the phylo 
topdown_phylo_ready <- topdown_master_table %>%
  left_join(name_map, by = "taxonID") %>%
  mutate(taxonID = coalesce(phyloName, taxonID))

topdown_matrix <- topdown_phylo_ready %>%
  group_by(siteID, taxonID) %>% summarize(total_area = sum(area_m2), .groups = "drop") %>%
  pivot_wider(names_from = taxonID, values_from = total_area, values_fill = 0) %>%
  column_to_rownames("siteID")

# Get hill number calcs for TD
td_hill_results <- data.frame(siteID = rownames(topdown_matrix)) %>%
  mutate(
    TD_Richness = specnumber(topdown_matrix),
    TD_Shannon_Eff = exp(diversity(topdown_matrix, index = "shannon")),
    TD_Simpson_Eff = 1/diversity(topdown_matrix, index = "simpson")
  )

# Get PD Rao Q 
pd_rao_results <- topdown_phylo_ready %>%
  group_by(siteID) %>%
  group_modify(~ {
    calc_rao_partition(.x %>% rename(total_area = area_m2), dist_phylo)
  })

# TD beta attempt (idk if this makes complete sense, but looking at avg difference btw plots on what's visible)
td_beta_values <- topdown_phylo_ready %>%
  group_by(siteID) %>%
  group_modify(~ {
    site_mat <- .x %>% group_by(plotID, taxonID) %>% summarize(area = sum(area_m2), .groups = "drop") %>%
      pivot_wider(names_from = taxonID, values_from = area, values_fill = 0) %>% column_to_rownames("plotID")
    if(nrow(site_mat) > 1) return(data.frame(TD_Beta = mean(vegdist(site_mat, method = "bray"))))
    return(data.frame(TD_Beta = 0))
  })

############################################
# Joiin results and save
final_topdown_diversity <- td_hill_results %>% 
  left_join(td_beta_values, by = "siteID") %>%
  left_join(pd_rao_results, by = "siteID") %>%
  mutate(across(where(is.numeric), ~round(., 4)))

write.csv(final_topdown_diversity, "data_out/NEON_TopDown_Diversity_Results.csv", row.names = FALSE)



############################################
############################################
############################################
# PLOT CHECK IF MAKES SENSE 

site_to_check <- "MOAB" 

# Filter data to that site
site_map_data <- canopy_visible_final %>% 
  filter(grepl(site_to_check, plotID))

# plot all the plots for that site with the mapped crowns
ggplot() +
  geom_sf(data = plot_bbox, fill = "gray95", color = "black", linetype = "dashed") +
  geom_sf(data = site_map_data, aes(fill = taxonID), alpha = 0.8, color = "white", size = 0.1) +
  facet_wrap(~plotID, ncol = 8) + 
  scale_fill_viridis_d(option = "plasma") + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = paste("'Top-Down' View of plots in:", site_to_check),
    subtitle = "Visible mapped canopy in each plot w/species ID",
    fill = "Taxon ID",
    x = NULL, y = NULL
  )
