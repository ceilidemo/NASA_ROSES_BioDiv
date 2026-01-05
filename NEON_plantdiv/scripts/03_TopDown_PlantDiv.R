## 03_TopDown_PlantDiv: Plant diversity partitioning at the site level with a top-down view (scaled abundance)

library(dplyr)
library(tidyr)
library(sf)

# Get mapped trees
mapped_trees <- structure_2023$vst_mappingandtagging %>%
  select(individualID, plotID, stemDistance, stemAzimuth, taxonID) %>%
  left_join(
    structure_2023$vst_apparentindividual %>%
      select(individualID, height, maxCrownDiameter, date),
    by = "individualID"
  ) %>%
  mutate(
    crown_radius = maxCrownDiameter / 2,
    theta = stemAzimuth * pi / 180,
    x = stemDistance * sin(theta),
    y = stemDistance * cos(theta)
  ) %>%
  filter(!is.na(x), !is.na(y), !is.na(crown_radius))

# Create circles for every tree crown
canopy_polys <- mapped_trees %>%
  st_as_sf(coords = c("x", "y")) %>%
  st_buffer(dist = .$crown_radius)

# For each plot, dissolve to get canopy foot
canopy_summary <- canopy_polys %>%
  group_by(plotID) %>%
  summarize(geometry = st_union(geometry)) %>%
  mutate(
    canopy_area = as.numeric(st_area(geometry)),
    # Scale to 400m2 plot, ensuring we don't exceed 100%
    gap_area = pmax(0, 400 - canopy_area),
    gap_fraction = gap_area / 400
  )

# Join gap fraction back to your ground data
ground_scaled <- ground_df %>%
  left_join(st_drop_geometry(canopy_summary), by = "plotID") %>%
  mutate(
    # If no canopy data exists for a plot, assume 100% gap (gap_fraction = 1)
    gap_fraction = ifelse(is.na(gap_fraction), 1, gap_fraction),
    area_m2_scaled = (mean_pct / 100) * (400 * gap_fraction)
  )

# Combine with woody data (using actual crown areas)
layered_scaled_final <- bind_rows(
  ground_scaled %>% 
    select(siteID, plotID, taxonID, area_m2 = area_m2_scaled) %>% 
    mutate(Source = "Scaled Herbaceous"),
  trees_df %>% 
    mutate(Source = "Direct Woody")
)