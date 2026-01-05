## 01_ProcessFiles: Process the downloaded data to get working files for diversity calculations

library(dplyr)
library(tidyr)

# figure which plots you want
plot_metadata <- structure_all$vst_perplotperyear %>%
  select(plotID, plotType) %>%
  distinct()

# process raw data for presence / percent cover to get the site, plot, and abiotic vs biotic
ground_df <- latest_cover_data %>%
  mutate(final_id = coalesce(taxonID, otherVariables)) %>% # combine the IDs (aka take taxon id but if gone then othervar)
  filter(!is.na(final_id)) %>% 
  group_by(siteID, plotID, final_id) %>%
  summarize(mean_pct = mean(percentCover, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    area_m2 = (mean_pct / 100) * 400,
    is_abiotic = grepl("litter|rock|soil|water|wood|other|scat|fungi|standing|biocrust|lichen|moss",  # pulled out some abiotic words i saw
                       tolower(final_id))
  ) %>% 
  rename(taxonID = final_id) %>% 
  mutate(source = "ground", functional_type = "Herbaceous")

# process the trees and shrubs 
taxon_map <- structure_all$vst_mappingandtagging %>%
  select(individualID, taxonID) %>%
  filter(!is.na(taxonID)) %>%
  distinct(individualID, .keep_all = TRUE)

trees_df <- structure_all$vst_apparentindividual %>%
  left_join(taxon_map, by = "individualID") %>%
  mutate(area_m2 = pi * (maxCrownDiameter / 2)^2) %>%
  group_by(siteID, plotID, taxonID) %>%
  summarize(area_m2 = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "canopy", functional_type = "Woody") %>% 
  mutate(is_abiotic = FALSE)

shrubs_df <- structure_all$vst_shrubgroup %>%
  group_by(siteID, plotID, taxonID) %>%
  summarize(area_m2 = sum(canopyArea, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "shrub", functional_type = "Woody") %>% 
  mutate(is_abiotic = FALSE)

# Merge them all 
merge_long <- bind_rows(ground_df, trees_df, shrubs_df) %>%
  left_join(plot_metadata, by = "plotID") %>%
  filter(plotType == "distributed" | is.na(plotType))

