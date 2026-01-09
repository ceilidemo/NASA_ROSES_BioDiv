## 01_ProcessFiles: Process the downloaded data to get working files for diversity calculations

library(dplyr)
library(tidyr)

load("data_in/veg_data.RData")

# get the taxa from all the products
taxon_lookup <- bind_rows(
  cover_all$div_1m2Data %>% select(taxonID, scientificName, family, taxonRank),
  structure_all$vst_mappingandtagging %>% select(taxonID, scientificName, family, taxonRank)
) %>% 
  distinct(taxonID, .keep_all = TRUE) %>%
  filter(!is.na(taxonID)) %>%
  mutate(phyloName = case_when(
    taxonRank == "species" ~ gsub(" ", "_", scientificName),
    TRUE ~ taxonID
  ))

# process the ground cover data (getting most recent visit)
ground_df <- cover_all$div_1m2Data %>%
  group_by(plotID, subplotID) %>%
  slice_max(endDate, n = 1, with_ties = FALSE) %>% # Latest per subplot
  ungroup() %>%
  mutate(final_id = coalesce(taxonID, otherVariables)) %>%
  filter(!is.na(final_id)) %>%
  group_by(siteID, plotID, final_id) %>%
  summarize(mean_pct = mean(percentCover, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    is_abiotic = grepl("litter|rock|soil|water|wood|other|scat|fungi|standing|biocrust|lichen|moss", 
                       tolower(final_id)),
    source = "ground"
  ) %>%
  rename(taxonID = final_id)

# now the shrubs (getting most recent visit)
shrubs_df <- structure_all$vst_shrubgroup %>%
  group_by(plotID, groupID) %>%
  slice_max(date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(siteID, plotID, taxonID) %>%
  summarize(area_m2 = sum(canopyArea, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "shrub", is_abiotic = FALSE)

# now the mapped trees
taxon_map <- structure_all$vst_mappingandtagging %>%
  select(individualID, taxonID) %>%
  distinct(individualID, .keep_all = TRUE)

trees_df <- structure_all$vst_apparentindividual %>%
  group_by(individualID) %>%
  slice_max(date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(taxon_map, by = "individualID") %>%
  filter(!is.na(maxCrownDiameter)) %>%
  mutate(area_m2 = pi * (maxCrownDiameter / 2)^2,
         source = "canopy", 
         is_abiotic = FALSE)

# and then the plot metadata so you can filter by distributed baseplots
plot_metadata <- structure_all$vst_perplotperyear %>% 
  select(plotID, plotType) %>% distinct()

# save this compilation
save(ground_df, shrubs_df, trees_df, taxon_lookup, plot_metadata, structure_all,
     file = "data_out/community.RData")

