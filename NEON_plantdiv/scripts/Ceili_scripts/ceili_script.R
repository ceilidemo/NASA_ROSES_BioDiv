#### Compile plant coer data myself using new plant code lookup... 
#### Load libraries
library(dplyr)
library(purrr)
library(V.PhyloMaker)
library(stringr)

###########################
#### Create NEON plant code lookup
dat <- read.csv("data_in/allyr/div_1m2Data.csv")

lookup <- dat %>%
  select(taxonID, scientificName, taxonRank, family, identificationQualifier, morphospeciesID) %>%
  filter(!is.na(taxonID) & taxonID != "") %>%   # keep only plant entries
  distinct() %>%
  mutate(
    phyloName = case_when(
      taxonRank == "species" & !is.na(scientificName) & scientificName != "" ~
        str_replace_all(scientificName, "\\s+", "_"),
      TRUE ~ taxonID
    )
  )

write.csv(lookup, "data_in/lookup/NEON_plantcodes_lookup.csv", row.names = FALSE)

###########################
#### Compile plant cover data

look <- read.csv("data_in/lookup/NEON_plantcodes_lookup.csv", as.is = TRUE)
fls <- list.files("data_out/plantdiv/", full.names = TRUE)
plants_all <- map(fls, read.csv, check.names = FALSE)

# Aggregate by plot
df_list <- map(plants_all, ~ {
  .x %>%
    mutate(plotID = substr(subID, 1, 13)) %>%
    select(plotID, everything()) %>%
    group_by(plotID) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
})

df <- bind_rows(df_list)
df[is.na(df)] <- 0

###########################
#### Clean column names
df_clean <- df

# Remove "cf.", "aff.", "spp."
replace_patterns <- c("cf\\.\\s*.*" = "", "aff\\.\\s*.*" = "", "spp\\.?$" = "")
for(pat in names(replace_patterns)){
  colnames(df_clean) <- gsub(pat, replace_patterns[pat], colnames(df_clean))
}

# Remove "/" and "." in names
colnames(df_clean) <- sapply(strsplit(colnames(df_clean), "/"), `[`, 1)
colnames(df_clean) <- gsub("\\.", "_", colnames(df_clean))

# Manual rename based on lookup
manual_rename <- look %>% filter(!is.na(taxonID)) %>%
  select(taxonID, scientificName)

for(i in seq_len(nrow(manual_rename))){
  old <- manual_rename$taxonID[i]
  new <- gsub(" ", "_", manual_rename$scientificName[i])
  colnames(df_clean)[colnames(df_clean) == old] <- new
}

# Sum duplicate columns
dup_base <- str_replace(colnames(df_clean), "\\..*", "")
dup_groups <- split(colnames(df_clean), dup_base)

df_summed <- lapply(dup_groups, function(cols){
  if("plotID" %in% cols) return(df_clean[, "plotID", drop=FALSE])
  rowSums(df_clean[, cols], na.rm = TRUE) %>% as.data.frame()
})

# Name columns properly and bind
for(i in seq_along(df_summed)){
  colnames(df_summed[[i]]) <- names(dup_groups)[i]
}
df_final <- bind_cols(df_summed) %>% select(plotID, sort(colnames(.)[-1]))

write.csv(df_final, "data_out/biodiv_compiled/herb_cover_allplots.csv", row.names = FALSE)

###########################
#### Prepare phylo-safe dataset
df <- read.csv("data_out/biodiv_compiled/herb_cover_allplots.csv", check.names = FALSE)

# Clean species names
species_cols <- colnames(df)[-1]
species_cols <- species_cols[!grepl("_NA$|Fungi", species_cols)]  # remove unknowns/fungi

clean_names <- species_cols %>%
  gsub("[()]", "", .) %>%
  gsub("__+", "_", .) %>%
  gsub("[/\\-]", "_", .) %>%
  gsub("\\s+", "_", .)

df_phylo <- df[, c("plotID", species_cols)]
colnames(df_phylo)[-1] <- clean_names

# Sum duplicates if any
dup_names <- colnames(df_phylo)[duplicated(colnames(df_phylo))]
for(d in unique(dup_names)){
  cols <- which(colnames(df_phylo) == d)
  df_phylo[, cols[1]] <- rowSums(df_phylo[, cols], na.rm = TRUE)
  df_phylo <- df_phylo[, -cols[-1]]
}

###########################
#### Generate phylogenetic tree
species_list <- data.frame(
  species = colnames(df_phylo)[-1]
) %>%
  mutate(genus = gsub("_.*", "", species)) %>%
  left_join(lookup %>% select(phyloName, family),
            by = c("species" = "phyloName"))

tre <- phylo.maker(
  sp.list = species_list,
  tree = GBOTB.extended,
  nodes = nodes.info.1,
  scenarios = "S3"
)

saveRDS(tre, "data_out/phylo/allplants_scenario3_tree.rds")
write.csv(df_phylo, "data_out/biodiv_compiled/herb_cover_allplots_phylomatch.csv", row.names = FALSE)


################
#Rough and dirty diversity
################

library(dplyr)
library(vegan)

df <- read.csv("data_out/biodiv_compiled/herb_cover_allplots_phylomatch.csv", check.names = FALSE)

# Extract site and year from plotID
# Assuming plotID format like "ABBY_2018_001"
df <- df %>%
  mutate(
    site = sub("_.*", "", plotID),
    year = sub(".*_(\\d{4})_.*", "\\1", plotID)
  )

species_data <- df %>% select(-plotID, -site, -year)

# Add species data back to main dataframe
df_long <- bind_cols(df %>% select(plotID, site, year), species_data)

species_cols <- setdiff(colnames(df_long), c("plotID", "site", "year"))

# Calculate alpha, gamma, beta per site-year
diversity_site_year <- df_long %>%
  group_by(site, year) %>%
  summarise(
    alpha_richness = mean(rowSums(across(all_of(species_cols)) > 0)),
    alpha_shannon = mean(apply(select(cur_data(), all_of(species_cols)), 1, vegan::diversity, index = "shannon")),
    gamma_richness = sum(colSums(select(cur_data(), all_of(species_cols))) > 0),
    gamma_shannon = vegan::diversity(colSums(select(cur_data(), all_of(species_cols))), index = "shannon"),
    beta_whittaker = gamma_richness / alpha_richness - 1,
    beta_bray = {
      mat <- as.matrix(select(cur_data(), all_of(species_cols)))
      if (nrow(mat) > 1) {
        mean(as.matrix(vegdist(mat, method = "bray")))
      } else {
        NA_real_  
      }
    },
    .groups = "drop"
  )


write.csv(diversity_site_year, "data_out/abg/rough_firstrun_notscaled_plantdiv.csv", row.names = FALSE)



#####################
# new woody position compilation
#####################

library(dplyr)
library(ggplot2)
library(ggforce)
library(tidyr)
library(vegan)
library(stringr)

vst_perplotperyear <- read.csv("data_in/allyr/vst_perplotperyear.csv")
vst_mappingandtagging <- read.csv("data_in/allyr/vst_mappingandtagging.csv")
vst_apparentindividual <- read.csv("data_in/allyr/vst_apparentindividual.csv")

# Get only distributed plots
plot_summary <- vst_perplotperyear %>%
  filter(str_detect(plotType, "distributed"))

# Filter mapped stems in those plots
stems_mapped <- vst_mappingandtagging %>%
  filter(plotID %in% plot_summary$plotID) %>%
  distinct(individualID, .keep_all = TRUE)

# Get structure info for those stems
vst_individual <- vst_apparentindividual %>%
  filter(plotID %in% plot_summary$plotID) %>%
  filter(individualID %in% stems_mapped$individualID) %>%
  mutate(date = as.Date(date)) %>%
  group_by(individualID) %>%
  slice_max(order_by = date, n = 1, with_ties = FALSE) %>%
  ungroup()

# Join mapping + structure
mapped_vst <- left_join(stems_mapped, vst_individual, by = "individualID")

# Compute coordinates relative to plot center
mapped_trees <- mapped_vst %>%
  filter(!is.na(height), !is.na(maxCrownDiameter)) %>%
  mutate(
    crown_radius = maxCrownDiameter / 2,
    azimuth_rad = stemAzimuth * pi / 180,
    x = stemDistance * sin(azimuth_rad),
    y = stemDistance * cos(azimuth_rad)
  ) %>%
  group_by(plotID.x) %>%
  arrange(height) %>%
  mutate(height_rank = row_number()) %>%
  ungroup()

mapped_trees <- mapped_trees %>%
  rename(plotID = plotID.x) %>%
  select(plotID, individualID, taxonID, height, crown_radius, x, y, everything())


# Example for a single plot
plot_id <- "DEJU_010"

ggplot(mapped_trees %>% filter(plotID == plot_id)) +
  geom_circle(aes(x0 = x, y0 = y, r = crown_radius, fill = taxonID), alpha = 0.4, color = "black") +
  coord_fixed() +
  theme_minimal() +
  labs(title = paste("Tree crowns in plot", plot_id),
       x = "x (m)", y = "y (m)")


shrubs_data <- read.csv("data_in/allyr/vst_shrub.csv")



##########
#Plotting sites EMIT with NEON
####

library(terra)
library(sf)

emit <- rast("boundary_sampling/2023/MOAB/hyperspectral_cube.tif")

plotRGB(emit, r = 43, g = 23, b = 10, stretch = "lin", main = "EMIT RGB")

NEON <- rast("NEON_specdiv/data_in/processed_plot_cubes/NEON_D12_YELL_DP3_539000_4979000_bidirectional_reflectance_BRDF_plot_YELL_018_processed.tif")

plot(NEON)
names(NEON)


