## 02_Layered_PlantDiv: "Layered" diversity calculation taking footprints from various levels of canopy

library(vegan)
library(tidyverse)
library(V.PhyloMaker2)

# load community data if not
load("data_out/community.RData")

# Diversity functions:

## Taxonomic (Hilll Numbers)
calc_t_div <- function(df) {
  comm <- df %>%
    group_by(plotID, taxonID) %>%
    summarize(total_area = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxonID, values_from = total_area, values_fill = 0) %>%
    column_to_rownames("plotID")
  
  site_totals <- colSums(comm)
  
  richness <- specnumber(site_totals)
  shannon_eff <- exp(diversity(site_totals, index = "shannon"))
  simpson_eff <- 1/diversity(site_totals, index = "simpson")
  
  beta_bray <- if(nrow(comm) > 1) mean(vegdist(comm, method = "bray"), na.rm = TRUE) else NA
  
  return(data.frame(TD_Richness = richness, TD_Shannon_Eff = shannon_eff, 
                    TD_Simpson_Eff = simpson_eff, TD_Beta = beta_bray))
}

## Phylo (Rao's Q Partitioning)
calc_p_div <- function(df, phylo_dist) {
  comm <- df %>%
    group_by(plotID, taxonID) %>%
    summarize(total_area = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxonID, values_from = total_area, values_fill = 0) %>%
    column_to_rownames("plotID")
  
  common_sp <- intersect(colnames(comm), rownames(phylo_dist))
  if(length(common_sp) == 0) return(data.frame(PD_Alpha = NA, PD_Beta = NA, PD_Gamma = NA))
  
  comm_matched <- as.matrix(comm[, common_sp, drop = FALSE])
  dist_matched <- as.matrix(phylo_dist[common_sp, common_sp])
  
  rao_alphas <- apply(comm_matched, 1, function(x) { 
    if(sum(x) <= 0) return(0)
    p <- x / sum(x)
    return(sum(outer(p, p) * dist_matched))
  })
  
  site_totals <- colSums(comm_matched)
  p_gamma <- site_totals / sum(site_totals)
  gamma_rao <- sum(outer(p_gamma, p_gamma) * dist_matched)
  
  return(data.frame(PD_Alpha = mean(rao_alphas, na.rm = TRUE), 
                    PD_Beta = gamma_rao - mean(rao_alphas, na.rm = TRUE), 
                    PD_Gamma = gamma_rao))
}

# Generate Phylogeny
sp_list <- taxon_lookup %>%
  filter(taxonRank == "species") %>%
  rename(species = phyloName) %>%
  mutate(genus = gsub("_.*", "", species)) %>%
  distinct(species, .keep_all = TRUE)

data("nodes.info.1.TPL", package = "V.PhyloMaker2")
sp_list_corrected <- sp_list %>%
  mutate(family_new = nodes.info.1.TPL$family[match(genus, nodes.info.1.TPL$genus)]) %>%
  mutate(family = coalesce(family_new, family)) %>%
  select(species, genus, family)

tre <- phylo.maker(sp.list = sp_list_corrected, tree = GBOTB.extended.TPL, 
                   nodes = nodes.info.1.TPL, scenarios = "S3")
dist_phylo <- as.matrix(cophenetic(tre$scenario.3))

# Get the names to match the phylogeny for calcs
layered_ready <- bind_rows(
  ground_df %>% mutate(area_m2 = (mean_pct / 100) * 400),
  shrubs_df,
  trees_df %>% select(siteID, plotID, taxonID, area_m2, is_abiotic)
) %>%
  filter(!is_abiotic) %>%
  left_join(taxon_lookup %>% select(taxonID, phyloName), by = "taxonID") %>%
  mutate(taxonID = coalesce(phyloName, taxonID)) %>%
  inner_join(plot_metadata %>% filter(plotType == "distributed") %>% select(plotID), by = "plotID") %>%
  group_by(siteID, plotID, taxonID) %>%
  summarize(area_m2 = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
  mutate(area_m2 = pmin(area_m2, 400))
  
# ok now actually calculate div
layered_td <- layered_ready %>%
  group_by(siteID) %>%
  group_modify(~ calc_t_div(.x))

layered_pd <- layered_ready %>%
  group_by(siteID) %>%
  group_modify(~ calc_p_div(.x, dist_phylo))

# join and save it
layered_master_table <- layered_td %>%
  left_join(layered_pd, by = "siteID") %>%
  mutate(across(where(is.numeric), ~round(., 4)))

write.csv(layered_master_table, "data_out/NEON_Layered_Diversity_Results.csv", row.names = FALSE)

# save phylogeny to use for the top down view :) 
save(dist_phylo, taxon_lookup, file = "data_out/phylo.RData")
