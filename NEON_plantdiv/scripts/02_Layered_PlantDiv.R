## 02_Layered_PlantDiv: Standardized Naming and Hill Number Conversion
library(vegan)
library(tidyverse)
library(FD)
library(V.PhyloMaker2)

# Function for abundance-weighted TD: 
calc_abundance_div <- function(df) {
  comm <- df %>% # comm matrix
    group_by(plotID, taxonID) %>%
    summarize(total_area = sum(total_area, na.rm = TRUE), .groups = "drop") %>% # summ tot area for each spc in each plot
    pivot_wider(names_from = taxonID, values_from = total_area, values_fill = 0) %>% # row -> plot, col -> species
    column_to_rownames("plotID")
  
  site_totals <- colSums(comm) # just to grab tot abundance 
  
  # Diversity calcs
  richness <- specnumber(site_totals) # count unique spec
  gamma_shannon_eff <- exp(diversity(site_totals, index = "shannon")) # shannon: weighted common spec (converted to eff)
  gamma_simpson_eff <- 1/diversity(site_totals, index = "simpson") # simpson: weighted dom spec (convert to eff)
  
  # Beta div ! (Bray-Curtis)
  if(nrow(comm) > 1) {
    beta_bray <- mean(vegdist(comm, method = "bray"), na.rm = TRUE)
  } else {
    beta_bray <- NA # if only 1 plot 
  }
  
  return(data.frame(
    TD_Richness = richness,
    TD_Shannon_Eff = gamma_shannon_eff,
    TD_Simpson_Eff = gamma_simpson_eff,
    TD_Beta = beta_bray
  ))
}

#######################
# Clean data 
layered <- merge_long %>%
  group_by(siteID, plotID, taxonID, is_abiotic) %>% 
  summarize(total_area = sum(area_m2, na.rm = TRUE), .groups = "drop") %>% 
  filter(total_area < 1000) # filter for any random mis-entries like a HUGE tree

layered_veg_results <- layered %>%
  filter(!is_abiotic) %>% # get rid of the abiotic tags
  group_by(siteID) %>%
  group_modify(~ calc_abundance_div(.x)) %>%
  mutate(Scope = "Veg Only")

#######################
# Function for PD (Rao's Q) 
calc_rao_partition <- function(df, phylo_dist) {
  comm <- df %>%
    group_by(plotID, taxonID) %>%
    summarize(total_area = sum(total_area, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxonID, values_from = total_area, values_fill = 0) %>%
    column_to_rownames("plotID")
  
  common_sp <- intersect(colnames(comm), rownames(phylo_dist)) #match the spec in site with those in tree
  if(length(common_sp) == 0) return(data.frame(PD_Alpha = NA, PD_Beta = NA, PD_Gamma = NA))
  
  comm_matched <- as.matrix(comm[, common_sp, drop = FALSE])
  dist_matched <- as.matrix(phylo_dist[common_sp, common_sp])
  
  # avg evol dist for each plto 
  rao_alphas <- apply(comm_matched, 1, function(x) { 
    if(sum(x) <= 0) return(0)
    p <- x / sum(x)
    return(sum(outer(p, p) * dist_matched)) # weighted dist btw pairs 
  })
  
  # gamma totals across site
  site_totals <- colSums(comm_matched)
  p_gamma <- site_totals / sum(site_totals)
  gamma_rao <- sum(outer(p_gamma, p_gamma) * dist_matched)
  
  # beta -> turnover (doing in additive way with gamma-alpha) 
  return(data.frame(
    PD_Alpha = mean(rao_alphas, na.rm = TRUE),
    PD_Beta = gamma_rao - mean(rao_alphas, na.rm = TRUE),
    PD_Gamma = gamma_rao
  ))
}

#######################
# Generate phylo tree
lookup <- cover_all$div_1m2Data %>%
  select(taxonID, scientificName, taxonRank, family) %>%
  filter(!is.na(taxonID) & taxonID != "") %>%
  distinct() %>%
  mutate(phyloName = case_when(
    taxonRank == "species" ~ str_replace_all(scientificName, "\\s+", "_"),
    TRUE ~ taxonID
  ))

name_map <- lookup %>% select(taxonID, phyloName) %>% distinct()

layered_phylo_ready <- layered %>%
  filter(!is_abiotic) %>% # get rid abiotic 
  left_join(name_map, by = "taxonID") %>%
  mutate(taxonID = coalesce(phyloName, taxonID))

sp_list <- layered_phylo_ready %>% # get just spc list for tree
  distinct(taxonID) %>%
  rename(species = taxonID) %>%
  mutate(genus = gsub("_.*", "", species)) %>%
  left_join(lookup %>% select(phyloName, family) %>% distinct(), by = c("species" = "phyloName")) %>%
  filter(!is.na(family))

data("nodes.info.1.TPL", package = "V.PhyloMaker2")
sp_list_corrected <- sp_list %>% # Correct the names with the ones that are in v.phylomaker2
  mutate(family_new = nodes.info.1.TPL$family[match(genus, nodes.info.1.TPL$genus)]) %>%
  mutate(family = coalesce(family_new, family)) %>%
  select(-family_new)

tre <- phylo.maker(sp.list = sp_list_corrected, tree = GBOTB.extended.TPL,  # generate tree using basic GBOTB backbone
                   nodes = nodes.info.1.TPL, scenarios = "S3")
dist_phylo <- as.matrix(cophenetic(tre$scenario.3)) # get the ditance between each in the tree 

#######################
# get PD for each site 
layered_rao_results <- layered_phylo_ready %>%
  group_by(siteID) %>%
  group_modify(~ calc_rao_partition(.x, dist_phylo))

#join TD with PD
layered_master_table <- layered_veg_results %>%
  left_join(layered_rao_results, by = "siteID") %>%
  select(siteID, Scope, starts_with("TD_"), starts_with("PD_")) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# save it 
write.csv(layered_master_table, "data_out/NEON_Layered_Diversity_Results.csv", row.names = FALSE)
