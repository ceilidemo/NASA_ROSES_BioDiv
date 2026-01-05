## 02_Layered_PlantDiv : Plant diversity partitioning at the site level with a 3 dimensional "layered" view

library(vegan)
library(tidyverse)
library(FD)
library(V.PhyloMaker2)

# calculating abundance-weighted diversity: 
calc_abundance_div <- function(df) {
  comm <- df %>%
    select(plotID, taxonID, total_area) %>%
    pivot_wider(names_from = taxonID, values_from = total_area, values_fill = 0) %>%
    column_to_rownames("plotID")
  
  # calc div for each plot, then take mean alpha for site
  alpha_shannon <- mean(diversity(comm, index = "shannon"))
  alpha_simpson <- mean(diversity(comm, index = "simpson"))
  
  # sum across sites
  site_totals <- colSums(comm)
  gamma_shannon <- diversity(site_totals, index = "shannon")
  gamma_simpson <- diversity(site_totals, index = "simpson")
  
  # beta div bray-curtis distance 
  if(nrow(comm) > 1) {
    dist_matrix <- vegdist(comm, method = "bray")
    beta_bray <- mean(dist_matrix)
  } else {
    beta_bray <- NA
  }
  
  return(data.frame(
    Alpha_Shannon = alpha_shannon,
    Alpha_Simpson = alpha_simpson,
    Gamma_Shannon = gamma_shannon,
    Beta_Distance = beta_bray
  ))
}

###############################################################
# "Layered" diversity calculations: BASIC abundance weighted...
layered <- merge_long %>%
  group_by(siteID, plotID, taxonID, is_abiotic) %>% # sum across all "layers"
  summarize(total_area = sum(area_m2, na.rm = TRUE), .groups = "drop")

# veg only
layered_veg_results <- layered %>%
  filter(!is_abiotic) %>%
  group_by(siteID) %>%
  group_modify(~ calc_abundance_div(.x)) %>%
  mutate(Scope = "Veg Only")

# full surface
layered_full_results <- layered %>%
  group_by(siteID) %>%
  group_modify(~ calc_abundance_div(.x)) %>%
  mutate(Scope = "Full Surface")

# comb
layered_abundance_diversity_table <- bind_rows(layered_veg_results, layered_full_results)


#############################################################
## Incorporating phylogenetic distance metrics ...
# Rao Q div function
calc_rao_partition <- function(df, phylo_dist) {
  # Aggregate to ensure one entry per taxon per plot
  comm_long <- df %>%
    group_by(plotID, taxonID) %>%
    summarize(total_area = sum(total_area, na.rm = TRUE), .groups = "drop")
  
  comm <- comm_long %>%
    pivot_wider(names_from = taxonID, values_from = total_area, values_fill = 0) %>%
    column_to_rownames("plotID")
  
  # Species matching
  common_sp <- intersect(colnames(comm), rownames(phylo_dist))
  
  # Handle sites with no phylogeny matches
  if(length(common_sp) == 0) {
    return(data.frame(Phylo_Alpha = NA, Phylo_Beta = NA, Phylo_Gamma = NA))
  }
  
  comm_matched <- as.matrix(comm[, common_sp, drop = FALSE])
  dist_matched <- as.matrix(phylo_dist[common_sp, common_sp])
  
  # Calculate Alpha 
  rao_alphas <- apply(comm_matched, 1, function(x) {
    if(sum(x) <= 0) return(0)
    p <- x / sum(x)
    # Rao formula: sum(pi * pj * dij)
    return(sum(outer(p, p) * dist_matched))
  })
  alpha_rao <- mean(rao_alphas, na.rm = TRUE)
  
  # Calculate Gamma 
  site_totals <- colSums(comm_matched)
  if(sum(site_totals) <= 0) {
    gamma_rao <- 0
  } else {
    p_gamma <- site_totals / sum(site_totals)
    gamma_rao <- sum(outer(p_gamma, p_gamma) * dist_matched)
  }
  
  # Beta = Gamma - Alpha 
  beta_rao <- gamma_rao - alpha_rao
  
  return(data.frame(
    Phylo_Alpha = alpha_rao,
    Phylo_Beta = beta_rao,
    Phylo_Gamma = gamma_rao
  ))
}

# create code lookup
lookup <- cover_all$div_1m2Data %>%
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
# clean up to use with phylogeny
name_map <- lookup %>% 
  select(taxonID, phyloName) %>% 
  distinct()

# get species list ready for phylogeny creation
layered_phylo_ready <- layered %>%
  filter(!is_abiotic) %>%
  left_join(name_map, by = "taxonID") %>%
  mutate(taxonID = coalesce(phyloName, taxonID)) %>% 
  select(siteID, plotID, taxonID, total_area)

sp_list_for_tree <- layered_phylo_ready %>%
  distinct(taxonID) %>%
  rename(species = taxonID) %>%
  mutate(genus = gsub("_.*", "", species)) %>%
  left_join(lookup %>% select(phyloName, family) %>% distinct(), 
            by = c("species" = "phyloName")) %>%
  select(species, genus, family) %>%
  filter(!is.na(family))

# Grab stuff (nodes) from phylomaker to fix the names
data("nodes.info.1.TPL", package = "V.PhyloMaker2")

# Create the correction list using the nodes info
sp_list_corrected <- sp_list_for_tree %>%
  mutate(family_new = nodes.info.1.TPL$family[match(genus, nodes.info.1.TPL$genus)]) %>%
  mutate(family = coalesce(family_new, family)) %>%
  select(-family_new)

# Check it with the acers where there are usually issues...
sp_list_corrected %>% filter(genus == "Acer")

# Generate that tree!
tre <- phylo.maker(
  sp.list = sp_list_corrected,
  tree = GBOTB.extended.TPL,
  nodes = nodes.info.1.TPL,
  scenarios = "S3"
)

my_tree <- tre$scenario.3
dist_phylo <- as.matrix(cophenetic(my_tree))

# Calculate phylogenetic diversity
layered_rao_results <- layered_phylo_ready %>%
  group_by(siteID) %>%
  group_modify(~ calc_rao_partition(.x, dist_phylo)) %>%
  mutate(Method = "Layered", Scope = "Veg Only (Phylo)")

print(layered_rao_results)

##########################################
# MErge output !
# Join Taxonomic and Phylogenetic results by siteID
taxonomic_metrics <- layered_abundance_diversity_table %>%
  rename(
    Tax_Alpha_Shannon = Alpha_Shannon,
    Tax_Alpha_Simpson = Alpha_Simpson,
    Tax_Gamma_Shannon = Gamma_Shannon,
    Tax_Beta_Bray = Beta_Distance
  )

phylo_metrics <- layered_rao_results %>%
  ungroup() %>%
  mutate(Scope = "Veg Only") %>% 
  select(siteID, Scope, Phylo_Alpha, Phylo_Beta, Phylo_Gamma)

layered_master_table <- taxonomic_metrics %>%
  left_join(phylo_metrics, by = c("siteID", "Scope")) %>%
  select(siteID, Scope, starts_with("Tax_"), starts_with("Phylo_")) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(siteID, desc(Scope))

head(layered_master_table, 10)

write.csv(layered_master_table, "data_out/NEON_Layered_Diversity_Results.csv", row.names = FALSE)
