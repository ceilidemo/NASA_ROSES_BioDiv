## 04_Metric_Comp: Comparisons of metrics...
library(tidyverse)
library(corrplot)
library(ggplot2)
library(ggrepel)

#Load in the results 
layered_df <- read.csv("data_out/NEON_Layered_Diversity_Results.csv")
topdown_df <- read.csv("data_out/NEON_TopDown_Diversity_Results.csv")

# Clean up and join
layered_clean <- layered_master_table %>%
  ungroup() %>%
  select(-Scope) %>%
  rename_with(~paste0("Layered_", .), -siteID)

topdown_clean <- final_topdown_diversity %>%
  ungroup() %>%
  rename_with(~paste0("TopDown_", .), -siteID)

comparison_final <- inner_join(layered_clean, topdown_clean, by = "siteID")


# Corr matrix 
cor_data <- comparison_final %>% 
  select(starts_with("Layered_"), starts_with("TopDown_")) %>%
  drop_na()

cor_matrix <- cor(cor_data, method = "pearson")

corrplot(cor_matrix, 
         method = "color", 
         type = "upper", 
         tl.col = "black", 
         addCoef.col = "black", 
         number.cex = 0.7,
         title = "Correlation: Layered Ground Truth vs. Top-Down Canopy",
         mar = c(0,0,2,0))

## Add in the biome info from the 2022 nature paper
site_biomes <- data.frame(
  siteID = c(
    "ABBY", "BARR", "BART", "BLAN", "BONA", "CLBJ", "CPER", "DCFS", "DEJU", "DELA", 
    "DSNY", "GRSM", "GUAN", "HARV", "HEAL", "JERC", "JORN", "KONA", "KONZ", "LAJA", 
    "LENO", "MLBS", "MOAB", "NIWO", "NOGP", "OAES", "ONAQ", "ORNL", "OSBS", "PUUM", 
    "RMNP", "SCBI", "SERC", "SJER", "SOAP", "SRER", "STEI", "STER", "TALL", "TEAK", 
    "TOOL", "TREE", "UKFS", "UNDE", "WOOD", "WREF", "YELL"
  ),
  System_Type = c(
    "Temperate rain forest", "Tundra", "Temperate seasonal forest", "Temperate seasonal forest", 
    "Boreal forest", "Temperate grassland/desert", "Temperate grassland/desert", "Boreal forest", 
    "Boreal forest", "Temperate seasonal forest", "Subtropical desert", "Temperate seasonal forest", 
    "Tropical rain forest", "Temperate seasonal forest", "Boreal forest", "Temperate seasonal forest", 
    "Subtropical desert", "Temperate grassland/desert", "Temperate grassland/desert", "Tropical seasonal forest/savanna", 
    "Tropical seasonal forest/savanna", "Temperate seasonal forest", "Subtropical desert", "Alpine", 
    "Temperate grassland/desert", "Temperate grassland/desert", "Subtropical desert", "Temperate seasonal forest", 
    "Subtropical desert", "Tropical rain forest", "Montane forest", "Temperate seasonal forest", 
    "Temperate seasonal forest", "Woodland/shrubland", "Woodland/shrubland", "Woodland/shrubland", 
    "Temperate grassland/desert", "Temperate grassland/desert", "Temperate grassland/desert", "Montane forest", 
    "Tundra", "Temperate seasonal forest", "Temperate seasonal forest", "Temperate seasonal forest", 
    "Temperate grassland/desert", "Temperate rain forest", "Montane forest"
  )
)

comparison_final <- comparison_final %>%
  left_join(site_biomes, by = "siteID") # bring it on in

# TD beta plot w/site
ggplot(comparison_final, aes(x = Layered_TD_Beta, y = TopDown_TD_Beta, 
                              label = siteID, color = System_Type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3.5, alpha = 0.7) +
  geom_text_repel(size = 3, 
                  fontface = "bold",
                  box.padding = 0.5, 
                  point.padding = 0.3,
                  force = 10,
                  max.overlaps = Inf,
                  show.legend = FALSE) + 
  # Using 'turbo' because it handles more than 9 categories easily
  scale_color_viridis_d(option = "turbo", name = "System Type") + 
  theme_minimal() +
  labs(title = "Taxonomic Beta Diversity",
       subtitle = "Colored site by biome type",
       x = "'Layered' Approach TD Beta",
       y = "Scaled Top-Down Approach TD Beta") +
  theme(legend.position = "right")

# PD
ggplot(comparison_final, aes(x = Layered_PD_Beta, y = TopDown_PD_Beta, 
                              label = siteID, color = System_Type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3.5, alpha = 0.7) +
  geom_text_repel(size = 3, 
                  fontface = "bold",
                  box.padding = 0.5, 
                  point.padding = 0.3,
                  force = 10,
                  max.overlaps = Inf,
                  show.legend = FALSE) + 
  scale_color_viridis_d(option = "turbo", name = "System Type") + 
  theme_minimal() +
  labs(title = "Phylogenetic Beta Diversity",
       subtitle = "Colored site by biome",
       x = "Layered' Approach PD Beta",
       y = "Scaled Top-Down Approach PD Beta") +
  theme(legend.position = "right")

