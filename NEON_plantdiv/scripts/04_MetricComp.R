## 04_Metric_Comp: Comparisons of metrics...

library(tidyverse)
library(corrplot)
library(ggrepel)
library(patchwork)

# load metrics/results
layered_results <- read.csv("data_out/NEON_Layered_Diversity_Results.csv") %>%
  mutate(Method = "Layered")
topdown_results <- read.csv("data_out/NEON_TopDown_Diversity_Results.csv") %>%
  mutate(Method = "Top-Down")

# join em with site info
layered_clean <- layered_results %>%
  rename_with(~paste0("Layered_", .), -siteID)

topdown_clean <- topdown_results %>%
  rename_with(~paste0("TopDown_", .), -siteID)

comparison_final <- inner_join(layered_clean, topdown_clean, by = "siteID")

# 3. Add Biome Info
site_biomes <- data.frame(
  siteID = c("ABBY", "BARR", "BART", "BLAN", "BONA", "CLBJ", "CPER", "DCFS", "DEJU", "DELA", 
             "DSNY", "GRSM", "GUAN", "HARV", "HEAL", "JERC", "JORN", "KONA", "KONZ", "LAJA", 
             "LENO", "MLBS", "MOAB", "NIWO", "NOGP", "OAES", "ONAQ", "ORNL", "OSBS", "PUUM", 
             "RMNP", "SCBI", "SERC", "SJER", "SOAP", "SRER", "STEI", "STER", "TALL", "TEAK", 
             "TOOL", "TREE", "UKFS", "UNDE", "WOOD", "WREF", "YELL"),
  System_Type = c("Temperate rain forest", "Tundra", "Temperate seasonal forest", "Temperate seasonal forest", 
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
                  "Temperate grassland/desert", "Temperate rain forest", "Montane forest")
)

comparison_final <- comparison_final %>%
  left_join(site_biomes, by = "siteID")


# corr matrix...?
cor_data <- comparison_final %>% 
  select(starts_with("Layered_"), starts_with("TopDown_")) %>%
  select(where(is.numeric)) %>% 
  drop_na()

cor_matrix <- cor(cor_data, method = "pearson")

png("figures/Corr_Matrix.png", width = 800, height = 800)
corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black", 
         addCoef.col = "black", number.cex = 0.7, mar = c(0,0,2,0),
         title = "Correlation: Layered Ground Truth vs. Top-Down Canopy")
dev.off()

# plot comparisons for top-down vs layered td and pd
p1 <- ggplot(comparison_final, aes(x = Layered_TD_Beta, y = TopDown_TD_Beta, 
                             label = siteID, color = System_Type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3.5, alpha = 0.7) +
  geom_text_repel(size = 3, fontface = "bold", box.padding = 0.5, max.overlaps = Inf) + 
  scale_color_viridis_d(option = "turbo", name = "System Type") + 
  theme_minimal() +
  labs(title = "Taxonomic Beta Diversity Comparison",
       subtitle = "Top-Down vs layered approacg",
       x = "'Layered' TD Beta (Total Understory + Canopy)",
       y = "'Top-Down' TD Beta ('Visible' Surfaces Only)")

p2 <- ggplot(comparison_final, aes(x = Layered_PD_Beta, y = TopDown_PD_Beta, 
                             label = siteID, color = System_Type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3.5, alpha = 0.7) +
  geom_text_repel(size = 3, fontface = "bold", box.padding = 0.5, max.overlaps = Inf) + 
  scale_color_viridis_d(option = "turbo", name = "System Type") + 
  theme_minimal() +
  labs(title = "Phylogenetic Beta Diversity Comparison",
       subtitle = "Top-Down vs layered approacg",
       x = "'Layered' PD Beta (Total Understory + Canopy)",
       y = "'Top-Down' PD Beta ('Visible' Surfaces Only)")

combined_plot <- p1 + p2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave("figures/Beta_Div_Comp.png", 
       combined_plot, 
       width = 14, 
       height = 7, 
       dpi = 300)
