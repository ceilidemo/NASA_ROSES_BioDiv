# playing around

emit_data <- read.csv("EMIT_specdiv/data_out/results/EMIT_SpecDiv_summary.csv") %>% 
  mutate(siteID = site) %>% select(siteID, gamma_ss, gamma_sdiv)

final_master <- comparison_final %>%
  left_join(emit_data, by = "siteID")


corr_data <- final_master %>%
  ungroup() %>%                     # CRITICAL: Remove the hidden siteID group
  select(gamma_ss, gamma_sdiv, 
         Layered_TD_Richness, Layered_TD_Beta, Layered_PD_Gamma,
         TopDown_TD_Richness, TopDown_TD_Beta, TopDown_PD_Gamma) %>%
  drop_na() %>% 
  select(where(is.numeric))

M <- cor(corr_data, method = "pearson")

library(corrplot)
png("Modeling/figures/EMITvPlants_corr.png", width = 800, height = 800, res = 120)
corr <- corrplot(M, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", 
         title = "EMIT & Plants Corr",
         mar=c(0,0,1,0))
dev.off()

