library(picante)
library(beepr)
library(tidyverse)

#### Phylogenetic diveristy metrics 
phy_list <- readRDS("./R_output/phylo/allplants_scenario2_tree_100.rds")
df_phylo <- read.csv("./R_output/phylo/veg_df_phylo.csv", row.names = 1, check.names = F)

tre <- phy_list[[1]]
names(df_phylo)[!(names(df_phylo) %in% tre$run.1$tip.label)] 

# phy <- phy_list$scenario.2

xx <- list()
for(i in 1:100){
  pphy <- drop.tip(phy[[i]], setdiff(phy[[i]]$tip.label, colnames(df_phylo)))
  psvobs <- psv(df_phylo, pphy,compute.var = F) 
  pseobs <- pse(df_phylo, pphy)
  psrobs <- psr(df_phylo, pphy,compute.var = F)
  
  faith1 <- pd(df_phylo, pphy, include.root=FALSE)
  faith1[is.na(faith1)] <- 0
  
  #combine data
  psv <- psvobs %>% rownames_to_column("plotID") %>% select(-SR)
  pse <- pseobs %>% rownames_to_column("plotID") %>% select(-SR)
  psr <- psrobs %>% rownames_to_column("plotID") %>% select(-SR)
  pd <- faith1 %>% rownames_to_column("plotID") %>% select(-SR)
  
  xx[[i]] <- psv %>% left_join(pse)%>% left_join(psr)%>% left_join(pd)
}
beep(2)
saveRDS(xx, "./R_output/phylo/phylo_metrics100.rds")


