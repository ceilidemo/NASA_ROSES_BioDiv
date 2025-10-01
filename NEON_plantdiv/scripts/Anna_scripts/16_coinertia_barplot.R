library(ggplot2)
# remotes::install_github("coolbutuseless/ggpattern")
# library(ggpattern)
library(tidyverse)

dat <- read.csv("./R_output/ordination_manova/coinertia_MonteCarlo.csv")
# test <-readRDS("./R_output/ordination_manova/coinertia_avspectra/ABBY_COIA.rds")

info <- read.csv("./R_input/sites_biom_color.csv")
names(info)[2] <- "siteID"

dati <- dat %>% 
  mutate(pclass=case_when(p_val <= 0.001 ~ 3, 
                          p_val <= 0.01 ~ 2,
                          p_val <= 0.05 ~ 4,
                          p_val >0.05~ 1))%>% 
  # mutate(patt=case_when(p_val <= 0.001 ~ "crosshatch", 
  #                         p_val <= 0.01 ~ "stripe",
  #                         p_val <= 0.05 ~ "point",
  #                         p_val >0.05~ "none"))%>% 
  mutate(densi=case_when(p_val <= 0.001 ~ 0.1, 
                         p_val <= 0.01 ~ 0.1,
                         p_val <= 0.05 ~ 0.1,
                         p_val >0.05~ 0))%>% 
  mutate(alph=case_when(p_val <= 0.001 ~ 1, 
                         p_val <= 0.01 ~ 1,
                         p_val <= 0.05 ~ 1,
                         p_val >0.05~ 0.3))%>% 
  mutate(stars=case_when(p_val <= 0.001 ~ "***", 
                        p_val <= 0.01 ~ "**",
                        p_val <= 0.05 ~ "*",
                        p_val >0.05~ ""))%>% 
  left_join(info)%>% 
  mutate(siteID=as.factor(siteID))%>%
  mutate(siteID=reorder(siteID, RV_obs))
  # mutate(siteID = reorder(siteID, Latitude))


##### plot with alpha
ggplot(data = dati, aes(y = RV_obs, x = siteID, width=0.85)) +
  geom_col(fill=dati$Biome_col, alpha=dati$alph)+
  coord_flip()+ theme_minimal()+
  ylim(0,0.8)+
  xlab("")+ylab("Covariance")+
  geom_text(aes(label = stars), hjust = -0.2, vjust=0.8, cex=7)+
  theme(axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        # legend.spacing.y = unit(1.5,"line"),
        # legend.key.height = unit(2,"line"),
        # legend.text = element_text(size = 14),
        # legend.title = element_text(size = 16),
        # legend.box.spacing = unit(3,"line"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))
ggsave("./R_output/figures/coinertia_stars.pdf")



#### pattern plot
dens <- dati$densi
patti <- as.character(dati$patt)

table(dati$patt, dati$pclass)

ggplot(data = dati, aes(y = RV_obs, x = siteID, width=0.85)) +
  # geom_col() +
  geom_col_pattern(aes(pattern=pclass),
                   colour  = 'black', fill=coli,
                   pattern_density=dens, show.legend = F)+
  coord_flip()+ theme_minimal()+
  xlab("")+ylab("Covariance")+
  theme(axis.title.x = element_text(size=26),
        axis.title.y = element_text(size=26),
        # legend.spacing.y = unit(1.5,"line"),
        # legend.key.height = unit(2,"line"),
        # legend.text = element_text(size = 14),
        # legend.title = element_text(size = 16),
        # legend.box.spacing = unit(3,"line"),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18))
ggsave("./R_output/figures/Coinertia.pdf")
