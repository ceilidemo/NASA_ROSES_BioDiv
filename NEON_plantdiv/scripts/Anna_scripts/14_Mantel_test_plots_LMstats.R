######## Mantel test ########
library(ade4)
library(ggplot2)
library(dichromat)
library(tidyverse)

veg_dist_mats <- readRDS("./R_input/model_input/veg_dist_hellinger01.rds") 
spec_dist_mats <- readRDS("./R_input/model_input/spec_dist_euclid01.rds")
abby <- as.matrix(spec_dist_mats[[1]])

mt <- list()
VD_SD_frame <- list()

for(i in 1:30){
  spec <- spec_dist_mats[[i]]
  veg <- veg_dist_mats[[i]]
  set.seed(1840)
  mt[[i]] <- mantel.randtest(spec,veg,nrepet=999)
  
  #### make data frame 
  VD <- as.data.frame(as.matrix(veg))
  abbrev <- row.names(VD)
  VD <- cbind(abbrev,VD)
  row.names(VD) <- NULL
  
  SD <- as.data.frame(as.matrix(spec))
  abbrev <- row.names(SD)
  SD <- cbind(abbrev,SD)
  row.names(SD) <- NULL
  
  VD <- VD[order(VD$abbrev),order(colnames(VD))]
  SD <- SD[order(SD$abbrev),order(colnames(SD))]
  
  table(VD$abbrev==SD$abbrev)
  
  namSD <- SD[,1]
  MSD <- SD[,-1]
  row.names(MSD) <- namSD
  
  namVD <- VD[,1]
  MVD <- VD[,-1]
  row.names(MVD) <- namVD
  
  ##### combine in one dataframe 
  colnames(MVD)==colnames(MSD)
  row.names(MVD)==row.names(MSD)
  
  nami <- rep(colnames(MVD),each=ncol(MVD))
  VVD <- as.vector(as.matrix(MVD))
  VSD <- as.vector(as.matrix(MSD))
  
  dat <- as.data.frame(cbind(VVD, VSD))
  dat <- cbind(nami, dat)
  
  as.character(namVD)==namSD
  plotID <- rep(namSD, times=ncol(MVD))
  VD_SD_frame[[i]] <- cbind(dat, plotID)
}

saveRDS(mt, "./R_output/Mantel/res_mantel_sites.rds")
saveRDS(VD_SD_frame, "./R_output/Mantel/VD_SD_frames.rds")

mt_df <- matrix(nrow=30,ncol=6, 
                dimnames = list(NULL,c("plotID","R2", "pval","Std.Obs", 
                                       "Expectation", "Variance")))
for(i in 1:30){
  mt_df[i,1]<- substr(VD_SD_frame[[i]]$plotID[1],1,4)
  mt_df[i,2]<- mt[[i]]$obs
  mt_df[i,3]<- mt[[i]]$pvalue
  mt_df[i,4:6]<- as.numeric(mt[[i]]$expvar)
}

mt_df <- as.data.frame(mt_df)
write.csv(mt_df, "./R_output/Mantel/res_mantel_compiled.csv",row.names = F)

#######################
#### Plots per Site
VD_SD_frame <- readRDS("./R_output/Mantel/VD_SD_frames.rds") 
# coliix <- dichromat(colours = rainbow(n=ncol(MFD)),type="deutan")
# coliix <- dichromat(colours = rainbow(n=ncol(MFD)),type="protan")
# plot(1:length(coliix), 1:length(coliix), col=coliix, pch=19, cex=3, xlab="", ylab="")

dat <- do.call(rbind,VD_SD_frame)
dat <- dat %>% mutate(siteID=substr(nami,1,4)) %>% group_by(siteID)

coliix <- dichromat(colours = rainbow(n=max(table(dat$plotID))),type="tritan")
nlin <- ceiling(max(table(dat$plotID))/3)
lini <- c(rep(1,times=nlin), rep(5, times=nlin),rep(4, times=nlin))

sel <- dat %>% filter(siteID=="LAJA")

#### one Site
ggplot(data=sel, aes(x=VVD, y=VSD)) +
  scale_y_continuous(limit= c(0,1),breaks=c(0,0.5,1))+
  scale_x_continuous(breaks=c(0,0.5,1))+
  geom_smooth(aes(color=plotID, linetype=plotID),method="lm",lwd=0.6,
              show.legend = T, se = F, fullrange=F)+
  geom_point(color="grey30", pch=1, size=2)+
  scale_linetype_manual(values=c(lini))+
  scale_color_manual(values=coliix)+
  # scale_color_viridis_d(option="viridis")+
  labs(x="Vegetation distance",y="Spectral distance")+
  # geom_smooth(method = "lm", formula = y~-1+x+I(x^2), col="black", lwd=0.7,
  #             lty=2, alpha=0.4, fullrange=T)
  geom_smooth(method = "lm", color="black",
              formula = y~x-1,  lwd=1.2,
              lty=1, alpha=0.4, fullrange=T, level=0.95)+
  # ggtitle(label = substr(namFD[1],1,4))+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major =  element_line(colour = "grey",size=0.2),
        # strip.text= element_text(size=12),
        # legend.text = element_text(size = 9),
        # legend.title = element_text(size = 12),
        # legend.box.spacing = unit(1.8,"line"),
        # legend.key.width = unit(2.2,"line"),
        axis.title.x = element_text(margin = margin(t = 15), size=18),
        axis.title.y = element_text(margin = margin(r = 15), size=18),
        plot.title =element_text(size=14, face="bold", hjust = -0.15),
        legend.key = element_rect(size = 0.1, fill = "white"),
        legend.key.size = unit(1.2,"line"),
        # text=element_text(size=12),
        axis.ticks = element_line(colour = "grey", size = 0.1)) 
ggsave("./R_output/figures/VD_SD_ABBY.pdf") 

#################
##### STRIP PLOT
ggplot(data=dat, aes(x=VVD, y=VSD)) +
  scale_y_continuous(limit= c(0,1),breaks=c(0,0.5,1))+
  scale_x_continuous(breaks=c(0,0.5,1))+
  geom_point(color="grey60", pch=1, size=1.2)+
  geom_smooth(aes(group = plotID),method="lm",color="grey10" ,lwd=0.4,
              show.legend = F, se = F, fullrange=F)+
  facet_wrap(~siteID,scales = "free")+
  # scale_linetype_manual(values=c(lini))
  # scale_color_manual(values=coliix)+
  # scale_color_viridis_d(option="viridis")+
  labs(x="Taxonomic distance",y="Spectral distance")+
  # geom_smooth(method = "lm", formula = y~-1+x+I(x^2), col="black", lwd=0.7,
  #             lty=2, alpha=0.4, fullrange=T)
  geom_smooth(method = "lm", formula = y~x-1, col="red3", lwd=1,
              lty=1, alpha=0.4, fullrange=T, level=0.95)+
  # ggtitle(label = substr(namFD[1],1,4))+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major =  element_line(colour = "grey",size=0.2),
        strip.text= element_text(size=14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 12),
        legend.box.spacing = unit(1.8,"line"),
        legend.key.width = unit(2.2,"line"),
        axis.title.x = element_text(margin = margin(t = 15), size=18),
        axis.title.y = element_text(margin = margin(r = 15), size=18),
        plot.title =element_text(size=14, face="bold", hjust = -0.15),
        legend.key = element_rect(size = 0.1, fill = "white"),
        legend.key.size = unit(1.2,"line"),
        text=element_text(size=12),
        axis.ticks = element_line(colour = "grey", size = 0.1)) 
ggsave("./R_output/figures/SD_TD/SD_TD_strips_ylim.pdf",width = 12, height = 8)  
  
####### Sites combined  #####

coliix <- dichromat(colours = rainbow(n=length(unique(dat$siteID))),type="tritan")
nlin <- ceiling(length(unique(dat$siteID))/3)
lini <- c(rep(1,times=nlin), rep(5, times=nlin),rep(4, times=nlin))

### color by latitute 
info <- read.csv("./pub/NEON_sites_overview.csv")
names(info)[1] <- "siteID"
ordi <- order(info$Latitude)
s <- unique(dat$siteID)
s_ord <- s[rev(ordi)]

dat_ord <- dat %>% left_join(info[,c(1,9)])%>% 
  mutate(siteID=factor(siteID, levels=s_ord))

levels(dat_ord$siteID)

ggplot(data=dat_ord, aes(x=VVD, y=VSD)) +
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  scale_x_continuous(breaks=c(0,0.5,1))+
  # stat_density_2d(aes(fill = ..level..), geom="polygon", show.legend = F, alpha=0.2)+
  geom_hex(bins=30,show.legend = F)+
  scale_fill_gradientn(name = "Number of observations", 
                       colours=rev(c("black","snow")),
                       n.breaks = 4,values = c(0,0.1,1))+
  # theme(legend.position = "bottom", 
  #       legend.direction = "horizontal",
  #       legend.key.width = unit(1,"cm"))
  # scale_fill_gradientn(name = "Number\nof observations\n", 
  #                      colours=c("skyblue1","lightgoldenrod1", "red1"),
  #                      n.breaks = 6)+
  # scale_fill_gradientn(name = "Number\nof observations\n", 
  #                      colours=rev(c(1,2,"gold","ivory")),n.breaks = 6)+
  # geom_point(color="grey30", pch=1, size=2, alpha=0.1)+
  geom_smooth(aes(group = siteID, color=siteID),method="lm",
              lwd=1,show.legend = T, se = F, fullrange=F)+
  # scale_linetype_manual(values=c(lini))+
  # scale_color_manual(values=coliix)+
  scale_color_viridis_d(option="magma")+
  labs(x="\nTaxonomic distance",y="Spectral distance\n")+
  # geom_smooth(method = "lm", formula = y~-1+x+I(x^2), col="black", lwd=0.7,
  #             lty=2, alpha=0.4, fullrange=T)
  geom_smooth(method = "lm", formula = y~x-1, lwd=2,
              lty=6, alpha=0.4, fullrange=T, level=0.95, color="gold")+
  # ggtitle(label = substr(namFD[1],1,4))+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major =  element_line(colour = "grey",size=0.2),
        strip.text= element_text(size=12),
        legend.text = element_text(size = 9),
        legend.title = element_blank(),
        legend.position = "right",
        legend.box.spacing = unit(1.8,"line"),
        legend.key.width = unit(2.2,"line"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        plot.title =element_text(size=14, face="bold", hjust = -0.15),
        # legend.key = element_rect(size = 0.1, fill = "white"),
        legend.key.size = unit(1.2,"line"),
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        # text=element_text(size=12),
        axis.ticks = element_line(colour = "grey", size = 0.1)) 
ggsave("./R_output/figures/SD_TD/SD_TD_noobs_bw.pdf",width = 8, height = 5) 

#################################################################
####### Correlations go through origin, when tax dist is zero, spec dist is zero

sel <- dat %>% filter(siteID=="ABBY")
seli <- dat %>% filter(nami=="ABBY_004")

#### one Site
ggplot(data=sel, aes(x=VVD, y=VSD)) +
  scale_y_continuous(limit= c(0,1),breaks=c(0,0.5,1))+
  scale_x_continuous(breaks=c(0,0.5,1))+
  geom_smooth(aes(color=plotID),method="lm", lwd=0.6,
              show.legend = T, se = F, fullrange=F,formula=y~x-1)+
  geom_point(color="grey30", pch=1, size=2)+
  scale_color_viridis_d(option="viridis")+
  labs(x="Vegetation distance",y="Spectral distance")+
  geom_smooth(method = "lm", color="black",
              formula = y~x-1,  lwd=1.2,
              lty=1, alpha=0.4, fullrange=T, level=0.95)+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major =  element_line(colour = "grey",size=0.2),
        axis.title.x = element_text(margin = margin(t = 15), size=18),
        axis.title.y = element_text(margin = margin(r = 15), size=18),
        plot.title =element_text(size=14, face="bold", hjust = -0.15),
        legend.key = element_rect(size = 0.1, fill = "white"),
        legend.key.size = unit(1.2,"line"),
        # text=element_text(size=12),
        axis.ticks = element_line(colour = "grey", size = 0.1)) 
ggsave("./R_output/figures/VD_SD_ABBY_origin.pdf")

#################
##### STRIP PLOT
ggplot(data=dat, aes(x=VVD, y=VSD)) +
  scale_y_continuous(limit= c(0,1),breaks=c(0,0.5,1))+
  scale_x_continuous(breaks=c(0,0.5,1))+
  geom_point(color="grey60", pch=1, size=1.2)+
  geom_smooth(aes(group = plotID),method="lm",color="grey10" ,lwd=0.4,
              show.legend = F, se = F, fullrange=F, formula=y~x-1)+
  facet_wrap(~siteID,scales = "free")+
  labs(x="Taxonomic distance",y="Spectral distance")+
  # geom_smooth(method = "lm", formula = y~-1+x+I(x^2), col="black", lwd=0.7,
  #             lty=2, alpha=0.4, fullrange=T)
  geom_smooth(method = "lm", formula = y~x-1, col="red3", lwd=1,
              lty=1, alpha=0.4, fullrange=T, level=0.95)+
  # ggtitle(label = substr(namFD[1],1,4))+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major =  element_line(colour = "grey",size=0.2),
        strip.text= element_text(size=14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 12),
        legend.box.spacing = unit(1.8,"line"),
        legend.key.width = unit(2.2,"line"),
        axis.title.x = element_text(margin = margin(t = 15), size=18),
        axis.title.y = element_text(margin = margin(r = 15), size=18),
        plot.title =element_text(size=14, face="bold", hjust = -0.15),
        legend.key = element_rect(size = 0.1, fill = "white"),
        legend.key.size = unit(1.2,"line"),
        text=element_text(size=12),
        axis.ticks = element_line(colour = "grey", size = 0.1)) 
ggsave("./R_output/figures/SD_TD/SD_TD_strips_ylim_origin.pdf",width = 12, height = 8)

### Combi plot all sites
ggplot(data=dat_ord, aes(x=VVD, y=VSD)) +
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  scale_x_continuous(breaks=c(0,0.5,1))+
  geom_hex(bins=30,show.legend = F)+
  scale_fill_gradientn(name = "Number of observations", 
                       colours=rev(c("black","snow")),
                       n.breaks = 4,values = c(0,0.1,1))+
  geom_smooth(aes(group = siteID, color=siteID),method="lm",
              lwd=1,show.legend = T, se = F, fullrange=F, formula=y~x-1)+
  scale_color_viridis_d(option="magma")+
  labs(x="\nTaxonomic distance",y="Spectral distance\n")+
  geom_smooth(method = "lm", formula = y~x-1, lwd=2,
              lty=6, alpha=0.4, fullrange=T, level=0.95, color="gold")+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major =  element_line(colour = "grey",size=0.2),
        strip.text= element_text(size=12),
        legend.text = element_text(size = 9),
        legend.title = element_blank(),
        legend.position = "right",
        legend.box.spacing = unit(1.8,"line"),
        legend.key.width = unit(2.2,"line"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        plot.title =element_text(size=14, face="bold", hjust = -0.15),
        # legend.key = element_rect(size = 0.1, fill = "white"),
        legend.key.size = unit(1.2,"line"),
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        # text=element_text(size=12),
        axis.ticks = element_line(colour = "grey", size = 0.1)) 
# ggsave("./R_output/figures/SD_TD/SD_TD_noobs_bw_origin.pdf",width = 8, height = 5)


#############################################
#### Model fit overall without intercept #### 
VD_SD_frame <- readRDS("./R_output/Mantel/VD_SD_frames.rds") 

dat <- do.call(rbind,VD_SD_frame)
dat <- dat %>% mutate(siteID=substr(nami,1,4)) %>% group_by(siteID)


dat2 <- dat[!(dat$nami==dat$plotID),]

R2noint<- function(m,y){
  as.numeric(1 - crossprod(residuals(m))/crossprod(y - mean(y)))
}

### Zeros don't matter for overall fit
mod1 <- lm(VSD ~ VVD-1, data=dat2)
R2noint(mod1, dat$VSD)
summary(mod1)

mod <- lm(VSD ~ VVD-1, data=dat)
R2noint(mod, dat$VSD)
summary(mod)

#### Tidy fit site models
#### model stats from model w/o zero

### TRY with ABBY
# selx <- sel[!(sel$nami==sel$plotID),]
# 
# abby_mod <- lm(VSD ~ VVD-1, data = selx)
# glance(abby_mod)
# summary(abby_mod)

by_site <- dat2 %>% 
  group_by(siteID) %>% 
  nest()

site_mod <- function(df) {
  lm(VSD ~ VVD-1, data = df)
}

mods <- map(by_site$data, site_mod)

tidy(mods[[1]])
glance(mods[[1]])

### Compile stats
sumi <- by_site %>% 
  add_column(tidy =map(mods, broom::tidy))%>% 
  unnest(tidy)%>%
  add_column(glance =map(mods, broom::glance))%>% 
  unnest(glance,.names_repair = "universal")

modi <- dat %>% group_by(siteID) %>% do(fit = lm(VSD ~ VVD-1, .))

R2 <- numeric()
for (i in 1:30){
  R2[i] <- as.numeric(1- crossprod(unlist(modi[i,2]$fit[[1]][2])) / 
                        crossprod(modi[i,2]$fit[[1]][12]$model[,1]-
                                    mean(modi[i,2]$fit[[1]][12]$model[,1])))
}

sumi$r2_noint <-R2
sumi <- sumi[,-2]

write.csv(sumi, "./R_output/Mantel/VD_SD_fit_persite_NOzero.csv", row.names = F)

#### END ########

