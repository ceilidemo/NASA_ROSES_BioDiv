#######
library(tidyverse)
library(ade4)
library(spectrolab)
library(stringr)
library(dave)
library(ggpubr)

### average spectra ###
ss <- readRDS("~/Documents/CABO/NEON/R_input/model_input/average_spectra/av_VNspectra.rds")

spec_divi <- read.csv("./R_input/model_input/NEON_diversity.csv")
spec_div <- spec_divi %>% replace_na(list(todo="NA")) %>%
  mutate(siteID=substr(plotID,1,4), .before=1)%>%
  filter(!(grepl("ut", todo)))%>%
  filter(n_pix_NDVI>=50)

### inventories
df <- read.csv("./R_output/inv_compiled/veg_cover_compiled_token.csv")

### Veg dataframe including NA species, moss, lichen, etc
div <- df %>% select(-c(2:8)) %>%
  select(-c(wood_NA, water_NA,standingDead_NA, litter_NA, rock_NA, soil_NA,scat_NA)) %>%
  filter(plotID %in% spec_div$plotID)

sites <- names(ss)
div_sites <- unique(substr(div$plotID,1,4))
good_sites <- unique(substr(spec_div$plotID,1,4))

sites[!(sites %in% good_sites)]
sites <- sites[sites %in% good_sites]

ss[(which(!(names(ss) %in% sites),arr.ind = T))] <- NULL

test <- ss[[1]]
test_spec <- as_spectra(test,name_idx=2, meta_idxs = c(1,3:5))
plot(test_spec)

table(sites==names(ss))

unique(spec_div$todo)

### index PC axes spec div
# ll <- list.files("//Volumes/Backup Plus/Anna/NEON/R_output/extract_specdiv_plants_noshade/",
#                  full.names = T, recursive = T)
# l <- list.files("//Volumes/Backup Plus/Anna/NEON/R_output/extract_specdiv_plants_noshade/",
#                       recursive = T)

### missing tree diameters 
# dd <- list.files("//Volumes/Backup Plus/Anna/NEON/R_output/figures/woody_diameter_missing/",
#                 pattern = "missing")
# ddd <- substr(dd,1,8)

coia_randu <- list()
eigen <- list()
#### start with ABBY
for (i in 1:length(sites)){
  spec <- ss[[i]]
  spec <- spec[order(spec$plotID),]
  spec <- spec %>% select_if(~all(!is.na(.)))
  
  divi <- div %>% filter(grepl(sites[i],plotID))
  divx <- divi[,c(1,which(as.numeric(colSums(divi[,-1]))> 0)+1)]
  divx <- divx[order(divx$plotID),]
  
  ### match dfs
  nam_v <- divx$plotID
  nam_s <- spec$plotID
  
  specx <- spec[spec$plotID %in% nam_v,]
  divxx <- divx[divx$plotID %in% nam_s,]
  
  table(specx$plotID==divxx$plotID)
  nami <- divxx$plotID
  
  ### site variables 
  site_vars <- spec_div %>% filter(plotID %in% nami)
  
  #################
  ### PCs if needed 
  # PC_dat <- read.csv(paste0("./R_input/model_input/PCs_sites/", sites[i],"_PCs.csv"))
 
  ### Hellinger transformation adjusts species vectors to equal sum of squares =1
  div_hell <- decostand(divxx[, -c(1)], method="hellinger") 
  div_hell <- div_hell %>% add_column(plotID = nami,.before = 1) 
  
  ### distance matrices 
  div_mat <- divxx[, -c(1)]
  row.names(div_mat) <- divxx[,1] 
  
  div_mat_hell <- div_hell[,-1] 
  row.names(div_mat_hell) <- div_hell[,1] 
  
  spec_mat <- specx[,-c(1:5)] 
  row.names(spec_mat) <- specx[,2] 
  
  row.names(div_mat_hell)==row.names(spec_mat)
  

  ### PCO vegetation ### 
  veg_dist_hell <- vegdist(div_mat_hell,method = "eucl") #### Hellinger Dist
  veg_dist_bray <-  vegdist(div_mat,method = "bray") ### Bray-Curtis
  
  veg_pco_hell <- pco(dis = veg_dist_hell)
  veg_pco_bray <- pco(dis = veg_dist_bray)

  # plot(veg_pco_hell, pch=16)
  # try(surf(veg_pco_hell, site_vars$SDiv_alpha, col=1, thinplate = F))

  #### PCO Spectra
  # spec_dist_manhattan <- vegdist(spec_mat, method = "manhatt")
  # spec_dist_mahala <- vegdist(spec_mat, method = "mahala")###
  spec_dist_euclid <- vegdist(spec_mat, method = "euclid")###
  # spec_pco_manhattan <- pco(dis = spec_dist_manhattan)
  # spec_pco_mahala <- pco(dis = spec_dist_mahala)
  # spec_pco_euclid <- pco(dis = spec_dist_euclid)
  
  # plot(spec_pco_manhattan, pch=16)
  # try(surf(spec_pco_manhattan, site_vars$SDiv_alpha, col=1, thinplate = F))

  # ### Procrustes analysis 1 
  # set.seed(1984)
  # (pro <- protest(veg_pco_hell, spec_pco_euclid,permutations = 999, scale=T,symetric=T))
  # pdf(paste0("./R_output/ordination_manova/procrustes_avspectra/", sites[i], 
  #            "_pro_errors_specEuclid_vegHell.pdf"),width = 4,
  #     height = 4)
  # plot(pro)
  # dev.off()
  # 
  # sink(paste0("./R_output/ordination_manova/procrustes_avspectra/", sites[i],
  #             "_specEuclid_vegHell.txt"))
  # print(pro)
  # sink()
  # 
  ### Procrustes analysis 2
  # set.seed(1984)
  # (pro <- protest(veg_pco_hell, spec_pco_mahala,permutations = 999, scale=T,symetric=T))
  # pdf(paste0("./R_output/ordination_manova/procrustes_avspectra/", sites[i], 
  #            "_pro_errors_specEuclid_vegHell.pdf"),width = 4,
  #     height = 4)
  # plot(pro)
  # dev.off()
  # 
  # sink(paste0("./R_output/ordination_manova/procrustes_avspectra/", sites[i],
  #             "_specEuclid_vegHell.txt"))
  # print(pro)
  # sink()
# }

#  (pro <- protest(veg_pco_hell, spec_pco_euclid,permutations = 999, scale=T,symetric=T))
# xx <- plot(pro, kind=1, choices=c(1,2), to.target = TRUE,  type = "p")

# fitted(pro)

# points(pro, display = c("target", "rotated"),
#        choices = c(1,2), truemean =F, col=c(2),pch=16)

############################
#### Coinertia analysis (Xavier's script)

### trafo dist into quasi-euclidean
# veg_dist_quasi <- quasieuclid(veg_dist_hell)
# spec_dist_quasi <- quasieuclid(spec_dist)

# dudi.veg <- dudi.pco(veg_dist_hell, full=TRUE, scannf = FALSE)
### this is the same as PCA Hellinger transformed data and Euclid dist. 
dudi.vegPCA <- dudi.pca(div_mat_hell,scale = F, scannf = FALSE, nf=4) ###
# scatter(dudi.vegPCA)
# dudi.vegPCA
# summary(dudi.vegPCA)


# dudi.spec <- dudi.pco(spec_dist_euclid, full=TRUE,scannf = FALSE)
## the same as PCA using Euclidean distances, only issue if other dist measure required
dudi.specPCA <- dudi.pca(spec_mat, scale=F, scannf = FALSE, nf=4)
# scatter(dudi.specPCA)

# relative variation of eigenvalues
dudi.specPCA$eig/sum(dudi.specPCA$eig)
dudi.vegPCA$eig/sum(dudi.vegPCA$eig)

# equal row weigths in 2 analysis ? (Answer should be TRUE)
# all.equal(dudi.vegPCA$lw,dudi.spec$lw)

# Perform Co-inertia analysis
# coia_veg.spec <-coinertia(dudi.veg, dudi.spec, scannf = FALSE, nf=3)
coia_veg.spec <-coinertia(dudi.vegPCA, dudi.specPCA, scannf = FALSE, nf=4)
  
# coia_veg.spec$tab
# summary(coia_veg.spec)
# plot(coia_veg.spec) ### PCA shows labs for canonical weights 


# relative variation of the first eigenvalues
Eig<-function(x){
  out<-matrix(0,nrow=5,ncol=1)
  for(i in 1:5){
    eigi<-x$eig[i]/sum(x$eig)
    out[i,]<-eigi
  }
  rownames(out)<-c('eigenV1','eigenV2','eigenV3','eigenV4','eigenV5')
  colnames(out)<-c('relative_variance')
  return(out)
}

# apply function
eigen[[i]]<-as.data.frame(Eig(coia_veg.spec))

# graphical representation of the relative variance of the first eigenvalues
# V.propre <-ggplot(data= eigen, aes(x=rownames(eigen), y=relative_variance))+
#   geom_bar(stat = 'identity', width = 0.7)+
#   xlab("")+
#   ylab("Proportion of variance explained")+
#   theme_bw()+
#   scale_y_continuous(expand = c(0,0), limits = c(0,1))

# V.propre

# view results
summary(coia_veg.spec)

coia_randu[[i]] <- RV.rtest(div_mat_hell, spec_mat,nrepet = 999)

# permutation test
sink(paste0("./R_output/ordination_manova/coinertia_avspectra/", sites[i],
            "_COIA_specEuclid_vegHell.txt"))
print(RV.rtest(div_mat_hell, spec_mat, nrepet=999))
print("####################################")
print(summary(coia_veg.spec))
sink()


# plot(randtest(coia_veg.spec,nrepet = 999))
#RV= 0.34 P= 0.001

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Extract result from coInertia
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# extract species procruste residuals
spec.coord<-data.frame(coia_veg.spec$lX)
veg.coord<-data.frame(coia_veg.spec$lY)

# merge both data frame and calculate de distance between the two ordinations
CoIax <-cbind(spec.coord,veg.coord)
CoIax$dist<-((CoIax$AxcY1-CoIax$AxcX1)^2 + (CoIax$AxcY2-CoIax$AxcX2)^2)^(0.5)

CoIax$plotID <-rownames(CoIax)

# add SDiv and site vars
CoIa <- CoIax %>%
  merge(site_vars[,c("plotID","SDiv_alpha","shannon", "meanLAI","vegtype_3",
                     "tree_cov", "herb_cov")], by = "plotID", all.x = T)%>%
  drop_na()

coli <- c("darkgreen","goldenrod2","brown")

# graphical representation of procruste residuals
p1 <-ggplot(data=CoIa, aes(x=reorder(plotID, -dist), y= dist))+
  geom_point(aes(colour=SDiv_alpha),size=4)+
  scale_color_viridis_c(option = "viridis",name="SDiv alpha")+
  # mycolors_Scale+
  xlab('')+
  ylab('Residuals')+
  guides(size = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1))+
  scale_shape_discrete(name="")

p2 <- ggplot(data=CoIa, aes(x=reorder(plotID, -dist), y= dist))+
  geom_point(aes(colour=meanLAI),size=4)+
  scale_color_viridis_c(name="LAI")+
  # mycolors_Scale+
  xlab('')+
  ylab('Residuals')+
  guides(size = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1))+
  scale_shape_discrete(name="")

p6 <- ggplot(data=CoIa, aes(x=reorder(plotID, -dist), y= dist))+
  geom_point(aes(colour=vegtype_3),size=4)+
  scale_color_manual(values=coli,name="Vegetation type")+
  xlab('')+
  ylab('Residuals')+
  guides(size = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1))+
  scale_shape_discrete(name="")


p3 <- ggplot(data=CoIa, aes(x=reorder(plotID, -dist), y= dist))+
  geom_point(aes(colour=tree_cov),size=4)+
  scale_color_viridis_c(option = "viridis",name="Tree cover %",na.value = "grey60")+
  xlab('')+
  ylab('Residuals')+
  guides(size = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1))+
  scale_shape_discrete(name="")

p5 <- ggplot(data=CoIa, aes(x=reorder(plotID, -dist), y= dist))+
  geom_point(aes(colour=herb_cov),size=4)+
  scale_color_viridis_c(option = "viridis",name="Herbaceous cover %",na.value = "grey60")+
  xlab('')+
  ylab('Residuals')+
  guides(size = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1))+
  scale_shape_discrete(name="")

p4 <- ggplot(data=CoIa, aes(x=reorder(plotID, -dist), y= dist))+
  geom_point(aes(colour=shannon),size=4)+
  scale_color_viridis_c(option="viridis", name="Shannon diversity")+
  xlab('')+
  ylab('Residuals')+
  guides(size = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1))+
  scale_shape_discrete(name="")

pp <- ggarrange(p1, p4, p2,p6, p3, p5, ncol = 2, nrow = 3)
ggexport(pp, filename = paste0("./R_output/ordination_manova/coinertia_avspectra/",
                               sites[i], "_COIA_resid.pdf"),width = 11,height = 8)


# extract plots trait data (trait loadings of the ordination conducted on the crossed covariance between spec and veg traits)
vegCoI_trait.coord <- data.frame(coia_veg.spec$co) %>%
  rownames_to_column("species")%>%
  mutate(species = make.unique(toupper(paste0(substr(species,1,3),
                                             substr(sapply(strsplit(species,"_"),"[",2),1,2)))))%>%
  mutate(leng= sqrt((abs(Comp1))^2+(abs(Comp2))^2))%>%
  arrange(desc(leng))

# rownames(vegCoI_trait.coord)[grep("Gal", rownames(vegCoI_trait.coord))]
# rownames(vegCoI_trait.coord) <- make.unique(toupper(paste0(substr(rownames(vegCoI_trait.coord),1,3),
#                substr(sapply(strsplit(rownames(vegCoI_trait.coord),"_"),"[",2),1,2))))

specCoI_trait.coord <- data.frame(coia_veg.spec$li) %>%
  rownames_to_column("wvl")%>%
  mutate(leng= sqrt((abs(Axis1))^2+(abs(Axis2))^2))%>%
  arrange(desc(leng))

saveRDS(coia_veg.spec, paste0("./R_output/ordination_manova/coinertia_avspectra/", sites[i],
                     "_COIA.rds"))

# sqrt((abs(0.0001097193))^2+(abs(-7.506199e-05))^2)

xy <- coia_veg.spec$lX[,1:2] %>% rownames_to_column("plotID")%>%
  cbind(coia_veg.spec$lY[,1:2]) %>%
  mutate(plotID= substr(plotID,nchar(plotID)-1,nchar(plotID)))

# Graph combining spec and veg trait covariation
# pdf(paste0("./R_output/ordination_manova/coinertia_avspectra/", sites[i], "_covar_plot.pdf"))
pp_co <- ggplot(data=xy, aes(label="plotID")) +
   geom_point(aes(x=AxcY1*4, y=AxcY2*4), col="goldenrod", cex=5)+
   geom_segment(aes(x=AxcX1, xend=AxcY1*4, y=AxcX2, yend=AxcY2*4),
                arrow=arrow(length=unit(0.2,"cm")), size=0.5) +
   geom_point(aes(x=AxcX1, y=AxcX2), col="springgreen4", cex=5)+
   geom_text(aes(x=AxcX1, y=AxcX2, label=plotID), cex=3, col="white")+
   geom_segment(data= specCoI_trait.coord[1:5,], aes(x=0, xend=Axis1*100, y=0, yend=Axis2*100),
                color="goldenrod", arrow=arrow(length=unit(0.3, "cm")), size=1) +
   geom_segment(data= vegCoI_trait.coord[1:5,], aes(x=0, xend=Comp1*50, y=0, yend=Comp2*50),
                color="springgreen4", arrow=arrow(length=unit(0.3,"cm")), size=1,) +
   geom_label(data=vegCoI_trait.coord[1:5,], aes(x= Comp1*50, y= Comp2*50,label=species),
              color="springgreen4", size=4, nudge_y=0.02)+
   geom_label(data=specCoI_trait.coord[1:5,], aes(x= Axis1*100, y= Axis2*100,label=wvl),
              color="goldenrod", size=4,nudge_x=0.15)+
   # geom_hline(yintercept=0, linetype="dotted") +
   # geom_vline(xintercept=0, linetype="dotted") +
   xlab(paste0("Axis1 (", round(eigen$relative_variance[1]*100,1),"%)"))+
   ylab(paste0("Axis2 (", round(eigen$relative_variance[2]*100,1),"%)"))+
   theme_bw()+ ggtitle(sites[i],subtitle = paste0("RV = ", round(coia_randu[[i]]$obs,2),
                              ", P = ", round(coia_randu[[i]]$pvalue,3)))
ggexport(pp_co, filename = paste0("./R_output/ordination_manova/coinertia_avspectra/",
                               sites[i], "_covar_plot.pdf"))

}

saveRDS(eigen, "./R_output/ordination_manova/coinertia_avspectra/rel_var_eigenvalues.rds")
saveRDS(coia_randu, "./R_output/ordination_manova/coinertia_avspectra/coinertia_MonteCarlo.rds")

coia_randu <- readRDS("./R_output/ordination_manova/coinertia_avspectra/coinertia_MonteCarlo.rds")
coia_randu[[1]]


### compile in table
coia_randu[[1]]$pvalue

xx <- glimpse(sapply(coia_randu, "[", 1))
y <- glimpse(sapply(coia_randu, "[", 4))
yy <- as.data.frame(do.call(rbind,y)) 
z <- glimpse(sapply(coia_randu, "[", 5))
zz <- as.data.frame(do.call(rbind,z)) %>% rename(p_val=V1)


xxx <- as.data.frame(do.call(rbind,xx)) %>% 
  bind_cols(siteID=sites)%>% 
  select(2,1)%>%
  rename(RV_obs=V1) %>%
  bind_cols(zz,yy)%>% remove_rownames()

write.csv(xxx, "./R_output/ordination_manova/coinertia_MonteCarlo.csv", row.names = F)

##### Rel Variance explained

var_dat <- bind_rows(eigen) %>% rownames_to_column("axis")%>%
  mutate(siteID = rep(good_sites, each=5))%>%
  mutate(axis=substr(axis,1,7))
write.csv(var_dat, "./R_output/ordination_manova/coinertia_variance_explained.csv", row.names = F)


### Number of obs per site
ls <- list.files("./R_output/ordination_manova/coinertia_avspectra/",pattern = "COIA.rds", full.names = T)
nobs <- data.frame(matrix(data = NA,nrow = length(ls),ncol = 2,dimnames = list(NULL, c("Site_ID","no_obs"))))

nobs$Site_ID <- substr(ls,nchar(ls)-12,nchar(ls)-9)

for(i in 1:length(ls)){
  da <-  readRDS(ls[i])
  nobs$no_obs[i] <- nrow(da$lX)
}

write.csv(nobs, "./pub/no_obs_coinertia.csv", row.names = F)



