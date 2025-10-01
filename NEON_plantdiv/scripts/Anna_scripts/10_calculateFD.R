library(FD)
library(tidyverse)
library(reshape2)

dat_save <- read.csv("./R_output/trait/traits_combined.csv")
dat <- unique(dat_save)
div_metrics <- read.csv("./R_input/model_input/NEON_diversity.csv")

### plants not in phylo (mosses, unknowns) unlikely to have traits measured
div_save <- read.csv("./R_output/inv_compiled/tree_cover_selectedplots_token.csv", row.names = 1)

div <- div_metrics[,c(1,4,5)] %>% left_join(div_save)

### LMA outliers from Samantha
out <- read.csv("~/Documents/CABO/NEON/R_input/plants/cfc_LMA_in_outliers.csv")
out$IDs <- paste(sapply(strsplit(as.character(out$plotID),"\\."),"[",1),
                 substr(out$sampleID,13,16), sep = "_") 

out$IDs <- gsub("CLIP","",out$IDs)

dat$IDs <- paste(dat$plotID,dat$taxonID, sep = "_")
dat[dat$IDs %in% out$IDs,]$leafMassPerArea <- NA

###### Sites needed 
(loci <- unique(substr(div_metrics$plotID,1,4)))

#### check remarks
unique(dat$remarks)
table(is.na(dat$remarks))
test <- dat %>% filter(!(is.na(remarks)), remarks>0)

#### select sites with diversity data
dat2 <- dat %>% mutate_if(is.factor, as.character) %>%
  select(-extractChlAConc, -extractChlBConc, -extractCarotConc, 
         -solventVolume, -PigmentFreshMass, -remarks) %>%
  mutate(species=paste(sapply(strsplit(scientificName," "),"[",1), 
                       sapply(strsplit(scientificName," "),"[",2),sep="_"))%>%
  select(1:5,36,37,6:35)%>%
  filter(siteID %in% loci)%>% 
  mutate(ChlAB_mg_g = ChlA_mg_g + ChlB_mg_g)%>% 
  mutate(ChlAB_mg_m2 = ChlA_mg_m2 + ChlB_mg_m2)

table(dat2$site %in% loci)
unique(div$site)[!(unique(div$site) %in% dat2$site)]

dat2[grepl("Herb", dat2$sampleType),]$taxonID <- "HERB"
dat2[grepl("Herb", dat2$sampleType),]$species <- "Herb"


####### Test value ranges
names(dat2)
testdf <- melt(dat2[,c(6,20:26,30,32:39)], id.vars = "IDs")
summary(dat2[,c(20:26,30,32:39)])

pdf("./R_output/trait/trait_ranges.pdf",
    height = 8, width = 10)
ggplot(data = testdf, aes(x=variable, y=value)) + 
  geom_boxplot() +  
  # geom_jitter()+
  facet_wrap( ~ variable, scales="free")
dev.off()

#############
### exclude implausible values... but mg.m2 implausibles are ok at mg.g
### eyeballed from boxplot for now, next calcualte maybe ymol
# names(dat2)

### Work with CHLa+b and CAR mg.g those are within range see Zhihui
# datx <- dat2 %>% 
#   mutate(ChlA_mg_m2=replace(ChlA_mg_m2, ChlA_mg_m2>500, NA))%>% 
#   mutate(ChlA_mg_g=replace(ChlA_mg_g, ChlA_mg_g>5.1, NA)) %>% 
#   mutate(ChlB_mg_g=replace(ChlB_mg_g, ChlB_mg_g>1.7, NA)) %>% 
#   mutate(ChlB_mg_m2=replace(ChlB_mg_m2, ChlB_mg_m2>150, NA)) %>% 
#   mutate(Car_mg_m2=replace(Car_mg_m2, Car_mg_m2>120, NA))%>% 
#   mutate(Car_mg_g=replace(Car_mg_g, Car_mg_g>1.5, NA))

# sort(dat2$ChlA_mg_m2[which(dat2$ChlA_mg_g>5.1)])

# summary(datx[,c(21:25,30,32:37)])

### End limit range

### Calc FD with limit range
### Each trait needs to be measured at least once per species and site 
### for the obs to be included
datx <- dat2
# loc <- sort(unique(datx$site))

# names(datx)
trai <- c("N_perc_mean","C_perc_mean", "d15N_mean","d13C_mean", "Lignin_perc_mean", 
          "Cell_perc_mean", "leafMassPerArea","Car_mg_g", "ChlAB_mg_g")
z <- list()
outi <- logical(length(loci))

for(i in 1:length(loci)){
  x <- datx [datx$site==loci[i],]
  y <- as.data.frame(t(rowsum(x[, trai], x$species, na.rm=T)))
  yy <- y %>% na_if(0)
  z[[i]] <- yy
  out <- colnames(yy)[is.na(colSums(yy))]
  x2 <- x %>% filter (!(species %in%out))
  outi[i] <- nrow(x2)==0
}

names(z) <- loci

### save overview for Samantha  
# saveRDS(z, "./R_input/plants/foliar_triats_sp_sites_missing_vals.rds")
# test <- readRDS("./foliar_traits_sp_sites_missing_vals.rds")
# test[1]

out_locs <- paste(loci[outi],collapse = "|")
names(datx)

dat3 <- datx %>% filter(!(grepl(out_locs, plotID)))%>% select(1:19,all_of(trai))

datx$plotID[!(datx$plotID %in% dat3$plotID)]

### standardize
dat_stand <- dat3 %>% 
  mutate_at(vars(trai), scale)

mean(dat_stand$C_perc_mean, na.rm=T)
sd(dat_stand$ChlAB_mg_g, na.rm=T)


#### Select sampled sites
div2 <- div %>% mutate(siteID=substr(plotID,1,4))%>%
  filter(siteID %in% unique(dat_stand$siteID))%>% arrange(plotID) %>% 
  dplyr::rename(Herb=herb_cov) %>% select(siteID,plotID,plotyear_tree,tree_cov, everything())

div2[is.na(div2)] <- 0

####
# clip <- datx %>% filter(sampleType=="Herbaceous clip strip")
# tree <- datx %>% filter(sampleType=="Woody individual")

# ss <- unique(tree$species)
# colnames(div)[!(colnames(div) %in% ss)]
# colnames(div)[grepl("Alnus",colnames(div))]

### Calculate functional div metrics
### trait values from clipstrip kept as background
### Loop through plots

ploti <- div2$plotID

FD <- list()
for(h in 1:length(ploti)){
  info <- div2[h,1:4]
  dd <- div2[h, -c(1:4)]
  ddx <- dd %>% select(where(~ sum(.) != 0))
  
  divi <- bind_cols(info,ddx)
  
  # divi <- div2[h,-(which(div2[h,]==0))]
  sp_ploti <- colnames(divi)[-c(1:4)]
  ### tree chem
  # if(ploti[i] %in% dat_c$plotiID){
  #   chemi <- dat_c[dat_c$plotiID==ploti[i], ] 
  #   t <- print("chem vals ploti")
  # } else {
  chemi <- dat_stand[grepl(substr(ploti[h],1,4),dat_stand$site), ] 
  # t <- print("chem vals site")
  # }
  
  # names(chemi)
  chemi2 <- chemi %>% select(plotID,species,trai)
  ss <- unique(chemi$species)
  # traits <- colnames(chemi2)[3:ncol(chemi2)]
  
  # chemix <- list()
  
  ## matrix traits per species per plot
  chemix <- matrix(nrow = length(ss),ncol=ncol(chemi2)-1) 
  colnames(chemix)<- colnames(chemi2)[-1]
  # u <- dd[dd$plotiID==ploti[i], traits[1]]
  # str(u)
  # length(as.vector(u))>0
  
  for (j in 1:length(ss)){
    dd <- chemi2[chemi2$species==ss[j],]
    for(k in 1:length(trai)){ ### was traits
      u <- as.vector(is.na(dd[dd$plotiID==ploti[h],trai[k]]))
      if(length(u)==1&u==TRUE||length(u)==0){ ### for NA vals
        chemix[j,k+1] <- mean(dd[,trai[k]], na.rm=T)
      } else{
        chemix[j,k+1] <- mean(dd[dd$plotiID==ploti[h],trai[k]], na.rm=T)  
      }
      chemix[j,1] <- ss[j]
    }
  }
  chemi2[order(chemi2$species),1:5] ### compare to chemix
  
  # if(ploti[h] %in% dd$plotiID){
  #   chemix[[j]] <- dd[dd$plotiID==ploti[i],-1]
  # }else{
  #   chemix[[j]] <- dd %>%
  #     group_by(species)%>%
  #     select(-plotiID)%>%
  #     summarise_all(mean,na.rm=TRUE)
  # }
  
  ### chemix is matrix for first ploti at first site
  chemixx <- as.data.frame(chemix) %>% filter(species %in% names(divi)) %>%
    arrange(species) %>% mutate_all(as.character) %>% mutate_at(-1, as.numeric)
  
  chemixx <- chemixx[complete.cases(chemixx), ]
  n <- names(divi)[names(divi)%in% chemixx$species]
  
  # div2[div2$plotID=="CLBJ_001",]$Celtis_laevigata ### test
  
  divxFD <- divi %>% remove_rownames()  %>%
    column_to_rownames("plotID") %>% select(all_of(n))
  divxFD <- divxFD %>% select(sort(colnames(divxFD))) 
  
  chemixFD <- chemixx %>% remove_rownames() %>%
    column_to_rownames("species")
  
  if(length(divxFD)==1|length(divxFD)==0){
    FD[[h]] <- list(NULL)
  } else {
    FD[[h]] <- dbFD(chemixFD,divxFD, w.abun=T)
  }
  
  # test <- weighted.var(divxFD, chemixFD)
  # test <- wtd.var(divxFD, chemixFD)
  
  divxFD2 <- divxFD %>% mutate(FD_total = as.numeric(rowSums(divxFD)))%>%
    bind_cols(divi[,c("tree_cov", "Herb")],.name_repair = "minimal")
  
  names(divxFD2)[ncol(divxFD2)] <- "Herb_cover"
  
  # row.names(divxFD2) <- row.names(divxFD)
  write.csv(divxFD2,paste0("./R_output/trait/FD_div_tables/", ploti[h], "_div4FDFD_.csv"))
}

names(FD) <- ploti
saveRDS(FD, "./R_output/trait/FD_list.rds")

names(FD)

######### END ###########
