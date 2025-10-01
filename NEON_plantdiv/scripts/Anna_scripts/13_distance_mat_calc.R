library(tidyverse)

### scaling function between 0 and 1 for distance matrices
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}

### Plots used 
divi <- read.csv("./R_input/model_input/NEON_diversity.csv")
div <- divi %>% replace_na(list(todo="NA")) %>%
  mutate(siteID=substr(plotID,1,4), .before=1)%>%
  filter(!(grepl("ut", todo)))%>%
  filter(n_pix_NDVI>=50)

p <- div$plotID
s <- unique(div$siteID)

### Calculate distance matrices

### Average spectra per plot and site ###
ss <- readRDS("~/Documents/CABO/NEON/R_input/model_input/average_spectra/av_VNspectra.rds")

ss[(which(!(names(ss) %in% s),arr.ind = T))] <- NULL
table(sites==names(ss))

#### start with ABBY
spec_dist_mats <- list()

for (i in 1:length(s)){
  spec <- ss[[i]] %>% filter(plotID %in%p) %>% arrange(plotID)%>%
    column_to_rownames("plotID")%>% select(-c(1:4))
  
  spec_mat <- as.matrix(spec)
  dist_euclid <- dist(spec_mat, method = "euclid") 
  # dist_spec <- dist(spec_mat, method = "manhattan") 
  
  dist_spec_mat <- as.matrix(dist_euclid) 
  dist_spec_mat01 <- scale01(dist_spec_mat)  
  dist_spec01 <- as.dist(dist_spec_mat01)
  
  spec_dist_mats[[i]] <- dist_spec01
  names(spec_dist_mats)[[i]]<- s[i]
}

saveRDS(spec_dist_mats, "./R_input/model_input/spec_dist_euclid01.rds")
    
##### Vegetation distances  
df <- read.csv("./R_output/inv_compiled/veg_cover_compiled_token.csv")

### Veg dataframe including NA species, moss, lichen, etc
dfi <- df %>% select(-c(2:8)) %>%
  select(-c(wood_NA, water_NA,standingDead_NA, litter_NA, rock_NA, soil_NA,scat_NA)) %>%
  filter(plotID %in% p)

###########
veg_dist_mats <- list()

for (i in 1:length(s)){
  veg <- dfi %>% filter(grepl(s[[i]], plotID))%>% filter(plotID %in%p) %>% 
    arrange(plotID) %>% column_to_rownames("plotID") %>% 
    select_if(~all(!is.na(.)))
  
  veg_mat <- as.matrix(veg)
  veg_mat_hell <- decostand(veg_mat, method="hellinger") 
  
  dist_euclid <- dist(veg_mat_hell, method = "euclid") 
  
  dist_veg_mat <- as.matrix(dist_euclid) 
  dist_veg_mat01 <- scale01(dist_veg_mat)  
  dist_veg01 <- as.dist(dist_veg_mat01)
  
  veg_dist_mats[[i]] <- dist_veg01
  names(veg_dist_mats)[[i]]<- s[i]
}

saveRDS(veg_dist_mats, "./R_input/model_input/veg_dist_hellinger01.rds") 

veg_dist_mats <- readRDS("./R_input/model_input/veg_dist_hellinger01.rds") 

### check no of obs
lapply(veg_dist_mats, FUN =  function(x) nrow(as.matrix(x)))


#### END ###


