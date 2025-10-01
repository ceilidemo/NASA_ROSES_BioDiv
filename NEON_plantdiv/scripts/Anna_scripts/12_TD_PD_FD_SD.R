library(ape)
library(vegan)
library(plyr)
library(picante)
library(tidyverse)

info <- read.csv("./R_input/NEON_plantcodes_lookup.csv", row.names = 1)
phy <- readRDS("./R_output/phylo/allplants_scenario2_tree_100.rds")


#### OLD SD" for update see end of the script

# ## Spectral diversity
# SD <- read.csv("~/Documents/CABO/NEON/R_input/model_input/NEON_diversity_bestmasks.csv")
# # 
# # SD[SD$plotID=="DSNY_003",]
# 
# SD_df <- SD %>% select(c(1:24,72:79)) %>%
#   select(-c(vegtype, site,State)) %>% drop_na(SDiv_alpha)


### Vegdata
df <-  read.csv("./R_output/inv_compiled/veg_cover_compiled_token.csv", check.names = F)

info[grepl("_NA", info$sp_matched),]$sp_matched
names(df)[grepl("_NA", names(df))]
names(df)[grepl("^X2", names(df))]

### Veg dataframe including NA species, moss, lichen, etc
row.names(df) <- df$plotID
veg_df <- df %>% select(-c(1:8)) %>%
  select(-c(wood_NA, water_NA,standingDead_NA, litter_NA, rock_NA, soil_NA,scat_NA)) 

### Veg dataframe phylomatch
names(veg_df)[!(names(veg_df)%in% phy$scenario.2$run.1$tip.label)]

veg_df_phylo <- veg_df %>% select(-names(veg_df)[!(names(veg_df)%in% phy$scenario.2$run.1$tip.label)])
names(veg_df_phylo)[!(names(veg_df_phylo)%in% phy$scenario.2$run.1$tip.label)]

# write.csv(veg_df_phylo, "./R_output/phylo/veg_df_phylo.csv")

#### Taxonomic diversity indices
# S <- specnumber(veg_df,MARGIN = 1)
SR <- apply(veg_df, 1, function(i) sum(i > 0))
shannon <- diversity(veg_df,index = c("shannon"))
simpson <- diversity(veg_df,index = c("simpson"))

SR_phylo <- apply(veg_df_phylo, 1, function(i) sum(i > 0))
shannon_phylo <- diversity(veg_df_phylo,index = c("shannon"))
simpson_phylo <- diversity(veg_df_phylo,index = c("simpson"))

res <- as.data.frame(SR) 
res <- as.data.frame(cbind(SR, shannon, simpson, SR_phylo,shannon_phylo, simpson_phylo))
res <- res %>% rownames_to_column(var="plotID") 

all_df <- df[,1:7] %>% left_join(res) %>% left_join(SD_df, by="plotID") %>% 
  drop_na(SDiv_alpha)

anyNA(SD_df$SDiv_alpha)

write.csv(all_df, "./R_input/model_input/NEON_diversity.csv", row.names = F)

######################
### FD metrics #######
all_df <- read.csv("./R_input/model_input/NEON_diversity.csv")

### check FD_div_tables if % veg in FD is representative
FD_tab <- readRDS("./R_output/trait/FD_list.rds")

FDs <- as.data.frame(unlist(FD_tab))
names(FDs) <- "value"
FDsx <- FDs %>% mutate(plotID=substr(row.names(FDs),1,8)) %>%
  mutate(metric=substr(row.names(FDs),10,nchar(row.names(FDs))))%>%
  mutate(metric= gsub("\\.[A-Z]{4}.[0-9]{3}$","", metric)) %>%
  spread(key = metric, value = value)


FDinfo <- list.files("./R_output/trait/FD_div_tables/", pattern= ".csv", full.names = T)
FDinfo_combi <- map(FDinfo, read.csv)

FD_infi <- bind_rows(FDinfo_combi)
FD_info <- FD_infi %>% select(X, FD_total, tree_cov,Herb_cover, everything())%>% 
  rename(plotID=X)
  
res <- FD_info[,1:4] %>% right_join(FDsx)
all_df2 <- all_df %>% left_join(res)

# write.csv(res5, "./R_input/model_input/NEON_diversity.csv", row.names = F)

##### PD METRICS RUN REMOTELY
#### BASED on 1m inventories ########## 
#### Phylogenetic diveristy metrics
xx <- readRDS("./R_output/phylo/phylo_metrics100.rds")
class(xx[[1]])

yy <- lapply(xx, as.matrix)
head(yy[[1]])

#### Mean phylo metrics
datx <- list()
for (i in 1:length(yy)){
  datx[[i]] <- yy[[i]][,-1]
  row.names(datx[[1]]) <- yy[[i]][,1]
  class(datx[[i]]) <- "numeric"
}

dati <- Reduce("+",datx)
dati <- dati/100

phylo_dat <- as.data.frame(dati) %>% rownames_to_column(var = "plotID")

phylo_dat$plotID[!(phylo_dat$plotID %in% all_df$plotID)]

all_df3 <- all_df2 %>% left_join(phylo_dat)
names(all_df3)

write.csv(all_df3, "./R_input/model_input/NEON_diversity.csv", row.names = F)

#####################################
#### Update Spectral diversity
all_df3 <- read.csv("./R_input/model_input/XX_NEON_diversity.csv")
df <- all_df3 %>% select(-c(SDiv_alpha_log, SDiv_alpha, LCSD_beta, LCSS_beta, 
                            n_pix, n_total))

### Spectral diversity 
fls <- list.files("./R_output/specdiv_goodplots/",pattern = "SDiv.csv",
                  full.names = T,recursive = T)
SD <- map(fls,read.csv)

siteID <- substr(fls, nchar(fls)-12, nchar(fls)-9)
SD_x <- mapply(cbind, SD, "siteID"= siteID, SIMPLIFY=F)

# test <- SD_x[[1]] ### thinking about PCs but number of PCs differs
# colnam_in <- colnames(test)[c(1:11, (ncol(test)-4): ncol(test))]

SD_xx <- lapply(SD_x, function(x) x[,!(grepl("PC",colnames(x)))]) ### without PCs
# SD_xx <- lapply(SD_x, function(x) x[,colnames(x) %in% colnam_in])
# lapply(SD_xx, ncol)

SD_df <- do.call(rbind.data.frame, SD_xx)
ids <- sapply(strsplit(as.character(SD_df$plotID),"_"),"[",2)
ids <- str_pad(ids, 3, pad = "0")

SD_df2 <- SD_df %>% 
  mutate(plotID=paste(siteID, ids,sep="_"))%>% 
  select(11,1:10) %>% arrange(plotID) 

names(SD_df2)
summary(SD_df2$n_pix_NDVI)

df2 <- merge(df, SD_df2[,-1],by = "plotID", all.x = T)

write.csv(df2, "./R_input/model_input/NEON_diversity.csv")
#### END####

