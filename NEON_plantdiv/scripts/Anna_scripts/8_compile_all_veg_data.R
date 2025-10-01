library(tidyverse)
library(plyr)

sites <- read.csv("./orig_data/NEON_struct-woody-plant/stackedFiles/vst_perplotperyear.csv")
info <- read.csv("./R_input/NEON_plantcodes_lookup.csv")
phy <- readRDS("./R_output/phylo/allplants_scenario2_tree_100.rds")

##### Trees 
df_trees <- read.csv("./R_output/inv_compiled/tree_cover_allplots_token.csv",row.names = 1)

## check names and phylogeny
table(names(df_trees) %in% info$sp_matched)
names(df_trees)[!(names(df_trees) %in% info$sp_matched)]

names(df_trees)[grepl("Betula_glandulosa.nana", names(df_trees))]
names(df_trees) <- gsub("Betula_glandulosa.nana", "Betula_nana", names(df_trees))

names(df_trees)[grepl("SARA", names(df_trees))]
names(df_trees) <- gsub("Sambucus_racemosa.nigra", "Sambucus_racemosa", names(df_trees))

df_trees <- df_trees %>% 
  mutate(tree_cov=rowSums(.))%>% select(tree_cov,everything())

### Phylogeny ######
names(df_trees)[!(names(df_trees) %in% phy$scenario.2$run.1$tip.label)]

#####################
#### Herbaceous match with df_trees and make new phylo 
df_herbs <- read.csv("./R_output/inv_compiled/herb_cover_allplots.csv", row.names = 1)

## check names and phylogeny
table(names(df_herbs) %in% info$taxonID)
miss <- names(df_herbs)[!(names(df_herbs) %in% info$taxonID)]

miss[!(grepl("^X2",miss))]

#### rename columns
nami <- names(df_herbs)[!(grepl("^X2",names(df_herbs)))]

### mutiple abbrevs for one taxon ID
# uni <- info[!(duplicated(info$sp_matched)),]
# uni <- uni[complete.cases(uni$taxonID),]
# uni_sel <- uni[uni$taxonID %in% nami,]

for (i in 1:length(nami)){
  names(df_herbs)[!(grepl("^X2",names(df_herbs)))][i] <- info[info$taxonID==names(df_herbs)[!(grepl("^X2",names(df_herbs)))][i],]$sp_matched 
}

names(df_herbs)[duplicated(names(df_herbs))]

test <- df_herbs[,256:257] 
test <- test[rowSums(test)>0,]
names(test) <- c("A", "A")
testi <- as.data.frame(sapply(split.default(test, names(test)), rowSums, na.rm = TRUE))

### could be saved 
df_herbi <- as.data.frame(sapply(split.default(df_herbs, names(df_herbs)), rowSums, na.rm = TRUE))

#############
### Total herb cover is way over 100 per plot -> sale to 100% cover
df_herbix <- df_herbi %>% mutate(tot_herb=rowSums(.)) 

df_herbix <- df_herbix %>% select(tot_herb, everything())%>%
  rownames_to_column(var = "plotyear") %>%
  mutate(plotyear=paste0(substr(plotyear,1,4),substr(plotyear,10,13), substr(plotyear,5,9)))

funx <- function(x) (x/(df_herbix$tot_herb/100))
df_herb_scale <- df_herbix %>% mutate_at(.vars = -c(1,2), .funs = funx)

df_herb_scale <- df_herb_scale %>% 
  mutate(plotID=substr(plotyear,1,8)) %>% select(plotID, everything())

rowSums(df_herb_scale[,-c(1:3)])


#### Select only plots where tree_cover is representative ### this can be done later
# veg_class <- sites %>% select(plotID, nlcdClass) %>% unique()
# anyDuplicated(veg_class$plotID)
# unique(veg_class$nlcdClass)

# ggplot(tree_tot, aes(tree_cov, nlcdClass, fill=factor(nlcdClass))) +
#   geom_boxplot(show.legend = F)+geom_point(show.legend = F)
# 
# outi <- tree_tot %>% filter(grepl("orest", nlcdClass), tree_cov<10) %>% arrange(plotyear)%>%
#   filter(plotID %in% herbIDs)
  

### Plots with spectral data
spec <- read.csv("~//Documents/CABO/NEON/R_input/model_input/NEON_diversity_bestmasks.csv")

###### DECIDE WHICH PLOTYEARS TO USE #######
h_2018 <- df_herb_scale[grepl("2018", df_herb_scale$plotyear),]

spec$plotID[!(spec$plotID %in% h_2018$plotID)]
# unique(df_herb_scale$plotID)[!(unique(df_herb_scale$plotID) %in% spec$plotID)]

df_trees <- df_trees %>% 
  rownames_to_column(var="plotyear") %>%
  mutate(plotID=substr(plotyear,1,8))%>% 
  mutate(year= substr(plotyear,10,13))%>% select(plotyear,plotID,year, everything())

### select tree inventories that match best
### start with exceptions
length(sort(unique(spec$site)) [sort(unique(spec$site)) %in%unique(substr(df_trees$plotID,1,4))])

alt <- df_trees %>% filter(plotyear %in% c("ABBY_012_2017","ABBY_023_2017"))  

ABBY <- df_trees %>% filter(grepl("ABBY", plotID)) %>%
  filter(year=="2019") %>% filter(!(plotyear %in% c("ABBY_012_2019","ABBY_023_2019")))%>%
  bind_rows(alt) %>% arrange(plotID)

BART <- df_trees %>% filter(grepl("BART", plotID)) %>%
  filter(year=="2016") %>% arrange(plotID)

BONA <- df_trees %>% filter(grepl("BONA", plotID)) %>%
  filter(year=="2017") %>% arrange(plotID)

CLBJ <- df_trees %>% filter(grepl("CLBJ", plotID)) %>%
  filter(year=="2016") %>% arrange(plotID)

DEJU <- df_trees %>% filter(grepl("DEJU", plotID)) %>%
  filter(year=="2019") %>% arrange(plotID)

DSNY <- df_trees %>% filter(grepl("DSNY", plotID)) %>%
  filter(year=="2018") %>% arrange(plotID)

GRSM <- df_trees %>% filter(grepl("GRSM", plotID)) %>%
  filter(year=="2016") %>% arrange(plotID)

alt <- df_trees %>% filter(plotyear %in% c("GUAN_003_2017","GUAN_015_2017","GUAN_016_2017")) 
GUAN <- df_trees %>% filter(grepl("GUAN", plotID)) %>%
  filter(year=="2018") %>% filter(!(plotyear %in% c("GUAN_003_2018","GUAN_015_2018","GUAN_016_2018")))%>%
  bind_rows(alt) %>% arrange(plotID)

alt <- df_trees %>% filter(plotyear %in% c("HARV_020_2015","HARV_024_2015")) 
HARV <- df_trees %>% filter(grepl("HARV", plotID)) %>%
  filter(year=="2019") %>% filter(!(plotyear %in% c("HARV_020_2019","HARV_024_2019")))%>%
  bind_rows(alt) %>% arrange(plotID)

alt <- df_trees %>% filter(plotyear %in% c("HEAL_021_2015","HEAL_022_2015","HEAL_026_2017")) 
HEAL <- df_trees %>% filter(grepl("HEAL", plotID)) %>%
  filter(year=="2016") %>% filter(!(plotyear %in% c("HEAL_021_2016","HEAL_022_2016","HEAL_026_2016")))%>%
  bind_rows(alt) %>% arrange(plotID)

JERC <- df_trees %>% filter(grepl("JERC", plotID)) %>%
  filter(year=="2019")

alt <- df_trees %>% filter(plotyear %in% c("JORN_013_2019","JORN_017_2019")) 
JORN <- df_trees %>% filter(grepl("JORN", plotID)) %>%
  filter(year=="2018") %>% filter(!(plotyear %in% c("JORN_013_2018","JORN_017_2018")))%>%
  bind_rows(alt) %>% arrange(plotID)

KONZ <- df_trees %>% filter(grepl("KONZ", plotID)) %>%
  filter(year=="2018")

LAJA <- df_trees %>% filter(grepl("LAJA", plotID)) %>%
  filter(year=="2016")

LENO <- df_trees %>% filter(grepl("LENO", plotID)) %>%
  filter(year=="2017")

MLBS <- df_trees %>% filter(grepl("MLBS", plotID)) %>%
  filter(year=="2018")

alt <- df_trees %>% filter(plotyear %in% c("MOAB_044_2016","MOAB_045_2016")) 
MOAB <- df_trees %>% filter(grepl("MOAB", plotID)) %>%
  filter(year=="2018") %>% filter(!(plotyear %in% c("MOAB_044_2018","MOAB_045_2018")))%>%
  bind_rows(alt) %>% arrange(plotID)

NIWO <- df_trees %>% filter(grepl("NIWO", plotID)) %>%
  filter(year=="2019")

alt <- df_trees %>% filter(plotyear %in% c("OSBS_027_2015")) 
OSBS <- df_trees %>% filter(grepl("OSBS", plotID)) %>%
  filter(year=="2016") %>% filter(!(plotyear %in% c("OSBS_027_2016")))%>%
  bind_rows(alt) %>% arrange(plotID)

RMNP <- df_trees %>% filter(grepl("RMNP", plotID)) %>%
  filter(year=="2018")

alt <- df_trees %>% filter(plotyear %in% c("SJER_012_2016","SJER_045_2016")) 
SJER <- df_trees %>% filter(grepl("SJER", plotID)) %>%
  filter(year=="2019") %>% filter(!(plotyear %in% c("SJER_012_2019","SJER_045_2019")))%>%
  bind_rows(alt) %>% arrange(plotID)

SOAP <- df_trees %>% filter(grepl("SOAP", plotID)) %>%
  filter(year=="2019")

alt <- df_trees %>% filter(plotyear %in% c("SRER_043_2016","SRER_047_2016","SRER_052_2017")) 
SRER <- df_trees %>% filter(grepl("SRER", plotID)) %>%
  filter(year=="2019") %>% filter(!(plotyear %in% c("SRER_043_2019","SRER_047_2019","SRER_052_2019")))%>%
  bind_rows(alt) %>% arrange(plotID)
  
TALL <- df_trees %>% filter(grepl("TALL", plotID)) %>%
  filter(year=="2015")

UKFS <- df_trees %>% filter(grepl("UKFS", plotID)) %>%
  filter(year=="2018")

WREF <- df_trees %>% filter(grepl("WREF", plotID)) %>%
  filter(year=="2019")

DELA <- df_trees %>% filter(grepl("DELA", plotID)) %>%
  filter(year=="2019")

alt <- df_trees %>% filter(plotyear %in% c("YELL_002_2018","YELL_006_2018","YELL_013_2018")) 
YELL <- df_trees %>% filter(grepl("YELL", plotID)) %>%
  filter(year=="2019") %>% filter(!(plotyear %in% c("YELL_002_2019","YELL_006_2019","YELL_013_2019")))%>%
  bind_rows(alt) %>% arrange(plotID)

tree_sel <- bind_rows(ABBY, BART,BONA, CLBJ, DEJU, DELA, DSNY,GRSM,GUAN,HARV,HEAL,JERC,JORN,KONZ, 
                      LAJA,LENO,MLBS,MOAB,NIWO, OSBS,RMNP, SJER,SOAP,SRER, TALL, UKFS, WREF, YELL)

table(rowSums(tree_sel[,-c(1:4)]) - tree_sel$tree_cov)

tree_seli <- tree_sel %>% dplyr::rename(plotyear_tree=plotyear)%>% select(-year)


### save selected tree inventories for FD calc
# tt_sel <- tree_seli %>% select(where(~ any(. != 0)))
# write.csv(tt_sel, "./R_output/inv_compiled/tree_cover_selectedplots_token.csv")


h_2018 <- h_2018 %>% rename(plotyear_herb=plotyear)
tree_tot <- tree_seli %>% select(plotID, tree_cov)

# tree_tot$plotID[!(tree_tot$plotID %in% h_2018$plotID)]
table(h_2018$plotID %in% tree_tot$plotID)

notrees <- as.data.frame(h_2018$plotID[!(h_2018$plotID %in% tree_tot$plotID)])
names(notrees) <- "plotID"

tree_tot2 <- tree_tot %>% bind_rows(notrees) %>% arrange(plotID)%>% replace_na(list(tree_cov=0)) %>% 
  filter(plotID %in% h_2018$plotID)

table(h_2018$plotID %in% tree_tot2$plotID)
table(tree_tot2$plotID %in% h_2018$plotID)

#################################################
#### Scale herb cover to area not covered by trees

df_herb_scale2 <- right_join(tree_tot2, h_2018)%>%
  mutate(herb_cov=100-tree_cov, .after=2) 

funi<- function(x) (x*(df_herb_scale2$herb_cov/100))
df_herb_scale3 <- df_herb_scale2 %>% mutate_at(.vars = -c(1:5), funi)

rowSums(df_herb_scale3[,-c(1:6)]) - df_herb_scale2$herb_cov

############
####### Combine with tree data
### https://stackoverflow.com/questions/16018863/combine-data-frames-summing-up-values-of-identical-columns-in-r

### test with HARV
anyDuplicated(df_herb_scale3$plotID)
# # 
# HARV_t <- tree_seli %>% filter(grepl("HARV", plotID)) %>%
#   select(-c(1,3)) %>% select(where(~ any(. != 0)))
# 
# HARV_h <- df_herb_scale3 %>% filter(grepl("HARV", plotID)) %>%
#   select(-c(2:5)) %>%select(where(~ any(. != 0)))
# 
# pp <- cbind(names=c(HARV_t$plotID, HARV_h$plotID),
#             rbind.fill(list(HARV_t[,-c(1)],HARV_h[,-c(1)])))
# pp_all <- ddply(pp, .(names), function(x) colSums(x[,-1], na.rm = TRUE))


### ALL
tt <- tree_seli %>% select(-c(1,3)) %>% select(where(~ any(. != 0)))
hh <- df_herb_scale3 %>% select(-c(2:5)) %>%select(where(~ any(. != 0)))

combi <- cbind(names=c(tt$plotID, hh$plotID),
            rbind.fill(list(tt[,-c(1)], hh[,-c(1)])))
combi_all <- ddply(combi, .(names), function(x) colSums(x[,-1], na.rm = TRUE))
any(colSums(combi_all[,-1])==0)

combi_all2 <- combi_all %>% dplyr::rename(plotID=names) %>%
  left_join(df_herb_scale3[,1:5]) %>% left_join(tree_seli[,1:2])%>%
  select(plotID, tree_cov,herb_cov, plotyear_herb, plotyear_tree,tot_herb,
         everything(.))

### add info about to exclude
### This could be checked again
info <- read.csv("./R_input/tree_inventories_missing.csv")
combi_all3 <- info %>% right_join(combi_all2) %>% arrange(plotID) %>%
  filter(plotID %in% spec$plotID) %>% select(1:8, sort(colnames(combi_all3)[-c(1:8)]))

write.csv(combi_all3, "./R_output/inv_compiled/veg_cover_compiled_token.csv", row.names=F)

### Test phylo 
veg <- read.csv("./R_output/inv_compiled/veg_cover_compiled_token.csv")
names(veg)[-c(1:8)][!(names(veg)[-c(1:8)]%in%info$sp_matched)]
names(veg)[-c(1:8)][!(names(veg)[-c(1:8)]%in%info$sp_matched)][1:10]

names(vegi)[-c(1:8)][!(names(vegi)[-c(1:8)]%in%info$sp_matched)][1:10]
vegi <- veg %>% dplyr::rename("Arctostaphylos_uva-ursi" = "Arctostaphylos_uva.ursi",
                       "Athyrium_filix-femina" ="Athyrium_filix.femina",
                       "Capsella_bursa-pastoris"="Capsella_bursa.pastoris",
                       "Clermontia_montis-loa" ="Clermontia_montis.loa",
                       "Crataegus_crus-galli" ="Crataegus_crus.galli",
                       "Hypericum_crux-andreae" ="Hypericum_crux.andreae",
                       "Macfadyena_unguis-cati" = "Macfadyena_unguis.cati",
                       "Pithecellobium_unguis-cati"="Pithecellobium_unguis.cati",
                       "Smilax_bona-nox"="Smilax_bona.nox",
                       "Vaccinium_vitis-idaea"= "Vaccinium_vitis.idaea")

# phy$scenario.2$run.1$tip.label[grepl("Arctostaphylos", phy$scenario.2$run.1$tip.label)]
names(vegi)[-c(1:8)][!(names(vegi)[-c(1:8)] %in% phy$scenario.2$run.1$tip.label)]
vegi <- veg_compiled %>% select(where(~ any(. != 0)))
write.csv(vegi, "./R_output/inv_compiled/veg_cover_compiled_token.csv", row.names=F)


#################################
### Test missing plotID from token
veg_compiled <- read.csv("./R_output/inv_compiled/veg_cover_compiled_token.csv")

veg_compiled[grepl("CLBJ_036",combi_all3$plotID),1:5]
combi_all3[grepl("CLBJ_036",combi_all3$plotID),1:5]


#### END 