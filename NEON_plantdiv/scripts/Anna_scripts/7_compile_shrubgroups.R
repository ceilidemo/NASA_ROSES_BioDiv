library(tidyverse)

### shrubs < 3m also have an entry in the herbaceous data, HARV 2 shrubgroups are < 3m
####
#### OTHER SITES: left out for now, species seem to all occur in the herb inventories

indiv_save <- read.csv("//Volumes/NEON_2021/NEON_SD/orig_data/NEON_struct-woody-plant/stackedFiles/vst_shrubgroup.csv")

unique(indiv_save$plotID)

### shrub groups are not mapped, subplot is recorded
# map <- read.csv("//Volumes/NEON_2021/NEON_SD/orig_data/NEON_struct-woody-plant/stackedFiles/vst_mappingandtagging.csv")

indiv <- indiv_save %>% filter(siteID=="HARV") 
indiv <- indiv[grepl("^2019", indiv$date),]


shrub <- indiv %>% mutate_at(1:15, as.character) %>%
  mutate(loc=paste(indiv$plotID,indiv$subplotID, indiv$nestedSubplotID, sep="_"),
         site = substr(as.character(indiv$plotID),1,4),
         cover_m2 = canopyArea*(volumePercent/100)*(livePercent/100),
         subplotyear = paste(loc,sapply(strsplit(as.character(date), "-"),"[",1),sep = "_"),
         plotyear = paste(plotID,sapply(strsplit(as.character(date), "-"),"[",1),sep = "_"),
         subplotdate = paste(loc,date,sep = "_"))
  # select(25:30,1:24)

length(unique(shrub$subplotdate))

# shrubi <- shrub %>%
#   pivot_wider(id_cols = c(site, subplotyear, loc, plotyear,subplotdate, plotID), names_from= taxonID,
#               values_from = cover_m2,values_fn = list(cover_m2=sum))


### This gets complex when there are more sp. per groups VolumePercent gives % per sp. per group
### calculate % cover, canopy area is m2 


shrub_perc <- shrub %>% 
  select (site, subplotyear, loc, plotyear,subplotdate, plotID,cover_m2,taxonID)%>% 
  group_by(plotyear,taxonID) %>% 
    summarise(tot_cov=sum(cover_m2))

shrub_perc <- shrub_perc %>% mutate(cover_perc = tot_cov/400*100)%>% arrange(plotyear)

shrub_perc2 <- shrub_perc %>%
  pivot_wider(id_cols = c(plotyear), names_from= taxonID,
              values_from = cover_perc)

write.csv(shrub_perc2, "./R_output/inv_compiled/shrub_groups_compiled.csv", row.names = F)

# ### Check shrubs
# plants <- read.csv("./R_input/NEON_plantcodes_lookup.csv")
# colnames(df_shrubs)[-1][!(colnames(df_shrubs)[-1] %in% plants$taxonID)]
# names(df_shrubs)[names(df_shrubs) %in% "SPDOD"] <- "SPDO"
# names(df_shrubs)[1] <- "plotID"
# df_shrubs <- df_shrubs %>% select(plotID,order(colnames(.))) %>% arrange(plotID)
# write.csv(df_shrubs,"./R_output/biodiversity_compiled/shrub_groups_compiled.csv", row.names = F)

######################################
### Transform from m2 cover to % cover
# shrubi <- read.csv("./R_output/biodiversity_compiled/shrub_groups_compiled.csv", check.names = F)
info <- read.csv("./R_input/plants/stackedFiles_woody/vst_perplotperyear.csv")

#### look up plots needed
old_div <- read.csv("./R_input/model_input/NEON_diversity_allplots_plants.csv")
test <- unique(as.character(old_div[grepl("2018", old_div$locID),3]))

shrubi <- shrubi %>% mutate_at(1:5, as.character) %>%
  mutate(locID=substr(loc,1,8)) %>%
  select(c(1:5,104,6:103))

shrubx <- shrubi[shrubi$locID %in% test,]

shrubx <- shrubx %>%  mutate(year=substr(plotyear,nchar(shrubi$plotyear)-3, nchar(shrubi$plotyear)))%>%
  select(c(1:6,105,7:104))

# ### select closest sampling date to 2018 from duplicates
# dubs <- shrubi$loc[which(duplicated(shrubi$loc))]
# 
# # shrubi$plotyear[which(duplicated(shrubi$plotyear))]
# 
# ### keep uniques
# unis <- shrubi[!(shrubi$loc %in% dubs), ]
# length(unique(dubs)) 
# 133 + 308
# 
# ### find closest dates
# dubi <- shrubi[shrubi$loc %in% dubs, 1:6]
# dubi <- dubi[order(dubi$loc),]
# dubi$diff <- abs(2018 - as.numeric(dubi$year))
# 
# sel <- dubi %>% group_by(loc) %>% filter(diff==min(diff))
# anyDuplicated(sel$loc)
# 
# xx <- sel$loc[duplicated(sel$loc)]
# yy <- sel[sel$loc %in%xx,]
# 
# seli <- shrubi[shrubi$subplotyear %in% sel$subplotyear, ]
# test <- seli[seli$loc %in%xx,] %>%discard(~all(is.na(.))) 
# 
# combi <- bind_rows(unis,seli) %>% arrange(loc)
# 
# combi<- cbind(combi[,1:6],combi[,sort(names(combi)[7:104])])

### Plotsize
info <- info %>% mutate(plotyear=paste(plotID,substr(date,1,4),sep = "_"),
                        plotdate=paste(plotID,date,sep = "_"))
infi <- unique(info[,c(7,29,35,36)])

anyDuplicated(infi$plotdate)
z <- infi[infi$plotdate %in% infi$plotdate[duplicated(infi$plotdate)],]

nami <- sort(names(shrubx)[8:105])


### calculate % cover
shrub_perc <- shrubx %>% 
  mutate(plotdate = paste(substring(subplotdate,1,8),
                          substring(subplotdate,nchar(subplotdate)-10,nchar(subplotdate)),sep = "")) %>%
  select(106,1:105) %>%
  group_by(plotdate) %>% summarise_at(nami,sum, na.rm=T)%>%
  left_join(infi,by = "plotdate", all.x=T) %>%
  mutate(plotID=substr(plotdate,1,8))%>%
  select(1,100:103,2:99) 


# ### no info about plot size
# noinfo <- shrub_perc[is.na(shrub_perc$totalSampledAreaShrubSapling),]
# tt <- cbind(noinfo[,1:4], noinfo[,names(which(colSums(noinfo[,6:103])>0))])
# 
# ### with info about plot size
# shurb_perc2 <- shrub_perc[complete.cases(shrub_perc$totalSampledAreaShrubSapling),]
# shrub_perc3 <- shurb_perc2 %>% mutate_at(vars(6:103),funs(. / totalSampledAreaShrubSapling*100))
# 
# xx <- shrub_perc3$plotID[duplicated(shrub_perc3$plotID)]
# dubs <- shrub_perc3[shrub_perc3$plotID %in%xx,]
# 
# unis <- shrub_perc3[!(shrub_perc3$plotID %in%xx),]
# 
# ### select later date for duplicated
# dubsi <- dubs %>% 
#   mutate(year=substr(dubs$plotdate,10,13), loc=substr(dubs$plotdate,1,4)) %>%
#   mutate(diff= abs(2018 - as.numeric(year))) %>% select(1:5,104:106,6:103)
#   
# dubs_sel <- dubsi %>% group_by(plotID) %>% filter(diff==min(diff))
# 
# dd <- dubs_sel$plotID[duplicated(dubs_sel$plotID)]
# dupi <- dubs_sel[dubs_sel$plotID %in%dd,]
# dubs_sel2 <- dubs_sel[!(dubs_sel$plotID %in%dd),]
# 
# dupi_sel <-dupi[c(2,4),]
# 
# dubs_selx <- bind_rows(dubs_sel2,dupi_sel)
# anyDuplicated(dubs_selx$plotID)
# 
# # zz <- cbind(dubs[,1:4], dubs[,names(which(colSums(dubs[,5:102])>0))])
# 
# combi <- bind_rows(unis,dubs_selx)
# anyDuplicated(combi$plotID)
# 
# 
# #### to add noinfo plots
# noinfo2 <- noinfo %>% mutate(plotyear=substr(plotdate,1,13)) %>% 
#   select(4:103) %>%
#   left_join(infi,by="plotyear", all.x=T) %>% select(1,101:103,2:100)
# 
# addix <- noinfo2[complete.cases(noinfo2$totalSampledAreaShrubSapling),]
# 
# #### search for data from other years
# noinfo3 <- noinfo2[!(complete.cases(noinfo2$totalSampledAreaShrubSapling)),]
# noinfo3$plotID[!(noinfo3$plotID %in% combi$plotID)] ## only SCBI_013 missing, 2018 data
# 
# test <- shrubi %>% mutate(plotID = substr(shrubi$loc,1,8)) %>% select(105,1:104)
# testi <- test[test$plotID %in% noinfo3$plotID,]
# 
# testx <- shrub_perc[shrub_perc$plotID %in% noinfo3$plotID,]
# noinfo3$plotID %in% testx$plotID
# unique(testx$plotID)
# 
# ### nothing found, no matches between sample year SCBI_013 = 2018 and plot size info
# 
# ### add data
# addi_perc <- addix %>% mutate_at(vars(6:103),funs(. / totalSampledAreaShrubSapling*100))
# combi2 <- rbind(combi[,-c(104:106)], addi_perc)
# 
# ####
# ploti <- unique(substr(infi$plotyear, 1,8))
# 
# 
# ploti[!(ploti %in% combi2$plotID)]
# 
# 
# # write.csv(combi2,"./R_output/biodiversity_compiled/shrub_cover.csv", row.names = F)
# 
#   
#   
#   