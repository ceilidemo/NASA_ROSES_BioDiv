library(vegan)
library(tidyverse)

### Data from NEON unzipped
plants_save <- read.csv("./orig_data/NEON_presence-cover-plant/stackedFiles/div_1m2Data.csv")

### plant codes
XX <- read.csv("./R_input/NEON_plantcodes_lookup.csv")

### var trafo
plants <- div_1m2Data %>% 
  mutate(subID = paste(plotID, subplotID, sep="_"), 
         year = sapply(strsplit(as.character(endDate), "-"),"[",1),
         month = sapply(strsplit(as.character(endDate), "-"),"[",2), 
         year_month = paste(year, month, sep="_"), 
         site_dat = paste(siteID, year_month, sep="_"), 
         site_year = paste(siteID, year, sep="_"),
         taxon_plus = paste(taxonID, morphospeciesID, sep="_"))
         
### replace empty taxonIDs
plants <- plants %>%
  mutate(taxonID = ifelse(taxonID == "" | is.na(taxonID), otherVariables, taxonID))
plants$taxon_plus[plants$taxon_plus=="_"] <- plants$otherVariables[plants$taxon_plus=="_"]

### check duplicates
plants <- plants %>% 
  mutate(dupstring = paste(subID, year, taxonID, sep="_"),
         dupstring2 = paste(subID, endDate, taxonID, sep="_"), 
         dupstring3 = paste(subID, year, taxon_plus, sep="_"))

#####################
#### TRY transpose
sy <- unique(plants$site_year)
colsin <- c("taxonID","subID","percentCover")

trans_funct <- function(i){
  spread(data = plants[plants$site_year==sy[i],colsin], key = taxonID,  
         value = percentCover)
}

poss_trans = possibly(trans_funct, otherwise = "some duplicates")

### Apply function
res <- map(1:length(sy), poss_trans)

### Save results for data w/o duplicates
res_sel <- Filter(function(x) is.data.frame(x) == T, res)

nam <-  sy[which(!(res=="some duplicates"))]

for (i in 1:length(res_sel)){
  res_sel[[i]]$subID <- paste0(nam[i], substr(res_sel[[i]]$subID,nchar(res_sel[[i]]$subID)-10, 
                                              nchar(res_sel[[i]]$subID)))
}

for (i in 1:length(res_sel)){
  write.csv(res_sel[[i]], file=paste0("data_out/plantdiv/",nam[i], "_plantdiv.csv"), row.names = F)
}


##################################################
### Inventories with duplicates
dupls <- sy[which(res=="some duplicates")]
dupls <- sort(dupls)

##### SEL datasets with duplicates ### site_year
# vars_to_keep <- c("uid","namedLocation","taxonID","endDate","subID","otherVariables","year_month","dupstring",
#                   "dupstring2","percentCover","identificationQualifier","boutNumber", "eventID", "samplingProtocolVersion")

##################
res2 <- list()

for(a in 1:length(dupls)){
  sel <- plants[plants$site_year==dupls[a],] 
  dupcases <- sel[sel$dupstring %in% sel$dupstring[duplicated(sel$dupstring)],] #### check end date
  unis <- sel[!(sel$dupstring %in% dupcases$dupstring),]

  #### rm NAs for %cover
  dupcases <- dupcases %>% drop_na(percentCover)
  unis <- unis %>% drop_na(percentCover)
  dupcases <- dupcases[order(dupcases$dupstring2),]
  unis <- unis[order(unis$dupstring2),]
  
  if (nrow(dupcases)>0){
    write.csv(dupcases,paste0("data_out/plant_div_dupl_entries/", dupls[a], "_dupl_entries.csv"), row.names = F)
    
    y <- split(dupcases,dupcases$dupstring3)
    
    ### remove true duplicates, everything but uid, person names are the same => keep only one
    z <- list()
  
    for(i in 1:length(y)){
      if (any(duplicated(y[[i]][,!(colnames(y[[i]]) %in% c("uid","remarks", "measuredBy","recordedBy"))]))==TRUE){
              z[[i]] <- y[[i]][1,]
            }else{
              z[[i]] <- y[[i]]
            }
    }
    
    yy <- do.call(rbind,z)
    
    ### Test true duplicates
    # pp <- dupcases[!(dupcases$uid %in% yy$uid),]
    # p <- dupcases[dupcases$dupstring3 %in% pp$dupstring3,]
    # table(p[1,-1]==p[2,-1])
      
    ### average plants per plot sample on multiple dates
    indi <- yy$dupstring3[duplicated(yy$dupstring3)]
    dups_multidate <- yy[yy$dupstring3 %in% indi,]
    
    dups_d <- dups_multidate %>% 
      mutate(taxID = paste0(taxonID, identificationQualifier)) %>% 
      group_by(taxID,subID,taxonID, dupstring3, taxon_plus) %>%
      summarise(percentCover = mean(percentCover))
    
    if(length(indi)>0){
      yy <-  yy[!(yy$dupstring3 %in% indi),] %>% 
        mutate(taxID = paste0(taxonID, identificationQualifier))%>% 
        select(taxID,subID,taxonID, dupstring3,percentCover) 
    } else{
      yy <-  yy %>% 
        mutate(taxID = paste0(taxonID, identificationQualifier))%>% 
        select(taxID,subID,taxonID, dupstring3,percentCover) 
    }

    yy <- bind_rows(yy,dups_d)
      
      ### Update: Dave Barnett: remove dupls
      ### 2 ways to deal with duplicates: sum and average
      
      taxID <- unique(yy$taxID)
  
      ##### for all taxonIDs: when there is an entry for taxID with cf ---> use sum 
      look <- as.data.frame(matrix(data = taxID,ncol=1, dimnames = list(NULL,"taxID")))
      look$what <- "sum"
      look$taxID <- as.character(look$taxID)
      look$taxonID <- sapply(strsplit(look$taxID,"cf"),"[",1)
      look <- look[order(look$taxID),]
      
      # look <- look %>%
      #   mutate(what=replace(what, grepl("^2", look$taxID),"mean"))
      
      yy$taxID <- gsub("cf. species", "",yy$taxID) 
      
      ## split cases
      dups_tosum <- yy[yy$taxID %in% look[look$what=="sum",]$taxID,]
      dups_toaverage <- yy[yy$taxID %in% look[look$what=="mean",]$taxID,] 
      
      dups_s <- dups_tosum %>% 
        group_by(taxID,subID, taxonID) %>%
        summarise(percentCover =sum(percentCover))
      
      dups_m <- dups_toaverage %>% 
        group_by(taxon_plus,subID, taxonID) %>%
        summarise(percentCover = mean(percentCover))%>%
        rename(taxID=taxon_plus)
      
      dups_sm <- rbind(dups_s, dups_m)
      
      dupsi <- dups_sm %>% 
        group_by(taxID,subID) %>%
        summarise(percentCover = sum(percentCover))%>%
        rename(taxonID=taxID)

    #### combine 
    unisx <- unis %>% 
      mutate(taxonID = ifelse(grepl("^2",taxonID),taxon_plus,taxonID)) %>% 
      select("taxonID", "subID", "percentCover") 
      
    sel_corr <- bind_rows(dupsi, unisx)
    
    sel_corr$taxonID[grepl("_$",sel_corr$taxonID)] <- gsub("_","",  sel_corr$taxonID[grepl("_$",sel_corr$taxonID)])
    
    dupi <- paste(sel_corr$subID, sel_corr$taxonID, sep="_")
    anyDuplicated(dupi)
    
    ##### TRANSPOSE AGAIN
    sel_corri <- as.data.frame(sel_corr)
    
    ########### spread complains about dups
    res2[[a]] <- spread(data = sel_corri, key = taxonID, value = percentCover)
    res2[[a]]$subID <- paste0(dupls[a], substr(res2[[a]]$subID,5,length(res2[[a]]$subID)))
    
  }else{
    sels <- unis %>% select("taxonID", "subID", "percentCover")
    
    res2[[a]] <- spread(data = sels, key = taxonID, value = percentCover)
    res2[[a]]$subID <- paste0(dupls[a], substr(res2[[a]]$subID,5,length(res2[[a]]$subID)))
  }
}

for (i in 1:length(res2)){
  write.csv(res2[[i]], file=paste0("data_out/plantdiv/",dupls[i], "_plantdiv.csv"), row.names = F)
}

####### END ##########


