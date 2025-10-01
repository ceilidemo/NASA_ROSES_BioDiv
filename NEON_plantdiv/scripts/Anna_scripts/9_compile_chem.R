library(tidyverse)
library(neonUtilities)

# stackByTable("/Volumes/NEON_2021/NEON_SD/orig_data/NEON_traits-foliar.zip")

####
field_orig <- read.csv("./orig_data/NEON_traits-foliar/stackedFiles/cfc_fieldData.csv")
chem_orig <- read.csv("./orig_data/NEON_traits-foliar/stackedFiles/cfc_chemistrySubsampling.csv") ### other sampling IDs
CN_orig <- read.csv("./orig_data/NEON_traits-foliar/stackedFiles/cfc_carbonNitrogen.csv")
CHL_orig <- read.csv("./orig_data/NEON_traits-foliar/stackedFiles/cfc_chlorophyll.csv")
LIG_orig <- read.csv("./orig_data/NEON_traits-foliar/stackedFiles/cfc_lignin.csv")
LMA_orig <- read.csv("./orig_data/NEON_traits-foliar/stackedFiles/cfc_LMA.csv")
elem_orig <- read.csv("./orig_data/NEON_traits-foliar/stackedFiles/cfc_elements.csv")

##################################################
### Filter data based on QF and combine repeats ##

unique(CN_orig$dataQF)
unique(CN$siteID) ### CN drying protocol errors everywhere except GUAN and PUUM

table(CN_orig$percentAccuracyQF)


# CN w/o isotopes
CN <- CN_orig %>%
  filter(cnPercentQF=="OK", percentAccuracyQF =="OK", !(remarks=="N amplitude below threshold for reliable data.")) %>% 
  select(c(4:7,9,11,16:18)) %>% mutate(collectDate=substr(collectDate,1, 10)) %>% 
  group_by_at(vars(-nitrogenPercent, -carbonPercent, -CNratio,-cnSampleID)) %>% 
  summarise(N_perc_mean = mean(nitrogenPercent, na.rm=T), 
            C_perc_mean = mean(carbonPercent, na.rm = T), 
            CN_ratio_mean = mean(CNratio, na.rm = T))%>% 
  ungroup() %>%
  mutate_if(is.factor,list(as.character))


anyDuplicated(CN_orig$sampleID)
dupsi <- as.character(CN_orig[which(duplicated(CN_orig$sampleID)),]$sampleID)
dup <- CN_orig[CN_orig$sampleID %in% dupsi,]


# CN Isotopes
table(CN_orig$isotopeAccuracyQF)

CNiso <- CN_orig %>%
  filter(cnIsotopeQF=="OK", isotopeAccuracyQF =="OK", !(remarks=="d15N values out of range for duplicate.")) %>% 
  select(c(4:7,9,11,14,15)) %>% mutate(collectDate=substr(collectDate,1, 10)) %>% 
  group_by_at(vars(-d15N, -d13C,-cnSampleID)) %>% 
  summarise(d15N_mean = mean(d15N, na.rm=T), 
            d13C_mean = mean(d13C, na.rm = T))%>% 
  ungroup() %>%
  mutate_if(is.factor,list(as.character))

anyDuplicated(CNiso$sampleID)


###
unique(CHL_orig$dataQF)
table(CHL_orig$chlCarotWavelength5)


CHL <- CHL_orig %>% 
  filter(!(dataQF %in%"deprecatedMethod"),handlingQF =="OK", measurementQF =="OK",
                           !(sampleCondition =="brown"), release=="RELEASE-2021") %>%
  filter(chlorophyllSampleID!="P3kunhL9imCoWAM/ZrwpAiYpU+Sgsb4xDZfeEcg13YyFh0IOSpGEmQ==")%>% ### two meas vary widely
  select(c(4:7,9,11,15,16,27:29)) %>% mutate(collectDate=substr(collectDate,1, 10)) %>%
  # group_by_at(vars(-extractChlAConc,-extractChlBConc, -extractCarotConc,-sampleID)) %>% 
  # summarise(ChlAConc_mean = mean(extractChlAConc, na.rm=T), 
  #           ChlBConc_mean = mean(extractChlBConc, na.rm=T), 
  #           CarotConc_mean = mean(extractCarotConc, na.rm=T))%>%
  mutate_if(is.factor,list(as.character))%>%
  rename("PigmentFreshMass"="freshMass")
  
names(CHL)
anyDuplicated(CHL$chlorophyllSampleID)

dupsi <- as.character(CHL[which(duplicated(CHL$chlorophyllSampleID)),]$chlorophyllSampleID)
dup <- CHL_orig[CHL_orig$chlorophyllSampleID %in% dupsi,]
CHL_orig[CHL_orig$chlorophyllSampleID %in% dupsi,]$uid

xx <- CHL_orig[CHL_orig$chlorophyllSampleID =="P3kunhL9imCoWAM/ZrwpAiYpU+Sgsb4xDZfeEcg13YyFh0IOSpGEmQ==",]
write.csv(dup, "./docu/CHL_duplicates_difference.csv",row.names = F)

dup_ex <- dup[dup$sampleID %in% dupsi[11],]

###
table(LIG_orig$measurementQF)

LIG <- LIG_orig %>% 
  filter(measurementQF=="OK", !(accuracyQF =="run QA criteria not met"))%>% 
  select(c(4:7,9,11,15:16))%>% mutate(collectDate=substr(collectDate,1, 10))%>% 
  group_by_at(vars(-ligninPercent, -cellulosePercent,-ligninSampleID)) %>% 
  summarise(Lignin_perc_mean = mean(ligninPercent, na.rm = T), 
            Cell_perc_mean = mean(cellulosePercent, na.rm=T)) %>%  
  ungroup() %>%
  mutate_if(is.factor,list(as.character))

anyDuplicated(LIG$sampleID)

###
unique(LMA_orig$dataQF)

LMA <- LMA_orig %>% 
  filter(lmaSampleCondition=="OK", percentGreen>=90,release=="RELEASE-2021") %>%
  mutate_if(is.factor,list(as.character)) %>%
  mutate(plotID = sapply(strsplit(namedLocation,"\\."),"[",1)) %>%
  select(c(4,27,5:6,8,9,14,17,21,22,23,28,25,29))%>% mutate(collectDate=substr(collectDate,1, 10))%>% 
  group_by_at(vars(-lmaSampleID, -freshMass, -dryMass, -percentGreen, 
  -leafMassPerArea, -dryMassFraction)) %>% 
  summarise(freshMass = mean(freshMass,na.rm=T), 
            dryMass = mean(dryMass,na.rm=T),
            percentGreen = mean(percentGreen,na.rm=T), 
            leafMassPerArea = mean(leafMassPerArea,na.rm=T),
            dryMassFraction = mean(dryMassFraction,na.rm=T)) %>%  
  ungroup() 

anyDuplicated(LMA$sampleID)
dupsi <- LMA[which(duplicated(LMA$sampleID)),]$sampleID
dup <- LMA[LMA$sampleID %in% dupsi,]


#####################
### Combine data ####
field <- field_orig %>% mutate_if(is.factor,list(as.character)) 
chem <- chem_orig %>% mutate_if(is.factor,list(as.character))

# dupsi <- as.character(LMA[which(duplicated(LMA$sampleID)),]$sampleID)
# dup <- LMA_orig[LMA_orig$sampleID %in% dupsi,]
# names(LMA_orig)[grepl("QF",names(LMA_orig))]

# test <- read.csv("//Volumes/Backup Plus/Anna/NEON/R_input/plants/stackedFiles_chem_Feb2020/cfc_fieldData.csv")
# names(test)[c(5,15:17,21:25,28:32,37)]

names(field)

dat <- field %>% 
  select (plotID,collectDate,subplotID,sampleType, tagID,individualID, taxonID,
          scientificName, plantStatus, clipLength, clipWidth, percentCoverClip, 
          toxicodendronPossible,sampleID, chlorophyllSampleID)%>%
  left_join(CN[,-5], by = c("plotID","sampleID"))%>%
  left_join(CNiso[,-c(1,3,5)], by = c("plotID","sampleID"))%>%
  left_join(CHL[,-c(1,3,6)], by = c("plotID","sampleID","chlorophyllSampleID")) %>%
  left_join(LIG[,-c(1,3,5)], by = c("plotID","sampleID")) %>%
  left_join(LMA[, -c(1,2,3,5,6)], by = c("plotID","sampleID"))

anyNA(dat2$siteID)

dat2 <- dat %>% mutate_if(is.factor,list(as.character)) %>%
  mutate(siteID = sapply(strsplit(plotID,"_"),"[",1))%>%
  select(1:29,31:35,30)

names(dat2)

### Calculate pigment vals 
dat2 <- dat2 %>% 
  mutate(ChlA_mg_g=extractChlAConc*solventVolume/PigmentFreshMass/dryMassFraction*0.001)%>% 
  mutate(ChlB_mg_g=extractChlBConc*solventVolume/PigmentFreshMass/dryMassFraction*0.001)%>% 
  mutate(Car_mg_g=extractCarotConc*solventVolume/PigmentFreshMass/dryMassFraction*0.001)%>% 
  mutate(ChlA_mg_m2=ChlA_mg_g*leafMassPerArea)%>% 
  mutate(ChlB_mg_m2=ChlB_mg_g*leafMassPerArea)%>%   
  mutate(Car_mg_m2=Car_mg_g*leafMassPerArea)

write.csv(dat2, "./R_output/chem/chemistry_combined.csv", row.names = F)

### compare to old chem data from FEB 2020

test <- read.csv("~/Documents/CABO/NEON/R_input/plants/chemistry_combined_May2020.csv")
names(test)
names(dat2)


 