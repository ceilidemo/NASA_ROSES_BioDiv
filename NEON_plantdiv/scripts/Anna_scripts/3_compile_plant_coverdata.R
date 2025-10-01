##################################################
##### KEEP ALL NOT PLANT COVER, some can be removed later
library(tidyverse)
library(neonPlantEcology)



look <- read.csv("./R_input/NEON_plantcodes_lookup.csv",as.is = T)
# plant_locs <- read.csv("./R_input/plant_locs.csv") ### location data

fls <- list.files("data_out/plantdiv/", full.names = T)

# test <- read.csv(fls[grepl("HARV_2018", fls)])
# testi <- test %>% mutate(plotID=substr(subID,1,13), .before=1) %>%
#   group_by(plotID) %>% summarise_at(-c(1:2),mean, na.rm=T)
# 
# rowSums(testi[,-1], na.rm=T) 


plants_all <- map(fls,read.csv, check.names=F)

# out <- as.character(look$taxonID[1:12]) ### variables to exclude, can be done later

zz <- list()
for (i in 1:length(plants_all)) {
  
  plot_dat <- plants_all[[i]] %>%
    mutate(plotID = substr(subID, 1, 13)) %>%
    select(plotID, everything())
  
  ### aggregate by plot
  plot_dat3 <- plot_dat %>%
    group_by(plotID) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  ### check names
  names(plot_dat3)[!(names(plot_dat3) %in% look$taxonID)]
  
  zz[[i]] <- plot_dat3
}


df <- bind_rows(zz)
df[is.na(df)] <- 0

### some tests ####
# hist(rowSums(df[,-1],na.rm=T), breaks = 20)

# testx <- df %>% filter(df$plotID=="ABBY_2018_006") %>% 
#   select(-moss) %>%
#   select_if(~sum(!is.na(.)) > 0)
# rowSums(testx[,-1], na.rm = T)

# testi <- plants_all[[1]]
# hist(rowSums(testi[,-1], na.rm = T), breaks = 20)


######## Rename df columns based on lookup table
miss <- sort(names(df)[!(names(df)%in%look$taxonID)])
miss <- miss[!(grepl("2PLA", miss))]

#### .cf family
miss[grepl("family", miss)]
look[grepl("Poacea", look$species),]

dfx <- df
colnames(dfx)[grepl("cf. fam",colnames(dfx))] <- gsub("cf. family", "",colnames(dfx)[grepl("cf. fam",colnames(dfx))])

#### .cf genus
miss[grepl("genus", miss)]
look[grepl("EUPAT", look$taxonID),]

colnames(dfx)[grepl("cf. genus",colnames(dfx))] <- gsub("cf. genus", "",colnames(dfx)[grepl("cf. genus",colnames(dfx))])

missx <- sort(names(dfx)[!(names(dfx)%in%look$taxonID)])
missx[!(grepl("2PLA", missx))]

### other cf
look[grepl("OPEN3", look$taxonID),]
look[grepl("Opuntia spp", look$species),]

colnames(dfx)[grepl("OPEN3cf. subspecies",colnames(dfx))] <- gsub("OPEN3cf. subspecies", "OPUNTSPP",colnames(dfx)[grepl("OPEN3cf. subspecies",colnames(dfx))])

look[grepl("ARTET", look$taxonID),]
look[grepl("Aristida ternipes", look$species),]

colnames(dfx)[grepl("ARTETcf. variety",colnames(dfx))] <- gsub("ARTETcf. variety", "ARTE3",colnames(dfx)[grepl("ARTETcf. variety",colnames(dfx))])

### aff. taxa -> keep info from lower level
look[grepl("CABI5", look$taxonID),]
look[grepl("Carex sp.", look$species),]

colnames(dfx)[grepl("CABI5aff. species",colnames(dfx))] <- gsub("CABI5aff. species", "CAREX",colnames(dfx)[grepl("CABI5aff. species",colnames(dfx))])
colnames(dfx)[grepl("CAREXaff. species",colnames(dfx))] <- gsub("CAREXaff. species", "CAREX",colnames(dfx)[grepl("CAREXaff. species",colnames(dfx))])

look[grepl("ELEOC", look$taxonID),]
look[grepl("Cypera", look$species),]

### keep family
colnames(dfx)[grepl("ELEOCaff. genus",colnames(dfx))] <- gsub("ELEOCaff. genus", "CYPERA",colnames(dfx)[grepl("ELEOCaff. genus",colnames(dfx))])

### keep genus
look[grepl("Mimosa", look$species),]
colnames(dfx)[grepl("MIDYaff. species",colnames(dfx))] <- gsub("MIDYaff. species", "MIMOS",colnames(dfx)[grepl("MIDYaff. species",colnames(dfx))])

#### deal with spp's 
colnames(dfx)[grepl("SPP$",colnames(dfx))] <- gsub("SPP$", "",colnames(dfx)[grepl("SPP$",colnames(dfx))])

## test
missx <- sort(names(dfx)[!(names(dfx)%in%look$taxonID)])
missx[!(grepl("2PLA", missx))]

look[grepl(pattern = "Wyeth",look$scientificName),]
colnames(dfx)[grepl("OPLON",colnames(dfx))] <- gsub("OPLON", "OPLONSPP",colnames(dfx)[grepl("OPLON",colnames(dfx))])
colnames(dfx)[grepl("SHEPH",colnames(dfx))] <- gsub("SHEPH", "SHEPHSPP",colnames(dfx)[grepl("SHEPH",colnames(dfx))])
colnames(dfx)[grepl("WYETH",colnames(dfx))] <- gsub("WYETH", "WYETHSPP",colnames(dfx)[grepl("WYETH",colnames(dfx))])

### deal with "/"
look[grepl("/",look$taxonID),]

colnames(dfx)[grepl("/",colnames(dfx))] <- sapply(strsplit(colnames(dfx)[grepl("/",colnames(dfx))],"/"),"[",1)
colnames(dfx)[grepl("CAEL3",colnames(dfx))] <- gsub("CAEL3", "KOMY",colnames(dfx)[grepl("CAEL3",colnames(dfx))])

look[grepl(pattern = "KOMY",look$taxonID),]

#### Test

missx <- sort(names(dfx)[!(names(dfx)%in%look$taxonID)])
missx[!(grepl("2PLA", missx))]


##########################
###### Summarize duplicates
dubs <- sort(colnames(dfx)[which(duplicated(colnames(dfx)))])

dubs_uni <- unique(dubs)
patt <- paste(dubs_uni, collapse = "$|")
patt <- paste0(patt,"$")

# dfxx <- dfx 

df_keep <- dfx[, !(grepl(patt,colnames(dfx)))]
anyDuplicated(colnames(df_keep))

df_change <- dfx[, grepl(patt,colnames(dfx))]
colnames(df_change) <- make.unique(colnames(df_change))
df_change <- df_change %>% select(sort(colnames(df_change)))


# sort(colnames(dfxx)[grepl(patt,colnames(dfxx))])
# colnames(dfxx)[grepl(patt,colnames(dfxx))] <- paste0(colnames(dfxx)[grepl(patt,colnames(dfxx))], "_X")

# patt <- paste(dubs_uni, collapse = "_X|")
# patt <- paste0(patt,"_X")
# 
# colnames(dfxx)[grepl(patt,colnames(dfxx))] <- make.unique(colnames(dfxx)[grepl(patt,colnames(dfxx))])

# sort(colnames(dfxx)[grepl(patt,colnames(dfxx))])
# sort(colnames(dfxx)[grepl("POA",colnames(dfxx))])

####
dubsi <- list()
for (i in 1:length(dubs_uni)){
  dubsi[[i]] <- colnames(df_change)[grepl(dubs_uni[i],colnames(df_change))]
} ### look at dubsi

uu <- unlist(dubsi)


### check dublicates with more than one entry
dubsi[lapply(dubsi,length)>2]
indi <- which(lapply(dubsi, length)>2)

### modify by hand
dubsi[indi[23]]

dubsi[[25]] <- dubsi[[25]][1:2]
dubsi[[62]] <- dubsi[[62]][1:3]
dubsi[[135]] <- dubsi[[135]][1:2]
dubsi[[141]] <- dubsi[[141]][1:2]
dubsi[[160]] <- dubsi[[160]][1:2]
dubsi[[172]] <- dubsi[[172]][1:2]
dubsi[[178]] <- dubsi[[178]][1:3]
dubsi[[189]] <- dubsi[[189]][1:2]

####
look[grepl("VIOLA",look$taxonID),]
dubsi[1:190]


## trying to automate the above part 
base_names <- str_replace(colnames(df_change), "\\..*", "")
dup_groups <- split(colnames(df_change), base_names)

sumi <- lapply(dup_groups, function(cols) {
  # subset df_change to those columns
  xx <- df_change %>% select(all_of(cols))
  # sum across them
  summed <- rowSums(xx, na.rm = TRUE)
  data.frame(summed)
})

# name the columns properly
for (i in seq_along(sumi)) {
  colnames(sumi[[i]]) <- names(dup_groups)[i]
}

df_sumi <- bind_cols(sumi)

df_all <- bind_cols(df_keep, df_sumi) %>%
  select(1, sort(colnames(.)[-1]))

#############
sumi <- list()

for(i in 1:length(dubsi)){
  xx <- df_change %>% select(dubsi[i][[1]])
  sumi[[i]] <- as.data.frame(rowSums(xx))
  colnames(sumi[[i]]) <- dubsi[i][[1]][1]
}

df_sumi <- do.call(bind_cols, sumi)

### test 
as.data.frame(df_change[1:30,grepl("VIOLA", colnames(df_change))])


############ DEL ########
#### loop through duplicates
# df2 <- dfxx
# 
# # # ### test
# # df2 <- as.data.frame(matrix(data=sample(1:60),nrow=10,
# #               dimnames = list(c(NULL),c("POACEA_X","POACEA_X.1","POAL2_X","POAR2_X","POARA2" ,"POA_X"))))
# # 
# # dubs_uni <- dubs_uni[grepl("POA",dubs_uni)]
# # dubs_uni <- dubs_uni[-1]
# 
# for(i in 1:length(dubs_uni)){
#   df2[, ncol(df2)+1] <- rowSums(df2[,grepl(dubs_uni[i], colnames(df2))])
#   colnames(df2)[ncol(df2)] <- dubs_uni[i]
# }
# 
# df3 <- df2 %>% select(-names(df2)[names(df2)%in%uu])
# df3 <- df3 %>% select(c(1, sort(colnames(df3)[-c(1)])))


###### combine 
df_all <- bind_cols(df_keep,df_sumi)
df_all <- df_all %>% select(1,sort(colnames(df_all)[-1]))


anyDuplicated(colnames(df_all))
colnames(df_all)[!(colnames(df_all)%in% look$taxonID)]


write.csv(df_all, "data_out/biodiv_compiled/herb_cover_allplots.csv", row.names = F)

#######################
#### GENERATE TREE
library(V.PhyloMaker)

species_list <- data.frame(
  species = colnames(df),   
  genus = gsub("_.*", "", colnames(df)),   
  family = NA               
)

tre <- phylo.maker(
  sp.list = species_list,
  tree = GBOTB.extended,  
  nodes = nodes.info.1,   
  scenarios = "S3"        
)

# Save tree object
saveRDS(tre, "data_out/phylo/allplants_scenario3_tree.rds")


#######################
#### PHYLO MATCHED
df_save <- read.csv("data_out/biodiv_compiled/herb_cover_allplots.csv", check.names = F)  
anyDuplicated(colnames(df))

tre <- readRDS("./R_output/phylo/allplants_scenario3_tree.rds") 
phy <- tre$scenario.3$tip.label

checked <- read.csv("./R_output/phylo/plants_checked.csv")
checked$species_nam <- gsub(" ", "_", checked$species_nam)
checked$species <- gsub(" ", "_", checked$species)

table(colnames(df) %in% checked$species_nam)
colnames(df)[!(colnames(df)%in% checked$species_nam)]

# colnames(df3)[grepl("Chenopodiaceae_sp.", colnames(df3))] <- "Chenopodium_sp." 
# colnames(df3)[grepl("Asclepi", colnames(df3))]

# tnrs <- read.csv("R_output/phylo/tnrs_results_missing_sp_bestresult_small.csv", row.names = 1)
# tnrs <- tnrs %>% mutate_if(is.factor, as.character)

# tnrs$Name_submitted <- gsub(" ", "_",tnrs$Name_submitted)
# tnrs$Name_matched <- gsub(" ", "_",tnrs$Name_matched)
# 
# tnrs2 <- read.csv("R_output/phylo/plants_checked_400m_edits_relatives.csv")
# tnrs2 <- tnrs2 %>% mutate_if(is.factor, as.character)
# 
# tnrs2$species <- gsub(" ", "_",tnrs2$species)
# tnrs2$species_nam <- gsub(" ", "_",tnrs2$species_nam)

table(checked$species %in% phy)

### Replacements 
mat <- as.data.frame(colnames(df)[-1])
names(mat) <- "species_nam"

dubs <- checked$species_nam[which(duplicated(checked$species_nam))]
dubsi <- checked[checked$species_nam %in% dubs,]
checked2 <- checked[-which(duplicated(checked$species_nam)),]

mat <- mat %>% left_join(checked2,all.x=T,sort=F)
table(mat$species_nam==colnames(df)[-1])

table(mat$species %in% phy)

###
sort(mat$species[!(mat$species %in% phy)]) #### known as missing species: ok.

xx <- mat[is.na(mat$species),] ### NAs

xx$species_nam
phy[grepl("Lycopo",phy)] 

mat[grepl("Chenopodiaceae_sp.",mat$species_nam),]$species <- "Chenopodiaceae_xx" 
mat[grepl("Asclepiadaceae_sp.",mat$species_nam),]$species <- "Asclepias_xx" 
mat[grepl("Lycopodium_hickeyi",mat$species_nam),]$species <- "Lycopodium_hickeyi" 

mat[is.na(mat$species),]$species <- mat[is.na(mat$species),]$species_nam

###### rename
colnames(df)[-1] <- mat$species
out <- colnames(df)[-1][!(colnames(df)[-1] %in% phy)]

out2 <- out[!(grepl("_NA|Unknown",out))]


### remove duplicates
dubs <- mat$species[duplicated(mat$species)]
mat$species[mat$species %in% dubs] <- paste0(mat$species[mat$species %in% dubs],"_X")

df1 <- df
colnames(df1)[-1] <- mat$species

df2 <- df1[,-c(which(colnames(df1)%in% out2))]
sort(colnames(df2)[-1][!(colnames(df2)[-1] %in% phy)])

####
dubs_uni <- unique(dubs)

####
dubsi <- list()
for (i in 1:length(dubs_uni)){
  dubsi[[i]] <- colnames(df2)[grepl(dubs_uni[i],colnames(df2))]
} ### look at dubsi

uu <- unlist(dubsi)

#### loop through duplicates
df3 <- df2
head(df3[,grepl(dubs_uni[i], colnames(df3))])

for(i in 1:length(dubs_uni)){
  df3[, ncol(df3)+1] <- rowSums(df3[,grepl(dubs_uni[i], colnames(df3))])
  colnames(df3)[ncol(df3)] <- dubs_uni[i]
}

names(df3)[grepl(uu[1], names(df3))]

df4 <- df3 %>% select(-names(df3)[names(df3)%in%uu])
df4 <- df4 %>% select(c(1, sort(colnames(df4)[-c(1)])))

anyDuplicated(colnames(df4))
write.csv(df4, "./R_output/biodiversity_compiled/herb_cover_allplots_phylomatch.csv", row.names = F)  



