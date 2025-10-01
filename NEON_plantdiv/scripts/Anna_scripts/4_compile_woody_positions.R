# library(devtools)
# install_github('NEONScience/NEON-geolocation/geoNEON', dependencies=TRUE)
library(geoNEON)
library(sp)
library(rgdal)
library(rgeos)
library(tidyverse)
library(viridis) 
library(neonUtilities)
source("./R_scripts/NEON_token.R")

### Unpack woody data zips
# stackByTable("./orig_data/NEON_struct-woody-plant.zip")

#### BELOW test if only one tree with data for JERC: JERC_054_2019

maptag <- read.csv("./orig_data/NEON_struct-woody-plant/stackedFiles/vst_mappingandtagging.csv") ### location and plant ID
indiv <- read.csv("./orig_data/NEON_struct-woody-plant/stackedFiles/vst_apparentindividual.csv") 

indiv2 <- indiv %>%
  mutate(year= sapply(strsplit(date, "-"),"[",1))%>% 
  mutate(plotyear = paste(plotID,year, sep="_"))

# tree_coords <- getLocTOS(maptag, "vst_mappingandtagging", token = NEON_TOKEN)
# write.csv(tree_coords, "./R_output/woody/tree_coords/getLocTOC_trees_token.csv", row.names = F)

# tab <- warnings()
# tabb <- as.data.frame(sapply(strsplit(names(tab),"was "),"[",2))
# colnames(tabb) <- "miss_from_getLocTOC"
# write.csv(tabb,"./R_output/woody/missing_coords.csv", row.names = F)

# tree_coords <- read.csv("./R_output/woody/getLocTOC_trees.csv")
# tree_coords2 <- tree_coords[,-c(2:5,7,11:13,17,21,26:32,36,40)]

# ### Merge with growth measurements
table(indiv$individualID %in% tree_coords$individualID)
# table(indiv$uid %in% tree_coords2$uidD)
xx <- indiv[!(indiv$individualID %in% tree_coords$individualID),]
# 
# unique(xx$date)
write.csv(xx,"./R_output/woody/tree_coords/indIDs_missing_coords_token.csv", row.names = F)

xx_old <- read.csv("./R_output/woody/tree_coords/indIDs_missing_coords.csv")
table(xx$uid %in% xx_old$uid) ### are the same

######### MORE ENTRIES IN TREE COORDS3 after using token
tree_coordsx <- indiv2 %>% left_join(tree_coords,by = "individualID", all.x=T)

tree_coords3 <- tree_coordsx %>% 
  mutate(maxCrownDiam_mean= rowMeans(cbind(ninetyCrownDiameter,maxCrownDiameter), na.rm=T))%>%
  mutate(indiv = paste(sapply(strsplit(individualID,"\\."),"[",5), taxonID, year,sep = "_")) %>% 
  drop_na(adjEasting,height, maxCrownDiam_mean) %>% arrange(individualID)
  # select(-colnames(tree_coords3)[grepl("\\.y",colnames(tree_coords3))])

colnames(tree_coords3)[grepl("\\.y",colnames(tree_coords3))] <- gsub("\\.y","_getloc",colnames(tree_coords3)[grepl("\\.y",colnames(tree_coords3))])
colnames(tree_coords3)[grepl("\\.x",colnames(tree_coords3))] <- gsub("\\.x","",colnames(tree_coords3)[grepl("\\.x",colnames(tree_coords3))])
write.csv(tree_coords3,"./R_output/woody/tree_coords/tree_coordinates_token.csv", row.names = F)

old <- read.csv("./R_output/woody/tree_coords/tree_coordinates.csv")

##################
tree_coords_save <- read.csv("./R_output/woody/tree_coordinates.csv")

# dubs <- tree_coords_save$individualID[duplicated(tree_coords3$individualID)]
# dubsi <- tree_coords_save[tree_coords_save$individualID %in% dubs,]

tree_coords3 <- tree_coords_save %>%
  filter(!(plantStatus=="Dead, broken bole"|plantStatus=="Standing dead"))%>%
  select(-c(basalStemDiameter, initialGapMeasurementDate, initialBandStemDiameter,
  initialDendrometerGap, dendrometerHeight,dendrometerGap, dendrometerCondition,
  bandStemDiameter, basalStemDiameterMsrmntHeight))%>% 
  distinct_at(-1, .keep_all = T)

tree_coords3[which(is.na(tree_coords3$ninetyCrownDiameter)),]$ninetyCrownDiameter <- tree_coords3[which(is.na(tree_coords3$ninetyCrownDiameter)),]$maxCrownDiameter

dubs <- tree_coords3$individualID[duplicated(tree_coords3$individualID)]
dubsi <- tree_coords_save[tree_coords_save$individualID %in% dubs,]%>%arrange(individualID)


###### JERC test
tree_coords <- read.csv("./R_output/woody/tree_coords/getLocTOC_trees_token.csv")
tree_coords <- tree_coords %>% filter(siteID=="JERC") %>%
  filter(adjEasting>0)  %>% filter(adjNorthing>0)

indiv2 <- indiv2 %>% filter(siteID=="JERC")%>% filter(height>0)%>%
  filter(maxCrownDiameter>0)
  
maptag <- maptag %>% filter(siteID=="JERC")

# ### Merge with growth measurements
table(indiv2$individualID %in% tree_coords$individualID)
# table(indiv$uid %in% tree_coords2$uidD)
xx <- indiv2[!(indiv2$individualID %in% tree_coords$individualID),]

# write.csv(xx,"./R_output/woody/indIDs_missing_coords.csv", row.names = F)

tree_coordsx <- indiv2 %>% left_join(tree_coords,by = "individualID", all.x=T)

# write.csv(tree_coords, "./R_output/JERC_tree_coords.csv",row.names=F)
# write.csv(indiv2, "./R_output/JERC_indiv.csv",row.names=F)

tree_coordsx <- read.csv("./R_output/JERC_tree_coords.csv")
indiv2x <- read.csv("./R_output/JERC_indiv.csv")
table(indiv2$individualID %in% tree_coords$individualID)
# 
# tree_coords <- tree_coords %>% select(individualID, adjEasting, adjNorthing)
# indiv2 <- indiv2 %>% select(individualID,height, maxCrownDiameter, ninetyCrownDiameter)
# test <- indiv2 %>% inner_join(tree_coords,by = "individualID")


tree_coords3 <- tree_coordsx %>% 
  mutate(maxCrownDiam_mean= rowMeans(cbind(ninetyCrownDiameter,maxCrownDiameter), na.rm=T))%>%
  mutate(indiv = paste(sapply(strsplit(individualID,"\\."),"[",5), taxonID, year,sep = "_"))%>%  
  drop_na(adjEasting,height, maxCrownDiam_mean) %>% arrange(individualID)
# select(-colnames(tree_coords3)[grepl("\\.y",colnames(tree_coords3))])

colnames(tree_coords3)[grepl("\\.y",colnames(tree_coords3))] <- gsub("\\.y","_getloc",colnames(tree_coords3)[grepl("\\.y",colnames(tree_coords3))])
colnames(tree_coords3)[grepl("\\.x",colnames(tree_coords3))] <- gsub("\\.x","",colnames(tree_coords3)[grepl("\\.x",colnames(tree_coords3))])
# write.csv(tree_coords3,"./R_output/woody/tree_coordinates.csv", row.names = F)


######################
### PLOTSIZE 
plotsize <- read.csv("//Volumes/NEON_2021/NEON_SD/orig_data/NEON_struct-woody-plant/stackedFiles/vst_perplotperyear.csv")
# plotsize$year <- sapply(strsplit(as.character(plotsize$date), "-"),"[",1)
# plotsize$plotyear <- paste(plotsize$plotID, plotsize$year, sep="_")

centroids <- plotsize %>%
  select(c(5,6,13,14,7,31:34)) %>%
  distinct(easting, northing,.keep_all = TRUE) %>%
  arrange(plotID)

anyDuplicated(centroids$plotID)

### check if all plots with coordinates in centroids
table(unique(tree_coords3$plotID) %in% centroids$plotID)
miss <- sort(unique(tree_coords3$plotID)[!(unique(tree_coords3$plotID)%in%centroids$plotID)])


#### Plot coordinates: size is not right take this for missing WOOD plots
centr <- read.csv("./orig_data/All_NEON_TOS_Plots_V8/All_NEON_TOS_Plot_Centroids_V8.csv")

centr_WOOD <- centr %>% filter(subtype=="basePlot") %>%
  filter(plotID %in% miss) %>%
  mutate(totalSampledAreaTrees=400, year="2014")%>%
  mutate(plotyear= paste(plotID, year, sep = "_"))
  
centr_WOOD <- centr_WOOD %>% select(colnames(centr_WOOD)[colnames(centr_WOOD) %in% colnames(centroids)])

# colnames(centroids) [!(colnames(centroids) %in% colnames(centr_WOOD))]

centroids <- centroids %>% bind_rows(centr_WOOD) ### add missing plots

radius <- sqrt(centroids$totalSampledAreaTrees)/2 # radius in meters
yPlus <- centroids$northing+radius
xPlus <- centroids$easting+radius
yMinus <- centroids$northing-radius
xMinus <- centroids$easting-radius
# calculate polygon coordinates for each plot centroid. 
square=cbind(xMinus,yPlus,  # NW corner
             xPlus, yPlus,  # NE corner
             xPlus,yMinus,  # SE corner
             xMinus,yMinus, # SW corner
             xMinus,yPlus)  # NW corner again - close ploygon
# Extract the plot ID information
ID <- as.character(centroids$plotID)
# crs <- as.character(unique(centroids$crsutm))
# crs <- crs[!is.na(crs)]
polys <- SpatialPolygons(mapply(function(poly, id){
  xy <- matrix(poly, ncol=2, byrow=TRUE)
  Polygons(list(Polygon(xy)), ID=id)
}, 
split(square, row(square)), ID))
# proj4string=CRS(crs))
polydf <- SpatialPolygonsDataFrame(Sr = polys, data = centroids,match.ID = F)
saveRDS(polydf, "./R_input/Plot_polygons.rds")


###################################
#### MAKE PLOTS 
# ###################################
# loc <- sort(unique(tree_coords3$plotyear))
# 
# for(i in 1:length(loc)){
#   sel <- tree_coords3[tree_coords3$plotyear==loc[i],]
#   sel <- sel[order(sel$height,decreasing = F),] ### order sel by height
#   # selx <- sel[is.na(sel$maxBaseCrownDiameter),]
#   # sel <- sel[!(is.na(sel$maxCrownDiam_mean)),]
#   
#   if(!(anyDuplicated(sel$individualID)==0)){
#     sel <- sel[-which(duplicated(sel$individualID)),] 
#   }
#   
#   if(nrow(sel)>0){
#     coordinates(sel) <-c("adjEasting","adjNorthing" ) ### X before y
#     
#     ### ellipse area and radius for circle of same area
#     elli_a <- ((sel@data$maxCrownDiameter)/2)*((sel@data$ninetyCrownDiameter)/2)*pi
#     elli_r <- sqrt(elli_a/pi)
#     
#     buff <- gBuffer(sel,byid = T,width = elli_r, quadsegs=100, id = sel@data$individualID) 
#     
#     # buff <- gBuffer(sel,byid = T,width = sel@data$maxCrownDiam_mean/2) ### buffer by crown radius
#     buff@plotOrder <- order(buff@data$height)
#     # select plotID
#     ploti <- subset(polydf,plotID==substr(loc[i],1,8))
#     
#     # Plot
#     pdf(file = paste0("./R_output/figures/woody_diameter_byheight_allyears/", loc[i],"_diam_perheight.pdf"),width = 9)
#     plot(ploti,border="green3", lwd=3)
#     # plot(buff, col=magma(nrow(sel),direction = -1, alpha=0.8,begin = 0.15),border="red",lwd=100,add=T)
#     plot(buff, col=magma(nrow(sel),direction = -1, alpha=0.8,begin = 0.15),add=T)
#     # plot(buff[buff$individualID=="LIST2_11070",], col=2, add=T)
#     # plot(buff[buff$individualID=="QUVE_11079",], col=3, add=T)
#     plot(sel, add=T)
#     text(x = sel@coords[,1], y=sel@coords[,2],pos = 2,cex=1.2,
#          labels = substr(sel@data$indiv,1,nchar(sel@data$indiv)-5))
#     title(loc[i], cex=1.2,sub = paste0("Total area sampled = ",
#                                        centroids[centroids$plotID==substr(loc[i],1,8),]$totalSampledAreaTrees, " m2"))
#     dev.off()
# 
# }else{
#     writeLines("no trees measured",paste0("./R_output/figures/woody_diameter_missing/", loc[i],"no_trees.txt"))}
# } ### 
# 
# 
# 
# 
