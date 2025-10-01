######## Calculate crown areas from trees with different heights ########
####### crowns can be overlapping, contained within other crowns, etc. ##
###### AS April 2020, based on NEON data  ##############################

library(neonUtilities)
library(geoNEON)
library(sf)
library(rgeos)
library(viridis)
library(tidyverse)
library(gdata)

# ### Start with downloading data from NEON (see https://www.neonscience.org/neonDataStackR)
# # zipsByProduct(dpID="DP1.10098.001", site=c("ABBY","BART"),  
# #               startdate="2017-01", enddate="2018-11",
# #               package="basic", check.size=T)
# # stackByTable(filepath="./filesToStack10098/", savepath = "./",saveUnzippedFiles=T)
# 
# ### Load data
# maptag <- read.csv("./orig_data/NEON_struct-woody-plant/stackedFiles/vst_mappingandtagging.csv") ### location and plant ID
# indiv <- read.csv("./orig_data/NEON_struct-woody-plant/stackedFiles/vst_apparentindividual.csv")
# 
# maptag <- maptag %>% mutate(year=sapply(strsplit(as.character(maptag$date), "-"),"[",1))
# maptag <- maptag %>%  mutate(plotyear= paste(maptag$plotID, maptag$year, sep="_"))

### Optional: Select desired site(s) and year(s)
# patt <- "JERC"
# look <- maptag[grepl(patt, maptag$plotyear), ]
# 
# ### OR: Run everything (can take a while depending on the no of site/years)
# # look <- maptag
# 
# ### Calculate tree coordinates
# tree_coords <- getLocTOS(look, "vst_mappingandtagging")

# ### test
# cc <- tree_coords %>% filter(plotID=="JERC_047")%>% filter(adjEasting>0)
# 
# ccx <- cc %>%
#   left_join(indiv[,c(9,12:27)],by = "individualID")

# ### Merge with growth measurements
# tree_coords_combi <- tree_coords %>%
#   left_join(indiv[,c(9,12:27)],by = "individualID")%>% #### date needs to come from indiv !!!
#   mutate_if(is.factor, as.character)

### no coords
# tt <- tree_coords_combi %>% filter(siteID=="JERC")%>% filter(maxCrownDiameter>0)

#### leave out for JERC
# tree_coords_save <- read.csv("./R_output/woody/tree_coordinates.csv")
# 
# xx <- tree_coords_save[grepl("JERC",tree_coords_save$siteID),]

# testi <- tree_coords_save[grepl("NEON.PLA.D01.HARV.01929",tree_coords_save$individualID),]
# write.csv(testi,"./shrub_dubs.csv", row.names = F)

tree_coords_save <- read.csv("./R_output/woody/tree_coords/tree_coordinates_token.csv")

tree_coords_combi <- tree_coords_save %>% 
  filter(!(plantStatus=="Dead, broken bole"|plantStatus=="Standing dead"))


### Replace missing crown90Diameter with maxDiameter
tree_coords_combi[which(is.na(tree_coords_combi$ninetyCrownDiameter)),]$ninetyCrownDiameter <- tree_coords_combi[which(is.na(tree_coords_combi$ninetyCrownDiameter)),]$maxCrownDiameter

### Reduce input (for easier overview later)
tree_combi2 <- tree_coords_combi %>% 
  select("siteID","year","plotID", "taxonID", "height", "adjNorthing","adjEasting", "maxCrownDiameter", 
         "ninetyCrownDiameter", "individualID", "indiv", "plotyear") %>%
  distinct() %>% arrange(plotyear)

table(tree_combi2$siteID, tree_combi2$year)
ID <- sort(unique(tree_combi2$plotyear))

### Spatial polygons for plots
# # plotsize <- read.csv("./stackedFiles/vst_perplotperyear.csv")
# plotsize <- read.csv("//Volumes/NEON_2021/NEON_SD/orig_data/NEON_struct-woody-plant/stackedFiles/vst_perplotperyear.csv")
# plotsize$year <- sapply(strsplit(as.character(plotsize$date), "-"),"[",1)
# plotsize$plotyear <- paste(plotsize$plotID, plotsize$year, sep="_")
# 
# centroids <- plotsize %>% 
#   select(c(5,6,13,14,7,31:34,40:42)) %>% 
#   arrange(plotID)
# 
# ### Optional: Select site and year if needed 
# centroids_sel <- centroids[grepl(patt, centroids$plotyear),]
# 
# ### OR: run all
# # centroids_sel <- centroids
# 
# ### Spatial polygons for sites
# radius <- sqrt(centroids_sel$totalSampledAreaTrees)/2 # radius (m) for plots
# yPlus <- centroids_sel$northing+radius
# xPlus <- centroids_sel$easting+radius
# yMinus <- centroids_sel$northing-radius
# xMinus <- centroids_sel$easting-radius
# 
# # Calculate polygon coordinates for each plot centroid. 
# square=cbind(xMinus,yPlus,  # NW corner
#              xPlus, yPlus,  # NE corner
#              xPlus,yMinus,  # SE corner
#              xMinus,yMinus, # SW corner
#              xMinus,yPlus)  # NW corner again - close ploygon
# # Extract the plot ID information
# ID <- sort(as.character(centroids_sel$plotyear))
# 
# polys <- SpatialPolygons(mapply(function(poly, id){
#   xy <- matrix(poly, ncol=2, byrow=TRUE)
#   Polygons(list(Polygon(xy)), ID=id)
# }, 
# split(square, row(square)), ID))
# polydf <- SpatialPolygonsDataFrame(Sr = polys, data = centroids_sel,match.ID = F)

polydf <- readRDS("./R_input/Plot_polygons.rds")
pp <- as.data.frame(polydf)

anyDuplicated(polydf$plotID)
unique(tree_combi2$plotID[!(tree_combi2$plotID %in% polydf$plotID)])

# polydf$plotID[grepl("JERC", polydf$plotID)]
# ID <- c("SRER_048_2016", "SRER_056_2017","UNDE_078_2016", "YELL_051_2019")
# ID <- "BART_047_2015"

### Start loop ###
### Calculate crown area per plot (ID) 

for(k in 1:length(ID)){  
  keep(tree_combi2,ID,k,polydf,sure=T)
  errfile <- file(paste0("./R_output/woody/errors/error_file_crowns_per_plot_",ID[k], ".txt")) ### useful for loop
  tryCatch({
  seli <- tree_combi2[tree_combi2$plotyear==ID[k]&!(is.na(tree_combi2$adjNorthing)),]
  seli <- seli[order(seli$height,decreasing = F),] ### order sel by height
  seli <- seli[!(is.na(seli$maxCrownDiameter)),]
  seli <- seli[!(is.na(seli$height)),]
  
  if(!(anyDuplicated(seli[,c("adjNorthing","adjEasting")])==0)){
  seli <- seli[-which(duplicated(seli[,c("adjNorthing","adjEasting")])),]
  }
  
  ### Spatial object
  sel <- seli
  sel$individualID <- as.character(sel$individualID)
  sel <- sel[order(sel$height,decreasing = T),]
  
  if(!(anyDuplicated(sel$individualID)==0)){
    sel <- sel[-which(duplicated(sel$individualID)),]
  }
  
  coordinates(sel) <- c("adjEasting", "adjNorthing")### X before y
  
  elli_a <- ((sel@data$maxCrownDiameter)/2)*((sel@data$ninetyCrownDiameter)/2)*pi
  elli_r <- sqrt(elli_a/pi)
  
  buff <- gBuffer(sel,byid = T,width = elli_r,
                  id = sel@data$individualID, quadsegs=100) ### buffer by ellipse radius
  
  # plot(buff)
  buff@plotOrder <- order(buff@data$height)
  buffx <- buff
  
  remx <- gContainsProperly(buff,byid = T,returnDense = F)
  remx[sapply(remx, is.null)] <- NULL 
  
  #### Plot
  bo <- st_as_sf(polydf[polydf$plotID==substr(ID[k],1,8),])
  
  pdf(file = paste0("./R_output/figures/woody_diameter_byheight/", ID[k],"_diam_perheight.pdf"),width = 9)
  plot(bo$geometry, border=4)
  plot(buff, col=magma(nrow(sel),direction = -1, alpha=0.9,begin = 0.15),add=T)
  plot(sel, add=T)
  polygonsLabel(buff, substr(buff$indiv,1, nchar(buff$indiv)-5),
                method = "centroid",doPlot = T, cex=1, pos=2)
  plot(bo$geometry, border=4, add=T)
  title(ID[k], cex=1.2,sub = paste0("Total area sampled = ",
                                     bo$totalSampledAreaTrees, " m2"))
  dev.off()
  
  if(is_empty(remx)==F){
    ### check for crowns fully contained within other crowns
    try(for(i in 1:20){ ### run this a couple of times
      remx <- gContainsProperly(buff,byid = T,returnDense = F)
      remx[sapply(remx, is.null)] <- NULL 
      
      if(is_empty(remx)==F){
      rem_dat <- as.data.frame(unlist(remx))
      rem_dat$big <- row.names(rem_dat)
      names(rem_dat)[1] <- "small_ID"
      
      ssx <- substr(sapply(strsplit(rem_dat$big, "\\."),"[",5),1,5)
      rem_dat$big <- paste0(substr(rem_dat$big,1,18), ssx)
      
      for(i in 1:nrow(rem_dat)){
        rem_dat$small[i] <-  buff@data$individualID[rem_dat$small_ID[i]]
        rem_dat$height_small[i] <-buff@data$height[rem_dat$small_ID[i]]
        rem_dat$height_big[i] <- buff@data$height[which(buff$individualID==rem_dat$big[i])]  
      }
      
      rem_dat <- rem_dat[,c(2,5,3,4,1)] 
      rem <- rem_dat[rem_dat$height_big>rem_dat$height_small,] ### remove trees with smaller diameter that are less tall
      
      #### remove smaller plots
      buff <- subset(buff, !(buff$individualID %in% rem$small))
      }
    }, silent = T)
  }
  
  # plot(subset(buff,buff$individualID %in%rem$small), add=T, col=2)
  buffi <- buff
  
  # plot(buffi,col=rev(rainbow(n=length(buffi),alpha=0.5)))
  # polygonsLabel(buffi, substr(buffi$indiv,1, nchar(buffi$indiv)-5),
  #               method = "centroid", cex=0.5,doPlot = T)
  
  buffi_sf <- st_as_sf(buffi)
  plot(buffi_sf$geometry, col=rainbow(n=length(buff),alpha=0.5))
  
  ###### WORK ONLY AFTER SECOND MANIPULATION
  ### TEST
  # buffi_sf$indiv
  # plot(buffi_sf$geometry[c(3,4,6)], col=rainbow(n=length(buff),alpha=0.5))
  # 
  # mini <- buffi_sf %>% filter(indiv=="01066_PIST_2019"|indiv=="01058_QURU_2019"|indiv=="01091_QURU_2019")
  # mini <- mini %>% arrange(height)
  # combi <- mini %>% st_set_precision(1000) %>% st_intersection() 
  ###### END TEST
  
  
  ### Set precision (some won't run): 1000 == 1/1000 m (mm level)
  buffi_sf <- buffi_sf %>% arrange(height) ### try

  # pdf(file = paste0("./R_output/figures/woody_diameter_byheight_ordered/", ID[k],"_diam_perheight_ordered.pdf"),width = 9)
  ggplot(buffi_sf)+
    geom_sf(aes(fill=height),fill=magma(n=nrow(buffi_sf),direction = 1), alpha=0.8,col=1, show.legend = F)+
    geom_sf_label(aes(label = substr(indiv,1,10)), cex=2)+
    geom_sf(data=bo, fill=NA, col=2, lwd=2)
  
  ggsave(paste0("./R_output/figures/woody_diameter_byheight_ordered/", ID[k],"_diam_perheight_ordered.pdf"))
  # dev.off()
  
  combi <- buffi_sf %>% st_set_precision(1000) %>% st_intersection() 
 
  ### plots
  # for(x in 1:nrow(combi)){
  #   # pdf(file = paste0("./R_output/figures/woody_diameter_pieces/", ID[k],
  #   #                   "_piece_", combi$indiv[x],".pdf"),width = 9)
  #   plot(bo$geometry) 
  #   plot(combi$geometry[x],col=x, add=T)
  #   # title(paste0(combi$indiv[x], "_",x))
  #   title(paste(combi$indiv[combi$origins%in% unlist(combi[x,]$origins)],x,
  #               sep="_"))
  #   # dev.off()
  # }
  
 ### Or: Snap features to themselves
 # xx <- st_snap(buffi_sf, buffi_sf, tolerance = 0.1)
 # combix <- xx %>% st_set_precision(1000) %>% st_intersection()
  
  ### add origins for all tree segments
  unis <- buffi_sf %>% mutate(n.overlaps=1, origins=1:nrow(buffi_sf))  
  uu <- unis[!(unis$origins %in% combi$origins),]
  
  if(nrow(uu)>0){
    combi <- rbind(combi,uu)
  }
  
 ### Calculate area
  combi$area <- st_area(combi)
  
  if(nrow(uu)>0){
  combi[combi$origins %in% uu$origins,]$area <- 0
  }
  
  ### Add info about height as classes for tree selection
  combi2 <- as.data.frame(combi)
  combi2$group <- seq(1:nrow(combi2))

  for (i in 1: nrow(combi2)){
    combi2$groupx[i] <- paste(combi2[i,]$origins[[1]], collapse = "|")
  }

  combi2$heightmax <- numeric(nrow(combi2))
  for (i in 1:nrow(combi2)){
    if(combi2[i,]$n.overlaps==1){
      combi2[i,]$heightmax <- combi2[i,]$height
    } else{
      # aa <- buffi[unlist(combi2[i,]$origins),]@data   
      # pp <- paste(paste0("^",unlist(strsplit(combi2[i,]$groupx,"\\|"))), collapse="|")
      #   selx <- combi2[grepl(pp, combi2$groupx) & combi$n.overlaps==1,] 
        combi2[i,]$heightmax <- max(combi2$height[combi2$origins %in% unlist(combi2[i,]$origins)])
    }
  }
  
  ### TESTS
  # for(x in 1:nrow(combi)){
  #   plot(bo$geometry) 
  #   plot(combi2$geometry[x],col=x, add=T)
  #   title(paste0(combi2$indiv[x], "_",x))
  # }
  # 
  combi2$indsel <- character(nrow(combi))
  combi2$IND <- character(nrow(combi))
  
  for(i in 1:nrow(combi2)){
    # aa <- buffi[unlist(combi2[i,]$origins),]@data #### plotnames
    combi2[i,]$indsel <- paste(combi2$individualID[combi2$origins %in% unlist(combi2[i,]$origins)],
                               collapse = ",")
    if(combi2[i,]$n.overlaps>1){
      pp <- paste(paste0("^",unlist(strsplit(combi2[i,]$groupx,"\\|")),"$"), collapse="|")
      # patt <- paste(aa$individualID,collapse = "|")
      patt <-  paste(combi2$individualID[combi2$origins %in% unlist(combi2[i,]$origins)],
            collapse = "|")
      nnh <- combi2[grepl(pp, combi2$groupx),] ### for regular overlaps
      groupi <- nnh[nnh$heightmax==max(nnh$heightmax),"group"]
      
      # ppp <- gsub("\\^", "", pp) ### for overlaps within crowns
      # ppp <- gsub("\\$", "", ppp)
      # ppp <- paste0("^", ppp, "$")
      # nnh2 <- combi2[grepl(ppp, combi2$groupx),]
      # nnh <- bind_rows(nnh2,nnh)
      # groupi <- nnh[nnh$heightmax==max(nnh$heightmax),"group"]
      # aa$heightmax <- max(aa$height)  ### add if clause for crown contained in crowns
  
      if(length(groupi)>1){
        combi2$IND[i] <- patt
      }else{
        combi2[i,]$group <- groupi

        indi <- nnh[nnh$heightmax==max(nnh$heightmax),grep(pattern = "individualID",colnames(combi2))]
        combi2$IND[i] <- indi[!(is.na(indi))]
        
        # 
        # indi <- aa[aa$height==max(aa$height),grep(pattern = "individualID",colnames(combi2))]
        # combi2$IND[i] <- indi[!(is.na(indi))] 
        
      }
    }else{
      indi <-  as.character(combi2[i,grep(pattern = "individualID",colnames(combi2))])
      combi2$IND[i] <- indi[!(is.na(indi))]
    }
  }

### And another round for crowns contained in crowns  
outx <- combi2
outxx <- st_as_sf(outx)

# plot(outxx$geometry)
st_agr(outxx) <- "constant"
st_agr(bo)  <- "constant"

outxxx <- st_intersection(outxx,bo)
plot(outxxx$geometry, col=rainbow(n=length(outx),alpha=0.5))
# polygonsLabel(buffi, buffi$individualID,method = "centroid", cex=0.5,doPlot = T)

outxxx$area <- st_area(outxxx)


if(nrow(uu)>0){
  if(any(outxxx$origins%in% uu$origins)==TRUE){  ### new addition
    outxxx[outxxx$origins %in% uu$origins,]$area <- 0 
  }
}


any(outxxx$origins%in% uu$origins)

### Add single tree info for trees containing trees
combi2 <- outxxx 

#### same height
if(exists("rem_dat")){
  sam <- rem_dat[rem_dat$height_big-rem_dat$height_small==0,]
  ss <- c(paste(sam$big, sam$small, sep=","),paste(sam$small, sam$big, sep=","))
  
  combi2[combi2$indsel %in% ss,]$IND <- gsub(",","\\|",combi2[combi2$indsel %in% ss,]$indsel)
  
  ### smaller tree is taller
  tall <- rem_dat[rem_dat$height_big-rem_dat$height_small<0,]
  tt <- c(paste(tall$big, tall$small, sep=","),paste(tall$small, tall$big, sep=","))
  
  xx <- combi2[combi2$indsel %in% tt,]
  
  for(i in 1:nrow(xx)){
    ll <- unlist(strsplit(xx$indsel[i],","))
    xx$IND[i] <- tall[tall$small%in%ll,]$small
  }
  
  # if(length(tt)>0 & any(combi2$indsel %in% tt)){
  #   combi2 <- combi2[-(which(combi2$indsel %in% tt)),]  
  #   combi2 <- rbind(combi2,xx) 
  # }
  if(length(tt)>0 & any(combi2$indsel %in% tt)){
    if(any(combi2$indsel %in%tt)){
      combi2 <- combi2[-(which(combi2$indsel %in% tt)),]
    }
  }
  if(nrow(xx)>0){
    for(u in 1:nrow(xx)){
      if(!(xx$IND[u] %in% combi2$IND))
        combi2 <- rbind(combi2,xx)
    }
  }
}

### Plot fragments
for(u in 1:nrow(combi2)){
  # pdf(file = paste0("./R_output/figures/woody_diameter_pieces/", ID[k],
  #                   "_piece_", u,".pdf"),width = 9)
  plot(bo$geometry) 
  plot(combi2$geometry[u],col=u, add=T)
  # title(paste0(combi$indiv[x], "_",x))
  title(paste(combi2$IND[u],u,
              sep="_"),sub = paste(combi2$area[u],"m2"))
  # dev.off()
}

#### split combi2 into 2 parts
combi2 <- as.data.frame(as_tibble(combi2))

### Save crown segments if needed (usually not)
saveRDS(combi2, paste0("./R_output/woody/segments/" ,ID[k],"_crown_segments.rds"))


### Part 1 with clear assignment[IND], Part 2 with trees sharing same height
combi_b <- combi2[nchar(combi2$IND)>23,]
combi_a <- combi2[nchar(combi2$IND)==23,]

# out_a <- combi_a %>% group_by(IND)%>% summarize(area_total=max(area))

ax <- combi_a[grep(pattern = "\\|",combi_a$groupx,invert = T),]

# for(u in 1:nrow(ax)){
#   # pdf(file = paste0("./R_output/figures/woody_diameter_pieces/", ID[k],
#   #                   "_piece_", u,".pdf"),width = 9)
#   plot(bo$geometry) 
#   plot(ax$geometry[u],col=u, add=T)
#   # title(paste0(combi$indiv[x], "_",x))
#   title(paste(ax$IND[u],u,
#               sep="_"),sub = paste(ax$area[u],"m2"))
#   # dev.off()
# }

# out <- combi2 %>% group_by(IND)%>% summarize(area_total=sum(area))
# combi2$origins <- as.character(combi2$origins)
# # combi2 <- subset(combi2 , select=-c(geometry))
# 
# plot(combi2$geometry)

### Overlaps with same height

bx <- list()
if(nrow(combi_b)>0){ ### combi_b overlapping areas
  for (i in 1:nrow(combi_b)){
    combi_b[i,]$IND <- paste(buffi_sf[buffi_sf$height==combi_b[i,]$heightmax,]$individualID,
          collapse ="|") 
  }
  ub <- unique(combi_b$IND)
  
  for(j in 1:length(ub)){ ### How many overlaps
    foo <- combi_b[combi_b$IND==ub[j],]
    mino <- min(foo$n.overlaps) 
    bx[[j]] <- foo[foo$n.overlaps==mino,]
  }
  
  bxx <- do.call(rbind,bx)
  idx <- paste(bxx$IND,collapse = "|")

  axx <- ax[!(grepl(idx,ax$IND)),]
  
  for(i in 1:nrow(bxx)){
    frag_size <- bxx[i,]$area
    frag_size_frac <- bxx[i,]$area/bxx[i,]$n.overlaps
    
    ff <- ax[grepl(bxx[i,]$IND, ax$IND),]
    ff$origins <- as.numeric(ff$origins)
    
    ff[ff$origins==max(ff$origins),]$area <- ff[ff$origins==max(ff$origins),]$area - frag_size + frag_size_frac
    ff[!(ff$origins==max(ff$origins)),]$area <-  ff[!(ff$origins==max(ff$origins)),]$area + frag_size_frac
  
    # plot(bo$geometry)
    # plot(ax[ax$origins==5,]$geometry, add=T)
    out  <- rbind(axx,ff)
    }
} else{
  out <- ax
}

outfin <- out %>% filter(area>0) %>% arrange(origins)

for(u in 1:nrow(outfin)){
  pdf(file = paste0("./R_output/figures/woody_diameter_pieces/", ID[k],
                    "_piece_", u,".pdf"),width = 9)
  plot(bo$geometry) 
  plot(outfin$geometry[u],col=u, add=T)
  # title(paste0(combi$indiv[x], "_",x))
  title(paste(outfin$indiv[u],u,
              sep="_"),sub = paste(outfin$area[u],"m2"))
  dev.off()
}

saveRDS(outfin, paste0("./R_output/woody/crown_area/" ,ID[k],"_crowns_geo.rds"))
# plot(bxx$geometry)
# cc <- out[grepl("\\|", out$IND),]$IND 

# if(length(cc)>0){
#   for (i in 1:length(cc)){ 
#     namx <- unlist(strsplit(cc[i],"\\|"))
#     # safe <- selx
#     selx <- bind_rows(combi2[combi2$IND ==cc[i],],combi2[combi2$IND %in% namx,]) ### plots to work with
#     # selx <- as_tibble(selx)
#     oo <- selx[grepl("\\|", selx$IND),] ### split plots
#     keep <- selx[!(grepl("\\|", selx$IND)),] ### decided plots to keep 
#     
#     ### three or more overlpas with same height: remove id's of trees that are not at max height from split
#     nammax <- unique(selx[(selx$height %in% max(selx$height,na.rm=T)) &selx$n.overlaps==1,"IND"])
#     
#     #### modify output
#     area_part <- sum(oo[,"area"]/length(nammax))
#     
#     for (j in 1:length(nammax)){
#       out[out$IND==nammax[j],"area_total"] <- out[out$IND==nammax[j],"area_total"]+area_part
#     }
#   }
# }

# out_fin <- out[!(grepl("\\|", out$IND)),]
# outfin2 <- merge(out_fin, seli[, c("individualID","indiv")],by.x="IND", 
#                  by.y="individualID", all.x=T)

### Save final result
# str(outfin)

outfin2 <- subset(outfin , select=-c(geometry, origins))

write.csv(outfin2, paste0("./R_output/woody/crown_area/" ,ID[k],"_crowns.csv"))

}, error=function(e){writeLines(as.character(e),errfile)})
close(errfile)
  }

