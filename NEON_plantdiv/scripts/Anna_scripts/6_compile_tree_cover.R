library(tidyverse)

#### check plots: tower plots too?
# herbi <- read.csv("./R_output/inv_compiled/herb_cover_allplots.csv")
look <- read.csv("./R_input/NEON_plantcodes_lookup.csv")

fls <- list.files("./R_output/woody/crown_area/", full.names = T,pattern = ".csv")
trees_all <- map(fls,read.csv, row.names=1)

# trees_all2 <- lapply(trees_all, function(x) transform(x,plotID=as.character(plotID)))
trees_all2 <- purrr::discard(trees_all, ~nrow(.) == 0)

trees <- bind_rows(trees_all2)
anyNA(trees$area_total)

# trees$taxonID[grepl("BEGL", trees$taxonID)]

### combine with sp matched

trees2 <- trees %>% left_join(look[,c(1,6)]) %>% select(c(1:4,27,5:26))

df_trees <- trees2 %>%
  group_by(plotID, year, sp_matched,easting, northing, plotType, totalSampledAreaTrees, 
           totalSampledAreaShrubSapling, plotyear)  %>%
  summarise(area_m2 = sum(area)) %>%  
  mutate(cover_perc = area_m2/totalSampledAreaTrees*100)

dada <- list()
loc <- unique(df_trees$plotyear)

for (i in 1:length(loc)){
  # errfile <- file(paste0("./R_output/woody/errors_compile/compile_trees",loc[i], ".txt"))
  # tryCatch({
  dada[[i]]<- spread(data = df_trees[df_trees$plotyear==loc[i],c("plotyear", "sp_matched", "cover_perc")],
                     key = sp_matched, value = cover_perc)
  #   },error=function(e){writeLines(as.character(e),errfile)})
  # close(errfile)
}

tree_cov <- bind_rows(dada) %>% select(plotyear,order(colnames(.)))
colnames(tree_cov)[!(colnames(tree_cov) %in% look$sp_matched)]

tree_cov[is.na(tree_cov)] <- 0
write.csv(tree_cov, "./R_output/inv_compiled/tree_cover_allplots_token.csv", row.names = F)


#############################
##### Compare both files
xx <- read.csv("./R_output/inv_compiled/tree_cover_allplots.csv")
tree_cov2 <- read.csv("./R_output/inv_compiled/tree_cover_allplots_token.csv")

names(tree_cov2)[!(names(tree_cov2)%in% names(tree_cov))]

names(tree_cov)[!(names(tree_cov)%in% names(xx))]
names(xx)[!(names(xx)%in% names(tree_cov))]

xx$plotyear[!(xx$plotyear %in%tree_cov$plotyear)]
tree_cov$plotyear[!(tree_cov$plotyear %in%xx$plotyear)]

miss <- tree_cov[!(tree_cov$plotyear %in%xx$plotyear),] 

missi <- miss %>% select(where(~ any(. != 0)))


