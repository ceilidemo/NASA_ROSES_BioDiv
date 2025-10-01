##Download data needed for plant div and RTM##

library(neonUtilities)
  library(dplyr)

# DP1.10058.001 = plant presence percent cover
# DP1.10026.001 = plant foliar traits
# DP1.10098.001 = vegetation structure
# DP3.30010.001 = lidar derived structure
# DP1.10033.001 = LAI 

# Loading all of the vegetation structure datasets for all sites all years
structure_all <- loadByProduct(
  dpID = "DP1.10098.001",
  site = "all",
  startdate = "2017-01",
  enddate = "2023-12",
  package = "basic",  
  check.size = FALSE
)

names(structure_all) # Pull out names of individual data

vst_apparentindividual <- structure_all$vst_apparentindividual  
head(vst_apparentindividual) # Apparentindividual has info on biomass and productivity of woody individuals... no geolocation other than plot
write.csv(vst_apparentindividual, "data_in/allyr/vst_apparentindividual.csv")

vst_mappingandtagging <- structure_all$vst_mappingandtagging
head(vst_mappingandtagging) # Mappingandtagging has info on ID and mapping of individual stems? I think it is stem distance and azimuth
write.csv(vst_mappingandtagging, "data_in/allyr/vst_mappingandtagging.csv")

vst_perplotperyear <- structure_all$vst_perplotperyear
head(vst_perplotperyear) # Perplotperyear just has info on the baseplots and when they were sampled
write.csv(vst_perplotperyear, "data_in/allyr/vst_perplotperyear.csv")


# Loading all of the plant presence data - hoping to find geolocation of individuals that were ID'd in vst
presence_all <- loadByProduct(
  dpID = "DP1.10058.001",
  site = "all",
  startdate = "2017-01",
  enddate = "2023-12",
  package = "basic",  
  check.size = FALSE
)

names(presence_all) # Pull out data names again

div_1m2Data <- presence_all$div_1m2Data
head(div_1m2Data)  #  Plant species identifications and cover, and cover of abiotic variables within 1 square meter SUBPLOT
write.csv(div_1m2Data, "data_in/allyr/div_1m2Data.csv")


div_10m2Data100m2Data <- presence_all$div_10m2Data100m2Data
head(div_10m2Data100m2Data) # Plant species identifications within 10 square meter and 100 square meter subplots
write.csv(div_10m2Data100m2Data, "data_in/allyr/div_10m2Data100m2Data.csv")



vst_apparentindividual <- read.csv("data_in/allyr/vst_apparentindividual.csv")
vst_mappingandtagging <- read.csv("data_in/allyr/vst_mappingandtagging.csv")
vst_perplotperyear <- read.csv("data_in/allyr/vst_perplotperyear.csv")
div_1m2Data <- read.csv("data_in/allyr/div_1m2Data.csv")
div_10m2Data100m2Data <- read.csv("data_in/allyr/div_10m2Data100m2Data.csv")


# The way that NEON does tree mapping / structure it would not make sense to filter to 2023
# So instead filter out to just get all the mapped stems for example
library(stringr)

# get just the distributed plots
plot_summary <- vst_perplotperyear%>%
  filter(str_detect(plotType, "distributed")) 

# Here we are gonna get the stems that are only in the distributed base plots
stems_mapped <- vst_mappingandtagging%>% 
  #filter(str_detect(recordType, "map and tag")) %>%
  filter(plotID %in% plot_summary$plotID) %>% 
  distinct(individualID, .keep_all = TRUE)

# get the structure info for the individuals that were mapped in the distributed plots
vst_individual <- vst_apparentindividual %>%
  filter(plotID %in% plot_summary$plotID) %>% 
  filter(individualID %in% stems_mapped$individualID) #%>% 
  mutate(date = as.Date(date)) %>% 
  group_by(individualID) %>%
  slice_max(order_by = date, n = 1, with_ties = FALSE) %>%
  ungroup()

# Join the mapped stems with their structure
mapped_vst <- left_join(stems_mapped, vst_individual, by = "individualID")


## Lettuce try and map these bad boys
library(ggplot2)
library(ggforce)
library(purrr)

# get dataset set with trees for mapping and their locations from stem azimuth and distance
mapped_trees <- mapped_vst %>%
  filter(!is.na(height), !is.na(maxCrownDiameter)) %>%
  mutate(
    crown_radius = maxCrownDiameter / 2,
    azimuth_rad = stemAzimuth * pi / 180,
    x = stemDistance * sin(azimuth_rad),
    y = stemDistance * cos(azimuth_rad)
  ) %>%
  group_by(plotID.x) %>%
  arrange((height)) %>%
  mutate(height_rank = row_number()) %>%
  ungroup() 


# List of unique plots
plots <- unique(mapped_trees$plotID.x)
plot_extent <- 20 # for 40x40 m baseplots

for(p in plots){
  
  plot_data <- mapped_trees %>% filter(plotID.x == p)
  
  p_map <- ggplot(plot_data) +
    geom_circle(aes(x0 = x, y0 = y, r = crown_radius, 
                    fill = height, color = height_rank), alpha = 0.4) +
    coord_fixed(xlim = c(-plot_extent, plot_extent),
                ylim = c(-plot_extent, plot_extent)) +
    scale_fill_viridis_c(option = "viridis") +
    scale_color_viridis_c(option = "plasma", direction = -1) +
    labs(title = paste("Tree Crowns - Plot", p),
         fill = "Height (m)", color = "Height Rank") +
    theme_minimal() +
    theme(panel.grid = element_blank())  # optional: remove grid
  
  ggsave(filename = paste0("tree_maps/", p, ".png"), plot = p_map,
         width = 6, height = 6, dpi = 300)
}




#foliar trait processing
lma <- plants_2022$cfc_LMA
cn <- plants_2022$cfc_carbonNitrogen
elements <- plants_2022$cfc_elements
field <- plants_2022$cfc_fieldData

lma_avg <- lma %>%
  group_by(sampleID) %>%
  summarize(
    # keep one of each of these — should be consistent within sampleID
    siteID = first(siteID),
    domainID = first(domainID),
    namedLocation = first(namedLocation),
    plotType = first(plotType),
    plotID = first(gsub("\\..*", "", namedLocation)),  # extract plotID from namedLocation
    tissueDescriptor = first(tissueDescriptor),
    collectDate = first(collectDate),
    # average all numeric columns
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

cn_avg <- cn %>%
  group_by(sampleID) %>%
  summarize(

    siteID = first(siteID),
    domainID = first(domainID),
    namedLocation = first(namedLocation),
    plotID = first(plotID),
    plotType = first(plotType),
    collectDate = first(collectDate),
    analysisDate = first(analysisDate),
    testMethod = first(testMethod),
    laboratoryName = first(laboratoryName),
    instrument = first(instrument),
    
    carbonPercent = mean(carbonPercent, na.rm = TRUE),
    nitrogenPercent = mean(nitrogenPercent, na.rm = TRUE),
    CNratio = mean(CNratio, na.rm = TRUE),
    d13C = mean(d13C, na.rm = TRUE),
    d15N = mean(d15N, na.rm = TRUE),
    
    .groups = "drop"
  )

elements_avg <- elements %>%
  group_by(sampleID) %>%
  summarize(
    siteID = first(siteID),
    domainID = first(domainID),
    namedLocation = first(namedLocation),
    plotID = first(plotID),
    plotType = first(plotType),
    collectDate = first(collectDate),
    analysisDate = first(analysisDate),
    laboratoryName = first(laboratoryName),
    testMethod = first(testMethod),
    
    foliarPhosphorusConc = mean(foliarPhosphorusConc, na.rm = TRUE),
    foliarPotassiumConc  = mean(foliarPotassiumConc,  na.rm = TRUE),
    foliarCalciumConc    = mean(foliarCalciumConc,    na.rm = TRUE),
    foliarMagnesiumConc  = mean(foliarMagnesiumConc,  na.rm = TRUE),
    foliarSulfurConc     = mean(foliarSulfurConc,     na.rm = TRUE),
    foliarManganeseConc  = mean(foliarManganeseConc,  na.rm = TRUE),
    foliarIronConc       = mean(foliarIronConc,       na.rm = TRUE),
    foliarCopperConc     = mean(foliarCopperConc,     na.rm = TRUE),
    foliarBoronConc      = mean(foliarBoronConc,      na.rm = TRUE),
    foliarZincConc       = mean(foliarZincConc,       na.rm = TRUE),
    
    .groups = "drop"
  )

field_avg <- field %>%
  group_by(sampleID) %>%
  summarize(
    domainID = first(domainID),
    siteID = first(siteID),
    namedLocation = first(namedLocation),
    plotID = first(plotID),
    plotType = first(plotType),
    collectDate = first(collectDate),
    individualID = first(individualID),
    taxonID = first(taxonID),
    scientificName = first(scientificName),
    plantStatus = first(plantStatus),
    canopyStatus = first(canopyStatus),
    canopyPosition = first(canopyPosition),
    leafSamplePosition = first(leafSamplePosition),
    decimalLatitude = first(decimalLatitude),
    decimalLongitude = first(decimalLongitude),
    elevation = first(elevation),
    samplingMethod = first(samplingMethod),
    .groups = "drop"
  )

traits_combined <- lma_avg %>%
  full_join(cn_avg, by = "sampleID") %>%
  full_join(elements_avg, by = "sampleID") %>%
  left_join(field_avg, by = "sampleID")

write.csv(traits_combined, "data_in/foliar_traits/foliar_traits_2022.csv")



