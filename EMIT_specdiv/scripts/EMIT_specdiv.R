## EMIT processing and Spectral Gamma Diversity ## 
library(sp)
library(raster)
library(tidyverse)
library(ncdf4)

#############################
#####DATAPROCESSING##########
#############################

#open netcdf files and explore info
dat <- nc_open("data_in/enveloped_sampling/2023/YELL/EMIT_L2A_RFL_001_20230622T175609_2317312_003_ortho_subset_aid0000.nc") 
print(names(dat$var))
wavelengths <- ncvar_get(dat, "wavelengths") #extract wavelengths
wavelengths_df <- data.frame(wavelengths = wavelengths) #create dataframe

wvl <- wavelengths_df
wvl <- wvl %>%
  mutate(band_name = 1:285)

hyperspectral_cube <- ncvar_get(dat, "reflectance") #extract reflectance values

bands <- dim(hyperspectral_cube)[3]
hyperspectral_stack <- stack(lapply(1:bands, function(i) {
  raster(hyperspectral_cube[, , i])
}))
writeRaster(hyperspectral_stack, "working/cubes/UKSF_hyperspectral_cube.tif", format = "GTiff", overwrite = TRUE) #write raster of reflectance
nc_close(dat)

# Hyperspectral data cube : assign wavelengths and plot to check cloud coverage and data
cube <- brick("working/cubes/KONZ_hyperspectral_cube.tif")
names(cube) <- wvl$band_name
plotRGB(cube,
        r = 36, g = 26, b = 14, stretch = 'lin')


#########################
#####GAMMADIVERSITY######
#########################

source("scripts/gamma_div.R") #source Gamma Diversity function

library(raster)
library(tidyverse)

cube <- rastercube <- brick('working/cubes/ABBY_hyperspectral_cube.tif')
#checks
#cube_points <- rasterToPoints(cube, spatial = FALSE) %>% as_tibble()
#print(head(cube_points))
#colnames(cube_points)

cube_norm <- bright_norm(cube)

result <- gamma_diversity(cube_norm)
result$gamma_ss ##Get Gamma Sum Squares ! 
result$gamma_sdiv ##Get Spectral Gamma Div 

#some error fixes for small sites
#res(cube) <- c(min(res(cube)), min(res(cube))) #ensure consistent resolution
#cube <- projectRaster(cube, crs = CRS("+proj=utm +zone=33 +datum=WGS84")) #reproject?
#extent(cube) <- c(0, 1, 0, 1)  # Update as necessary #for tiny tiny extent ex: 0 to 1

