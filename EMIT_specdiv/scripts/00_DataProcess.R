## 00_DataProcess
# Processing EMIT NetCDF file to tif 

library(sp)
library(raster)
library(tidyverse)
library(ncdf4)

# --- Input file ---
nc_file <- "data_in/enveloped_sampling/2023/YELL/EMIT_L2A_RFL_001_20230622T175609_2317312_003_ortho_subset_aid0000.nc"

# --- Open NetCDF ---
dat <- nc_open(nc_file)
print(names(dat$var))

# --- Extract wavelengths ---
wavelengths <- ncvar_get(dat, "wavelengths")
wvl <- data.frame(wavelengths = wavelengths) %>%
  mutate(band_name = paste0("Band_", 1:nrow(.)))

# --- Extract reflectance cube ---
hyperspectral_cube <- ncvar_get(dat, "reflectance")
bands <- dim(hyperspectral_cube)[3]

# --- Build raster stack ---
hyperspectral_stack <- stack(lapply(1:bands, function(i) {
  raster(hyperspectral_cube[, , i])
}))
names(hyperspectral_stack) <- wvl$band_name

# --- Save to GeoTIFF ---
out_tif <- "working/cubes/YELL_hyperspectral_cube.tif"
writeRaster(hyperspectral_stack, out_tif, format = "GTiff", overwrite = TRUE)

# --- Close NetCDF ---
nc_close(dat)

# --- Quick RGB check ---
cube <- brick(out_tif)
plotRGB(cube, r = 36, g = 26, b = 14, stretch = 'lin')