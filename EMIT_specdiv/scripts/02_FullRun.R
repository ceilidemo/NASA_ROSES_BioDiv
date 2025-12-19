## 02_FullRun
#Get spectral diversity for all sites

library(raster)
library(tidyverse)
library(ncdf4)

# source the div functions (simplified gamma only)
source("scripts/gamma_div_functions.R")

# load the tile info
emit_csv <- read_csv("data_in/EMIT_files2023.csv")

# base path
base_path <- "C:/Users/ceilidemarais/Projects/NASA_ROSES/EMIT_specdiv/data_in/enveloped_sampling/2023"

# results list
results <- list()

for (i in 1:nrow(emit_csv)) {
  site <- emit_csv$site[i]
  file <- emit_csv$file[i]
  
  if (is.na(file)) {
    message("Skipping site with missing file: ", site)
    next
  }
  
  nc_file <- file.path(base_path, site, file)
  if (!file.exists(nc_file)) {
    warning("File not found for site: ", site)
    next
  }
  
  message("Processing site: ", site)
  
  tryCatch({
    # read NetCDF
    dat <- nc_open(nc_file)
    wavelengths <- ncvar_get(dat, "wavelengths")
    wvl <- data.frame(wavelengths = wavelengths) %>%
      mutate(band_name = paste0("Band_", 1:nrow(.)))
    
    hyperspectral_cube <- ncvar_get(dat, "reflectance")
    bands <- dim(hyperspectral_cube)[3]
    
    hyperspectral_stack <- stack(lapply(1:bands, function(b) {
      raster(hyperspectral_cube[, , b])
    }))
    names(hyperspectral_stack) <- wvl$band_name
    
    out_tif <- file.path("data_work/new_cubes", paste0(site, "_hyperspectral_cube.tif"))
    dir.create(dirname(out_tif), recursive = TRUE, showWarnings = FALSE)
    writeRaster(hyperspectral_stack, out_tif, format = "GTiff", overwrite = TRUE)
    nc_close(dat)
    
    # gamma diversity
    cube <- brick(out_tif)
    cube_norm <- bright_norm(cube)
    result <- gamma_diversity(cube_norm)
    
    results[[site]] <- result
  }, error = function(e) {
    warning("Error processing site ", site, ": ", conditionMessage(e))
  })
}

# export results
results_df <- bind_rows(lapply(names(results), function(site) {
  data.frame(site = site,
             gamma_ss = results[[site]]$gamma_ss,
             gamma_sdiv = results[[site]]$gamma_sdiv)
}))

dir.create("data_out/results", recursive = TRUE, showWarnings = FALSE)
write_csv(results_df, "data_out/results/EMIT_SpecDiv_summary.csv")

print(results_df)
