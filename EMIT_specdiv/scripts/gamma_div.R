####GAMMA Diversity Function Adapted from Ettienne's specdiv_functions.R####

gamma_diversity <- function(cube, fact = 40) {
  require(raster)
  require(tidyverse)
  
  # Get cube with plots (i.e., count pixels)
  n_pixels <- count_pixels(cube, fact = fact)
  
  # Prepare to store results
  gamma_ss <- double()
  gamma_sdiv <- double()
  gamma_fcsd <- matrix(nrow = 1, ncol = dim(cube)[3], dimnames = list(1, names(cube)))
  
  # Convert raster to points and extract only the values (skip x and y)
  cube_points <- rasterToPoints(cube, spatial = FALSE) %>% 
    as_tibble() 
  
  # Check if the 'layer' column exists
  # If not, we just assume the values are in the columns beyond the first two (x, y)
  value_columns <- colnames(cube_points)[3:ncol(cube_points)]  # assuming x, y are the first two columns
  
  # Now, cube_points should only contain the pixel data (excluding x, y)
  cube_points_sel_gamma <- cube_points %>%
    dplyr::select(all_of(value_columns))  # Select only the pixel values
  
  # Compute gamma diversity using sum of squares function
  sdiv_gamma <- sum_squares(cube_points_sel_gamma)  # Calculate sum of squares for gamma diversity
  gamma_ss[1] <- sdiv_gamma$ss
  gamma_sdiv[1] <- sdiv_gamma$sdiv
  gamma_fcsd[1, ] <- sdiv_gamma$fcsd
  
  # Prepare output
  output <- list(
    gamma_ss = gamma_ss,
    gamma_sdiv = gamma_sdiv,
    gamma_fcsd = gamma_fcsd
  )
  
  return(output)
}

### SS gamma and alpha ----
sum_squares <- function(Y) {
  n <- nrow(Y)
  Y.cent <- scale(Y, center = T, scale = F)
  sij <- Y.cent^2
  SS.total <- sum(sij)
  SS.row <- rowSums(sij)
  SS.col <- colSums(sij)
  fcsd <- SS.col / SS.total
  lcsd <- SS.row / SS.total
  sdiv <- SS.total / (n - 1)
  out <- list(ss = SS.total, sdiv = sdiv, lcsd = lcsd, fcsd = fcsd)
  return(out)
}

### Count masked/missing pixels -----
count_pixels <- function(cube, fact = 40) {
  require(raster)
  require(tidyverse)
  # Get plots (or communities)
  new_res <- res(cube) * fact
  cube_plots <- raster(crs = proj4string(cube))
  extent(cube_plots) <- extent(cube)
  res(cube_plots) <- new_res
  cube_plots <- setValues(cube_plots, values = 1:ncell(cube_plots))
  plot_xy <- rasterToPoints(cube_plots)[,1:2]
  cube_pixels <- disaggregate(cube_plots, fact = fact)
  # Convert to points
  plot_points <- rasterToPoints(cube_pixels) %>% 
    as_tibble() %>% 
    rename(group = layer)
  cube_points <- rasterToPoints(cube) %>% 
    as_tibble() %>% 
    right_join(plot_points, by = c('x', 'y'))
  n_pixels <- cube_points %>% 
    group_by(group) %>% 
    do(count_noNA(.)) %>% 
    mutate(n_total = fact^2,
           prop = n / n_total) %>% 
    ungroup()
  group_df <- SpatialPixelsDataFrame(plot_xy, dplyr::select(n_pixels, n:prop), proj4string = CRS(proj4string(cube) ))
  cube_count <- brick(group_df)
  return(cube_count)
}

#count no NA rows
count_noNA <- function(x) {
  n <- nrow(na.omit(x))
  n <- as.data.frame(n)
}