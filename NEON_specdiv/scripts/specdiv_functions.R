#### Spectral diversity functions
# Etienne Laliberté, February 2019


### Brightness normalization -----
bright_norm <- function(x) {
  x_norm <- x / sqrt(sum(x^2))
  return(x_norm) 
}


### Brightness normalization for data frame ----
bright_norm_df <- function(x) {
  x$refl_norm <- x$refl / sqrt(sum(x$refl^2))
  return(x)  
}


### Short PCA for matrix-----
pca_mat <- function(x, scaling = c(1, 2), p = 0.99) {
  Y <- scale(x, center = TRUE, scale = FALSE)
  n <- nrow(Y)
  Y.svd = svd(Y)
  values = (1/(n - 1)) * Y.svd$d^2
  epsilon = sqrt(.Machine$double.eps)
  k <- sum(values > epsilon)
  values <- values[1:k]
  prop <- values / sum(values)
  cumprop = cumsum(prop)
  if (p < cumprop[1]) which.values <- c(1, 2) else which.values <- which(cumprop < p)
  values.sel <- values[c(which.values, length(which.values) + 1)]
  n.pcs <- length(values.sel)
  U <- as.matrix(Y.svd$v[, 1:n.pcs])
  if (scaling == 1) {
    obj <- Y %*% U
    descript <- U
  }
  else {
    obj <- sqrt(n - 1) * as.matrix(Y.svd$u[, 1:n.pcs])
    descript <- U %*% diag(values.sel^(0.5))
  }
  colnames(obj) <- paste0('PC', 1:n.pcs)
  colnames(descript) <- paste0('PC', 1:n.pcs)
  rownames(descript) <- colnames(x)
  prop.sel <- prop[1:n.pcs]; names(prop.sel) <- colnames(descript)
  cumprop.sel <- cumprop[1:n.pcs]; names(cumprop.sel) <- colnames(descript)
  out <- list(obj = obj, descript = descript, prop = prop.sel, cumprop = cumprop.sel)
  return(out)
}



### PCA for image cube -----
pca <- function(cube, scaling = 1, p = 0.99) {
  
  pts_df <- rasterToPoints(cube, spatial = FALSE) 
  xy <- pts_df[, 1:2]
  pixels <- pts_df[, -c(1, 2)]
  
  complete_idx <- complete.cases(pixels)
  xy <- xy[complete_idx, ]
  pixels <- pixels[complete_idx, ]
  
  Y <- scale(pixels, center = TRUE, scale = FALSE)
  n <- nrow(Y)
  
  Y.svd <- svd(Y)
  values <- (1 / (n - 1)) * Y.svd$d^2
  epsilon <- sqrt(.Machine$double.eps)
  k <- sum(values > epsilon)
  values <- values[1:k]
  prop <- values / sum(values)
  cumprop <- cumsum(prop)
  
  if (p < cumprop[1]) {
    which.values <- c(1, 2)
  } else {
    which.values <- which(cumprop < p)
  }
  values.sel <- values[c(which.values, length(which.values) + 1)]
  n.pcs <- length(values.sel)
  
  U <- as.matrix(Y.svd$v[, 1:n.pcs])
  
  if (scaling == 1) {
    obj <- Y %*% U
    descript <- U
  } else {
    obj <- sqrt(n - 1) * as.matrix(Y.svd$u[, 1:n.pcs])
    descript <- U %*% diag(values.sel^(0.5))
  }
  
  colnames(obj) <- paste0('PC', 1:n.pcs)
  colnames(descript) <- paste0('PC', 1:n.pcs)
  rownames(descript) <- colnames(pixels)
  
  print(raster::crs(cube))
  cube_pc_df <- cbind(x = xy[, 1], y = xy[, 2], as.data.frame(obj))
  crs_string <- as.character(raster::crs(cube))
  cube_pc <- raster::stack(raster::rasterFromXYZ(cube_pc_df, crs = raster::crs(cube)))
  
  prop.sel <- prop[1:n.pcs]
  names(prop.sel) <- colnames(descript)
  cumprop.sel <- cumprop[1:n.pcs]
  names(cumprop.sel) <- colnames(descript)
  
  out <- list(cube_pc = cube_pc, band_contrib = descript, prop = prop.sel, cumprop = cumprop.sel)
  return(out)
}





### Make PC plots ----
make_plot <- function(df) {
  x <- ggplot(df, aes(x = x, y = y) ) +
    geom_raster(aes(fill = value)) +
    scale_fill_gradientn(colors = rainbow(20)) +
    ggtitle(label = unique(df$PC)) +
    theme_void() +
    theme(legend.position = 'none',
          plot.title = element_text(face = 'bold', hjust = 0.5, size = 15) ) +
    coord_equal()
  return(x)
}


### Build the PC plots ----
build_pc_plot <- function(cube) {
  require(tidyverse)
  points <- rasterToPoints(cube, spatial = F) %>% 
    as_tibble() %>% 
    gather(key = PC, value = value, -x, -y) %>% 
    mutate(PC = factor(PC, levels = paste0('PC', 1:n_distinct(PC)) ) ) 
  cube_plots <- points %>%
    group_by(PC) %>% 
    do(plots = make_plot(df = .))
  return(cube_plots)
}


### Show the PC plots ----
show_pc_plot <- function(x, ...) {
  require(ggplot2)
  require(gridExtra)
  grid.arrange(grobs = x$plots, ...) 
}


# Count non-NA rows
count_noNA <- function(x) {
  n <- nrow(na.omit(x))
  n <- as.data.frame(n)
}

### Count masked/missing pixels -----
count_pixels <- function(cube, fact = 40) {
  require(raster)
  require(tidyverse)
  require(sp)
  
  new_res <- res(cube) * fact
  cube_plots <- raster(crs = proj4string(cube))
  extent(cube_plots) <- extent(cube)
  res(cube_plots) <- new_res
  cube_plots <- setValues(cube_plots, values = 1:ncell(cube_plots))
  
  plot_xy <- rasterToPoints(cube_plots)[, 1:2] %>%
    matrix(ncol = 2) %>%
    as.data.frame()
  colnames(plot_xy) <- c("x", "y")
  
  cube_pixels <- disaggregate(cube_plots, fact = fact)

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
  
  plot_xy_sp <- SpatialPoints(plot_xy, proj4string = CRS(proj4string(cube)))
  

  coords_df <- as.data.frame(plot_xy)
  coordinates(plot_xy) <- ~ x + y
  proj4string(plot_xy) <- CRS(proj4string(cube))
  
  group_df <- SpatialPixelsDataFrame(
    coords_df,
    dplyr::select(n_pixels, n:prop),
    tolerance = 0.5, # helps avoid tiny floating point misalignments
    proj4string = CRS(proj4string(cube))
  )
  
  cube_count <- brick(group_df)
  
  return(cube_count)
}


### SS gamma and alpha ----
sum_squares <- function(Y) {
  Y <- as.matrix(Y)
  Y <- Y[complete.cases(Y), , drop = FALSE]  # remove NA rows if any
  
  n <- nrow(Y)
  if (n < 2) {
    return(list(ss = NA_real_, sdiv = NA_real_, lcsd = NA, fcsd = NA))
  }
  
  Y.cent <- scale(Y, center = TRUE, scale = FALSE)
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



### SS beta ----
sum_squares_beta <- function(Y, m) {
  n <- nrow(Y)
  Y.cent <- bind_cols(dplyr::select(Y, group), as.data.frame(scale(dplyr::select(Y, -group), scale = F)))
  mskj <- Y.cent %>% 
    group_by(group) %>% 
    mutate_at(vars(-group_cols()), function(x) (mean(x))^2) %>% 
    summarise_at(vars(-group_cols()), sum) %>% 
    ungroup() %>% 
    dplyr::select(-group)
  SSbk <- rowSums(mskj)
  SSbj <- colSums(mskj)
  SSb <- sum(SSbk)
  sdiv <- SSb / (n - 1)
  fcsd <- SSbj / SSb
  lcsd <- SSbk / SSb
  out <- list(ss = SSb, sdiv = sdiv, lcss = SSbk, lcsd = lcsd, fcsd = fcsd)
  return(out)
}


### Partitioning spectral diversity ----
specdiv <- function(cube, fact = 40, prop = 0.5, n = 1) {
  require(raster)
  require(tidyverse)
  require(sf)
  require(stars)
  
  # Extract CRS from raster cube and convert to sf crs
  crs_proj4 <- raster::crs(cube)
  sf_crs <- sf::st_crs(crs_proj4@projargs)
  
  # Find minimum number of pixels for resampling
  n_pixels <- count_pixels(cube, fact = fact)
  
  # Remove plots with proportion of pixels less than prop
  plot_mask <- subset(n_pixels, 'prop') < prop
  n_pixels_mask <- raster::mask(n_pixels, plot_mask, maskvalue = 1)
  min_pixels <- minValue(n_pixels_mask)
  
  # Prepare a raster with plot IDs
  nlayers <- raster::nlayers(cube)
  cube_plots <- raster::raster(crs = crs_proj4)
  raster::extent(cube_plots) <- raster::extent(n_pixels_mask)
  raster::res(cube_plots) <- raster::res(n_pixels_mask)
  cube_plots <- raster::setValues(cube_plots, values = 1:raster::ncell(cube_plots))
  cube_plots_masked <- raster::mask(cube_plots, plot_mask, maskvalue = 1)
  cube_pixels <- raster::disaggregate(cube_plots_masked, fact = fact)
  
  # Convert masked plots and pixels rasters to data frames
  plots_points <- raster::rasterToPoints(cube_plots_masked) %>% 
    as_tibble() %>% 
    rename(group = layer)
  pixels_points <- raster::rasterToPoints(cube_pixels) %>% 
    as_tibble() %>% 
    rename(group = layer)
  
  # Convert plots_points to sf points
  xy_plots_sf <- sf::st_as_sf(plots_points %>% select(x, y), coords = c("x", "y"), crs = sf_crs)
  
  # Prepare storage for results
  gamma_ss <- numeric(n)
  gamma_sdiv <- numeric(n)
  gamma_fcsd <- matrix(nrow = n, ncol = nlayers, dimnames = list(NULL, names(cube)))
  alpha_sdiv <- tibble()
  alpha_fcsd <- tibble()
  alpha_ss <- tibble()
  beta_ss <- numeric(n)
  beta_sdiv <- numeric(n)
  beta_fcsd <- matrix(nrow = n, ncol = nlayers, dimnames = list(NULL, names(cube)))
  beta_lcsd <- tibble()
  beta_lcss <- tibble()
  
  # Main loop for sampling and calculations
  for (i in seq_len(n)) {
    cube_points <- raster::rasterToPoints(cube) %>% 
      as_tibble() %>% 
      inner_join(pixels_points, by = c("x", "y")) %>% 
      group_by(group) %>% 
      sample_n(size = min_pixels) %>% 
      ungroup()
    
    # Gamma diversity
    sdiv_gamma <- sum_squares(dplyr::select(cube_points, -x, -y, -group))
    gamma_ss[i] <- sdiv_gamma$ss
    gamma_sdiv[i] <- sdiv_gamma$sdiv
    gamma_fcsd[i, ] <- sdiv_gamma$fcsd
    
    # Alpha diversity
    sdiv_alpha <- cube_points %>%
      dplyr::select(-x, -y) %>%
      group_by(group) %>%
      summarise(res = list(sum_squares(dplyr::select(., -group)))) %>%
      ungroup()
    
    alpha_sdiv_tmp <- tibble(rep = i, group = sdiv_alpha$group, sdiv = map_dbl(sdiv_alpha$res, "sdiv"))
    alpha_fcsd_tmp <- tibble(rep = i, group = sdiv_alpha$group) %>% 
      bind_cols(as_tibble(do.call(rbind, map(sdiv_alpha$res, "fcsd"))))
    alpha_ss_tmp <- tibble(rep = i, group = sdiv_alpha$group, ss = map_dbl(sdiv_alpha$res, "ss"))
    
    alpha_sdiv <- bind_rows(alpha_sdiv, alpha_sdiv_tmp)
    alpha_fcsd <- bind_rows(alpha_fcsd, alpha_fcsd_tmp)
    alpha_ss <- bind_rows(alpha_ss, alpha_ss_tmp)
    
    # Beta diversity
    sdiv_beta <- sum_squares_beta(dplyr::select(cube_points, -x, -y), m = min_pixels)
    beta_ss[i] <- sdiv_beta$ss
    beta_sdiv[i] <- sdiv_beta$sdiv
    beta_fcsd[i, ] <- sdiv_beta$fcsd
    
    beta_lcsd_tmp <- tibble(rep = i, group = seq_along(sdiv_beta$lcsd), lcsd = sdiv_beta$lcsd)
    beta_lcss_tmp <- tibble(rep = i, group = seq_along(sdiv_beta$lcss), lcss = sdiv_beta$lcss)
    beta_lcsd <- bind_rows(beta_lcsd, beta_lcsd_tmp)
    beta_lcss <- bind_rows(beta_lcss, beta_lcss_tmp)
  }
  
  # Summarize and convert beta LCSD to sf + stars raster
  lcsd_beta_values <- beta_lcsd %>%
    group_by(group) %>%
    summarise(lcsd_beta = mean(lcsd)) %>%
    ungroup()
  
  lcsd_beta_sf <- sf::st_as_sf(cbind(sf::st_coordinates(xy_plots_sf), lcsd_beta_values %>% select(-group)),
                               coords = c("X", "Y"), crs = sf_crs)
  lcsd_beta_raster <- stars::st_as_stars(lcsd_beta_sf)
  
  # Summarize and convert beta LCSS to sf + stars raster
  lcss_beta_values <- beta_lcss %>%
    group_by(group) %>%
    summarise(lcss_beta = mean(lcss)) %>%
    ungroup()
  
  lcss_beta_sf <- sf::st_as_sf(cbind(sf::st_coordinates(xy_plots_sf), lcss_beta_values %>% select(-group)),
                               coords = c("X", "Y"), crs = sf_crs)
  lcss_beta_raster <- stars::st_as_stars(lcss_beta_sf)
  
  # Summarize and convert FCSD alpha to sf + stars raster
  fcsd_alpha_values <- alpha_fcsd %>%
    group_by(group) %>%
    summarise(across(where(is.numeric), mean)) %>%
    ungroup()
  
  fcsd_alpha_sf <- sf::st_as_sf(cbind(sf::st_coordinates(xy_plots_sf), fcsd_alpha_values %>% select(-group)),
                                coords = c("X", "Y"), crs = sf_crs)
  fcsd_alpha_brick <- stars::st_as_stars(fcsd_alpha_sf)
  fcsd_alpha_mean <- colMeans(fcsd_alpha_values %>% select(-group))
  
  # SS summaries
  ss_alpha_sum <- alpha_ss %>%
    group_by(rep) %>%
    summarise(ss = sum(ss)) %>%
    summarise(mean_ss = mean(ss)) %>%
    pull(mean_ss)
  
  ss_beta <- mean(beta_ss)
  ss_gamma <- mean(gamma_ss)
  
  # SDiv summaries
  sdiv_beta <- mean(beta_sdiv)
  sdiv_gamma <- mean(gamma_sdiv)
  sdiv_alpha_values <- alpha_sdiv %>%
    group_by(group) %>%
    summarise(mean_sdiv = mean(sdiv)) %>%
    ungroup()
  
  sdiv_alpha_sf <- sf::st_as_sf(cbind(sf::st_coordinates(xy_plots_sf), sdiv_alpha_values %>% select(-group)),
                                coords = c("X", "Y"), crs = sf_crs)
  sdiv_alpha_raster <- stars::st_as_stars(sdiv_alpha_sf)
  sdiv_alpha_mean <- mean(sdiv_alpha_values$mean_sdiv)
  
  # Output lists
  ss <- tibble(source = c("alpha", "beta", "gamma"),
               sum_squares = c(ss_alpha_sum, ss_beta, ss_gamma)) %>%
    mutate(prop_gamma = sum_squares / ss_gamma)
  
  sdiv <- c(mean_alpha = sdiv_alpha_mean, beta = sdiv_beta, gamma = sdiv_gamma)
  
  fcsd <- bind_rows(
    tibble(source = "mean_alpha", t(fcsd_alpha_mean)),
    tibble(source = "beta", t(colMeans(beta_fcsd))),
    tibble(source = "gamma", t(colMeans(gamma_fcsd)))
  )
  
  rasters <- list(
    beta_lcsd = lcsd_beta_raster,
    beta_lcss = lcss_beta_raster,
    alpha_sdiv = sdiv_alpha_raster,
    alpha_fcsd = fcsd_alpha_brick
  )
  
  list(ss = ss, sdiv = sdiv, fcsd = fcsd, rasters = rasters)
}

