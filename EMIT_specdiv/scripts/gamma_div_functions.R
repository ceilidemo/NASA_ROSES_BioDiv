## Gamma Diversity Functions (simplified for EMIT cubes at 60m)

# --- Brightness normalization for raster cubes ---
bright_norm <- function(cube) {
  require(raster)
  cube_norm <- calc(cube, fun = function(x) {
    if (all(is.na(x))) return(NA)
    x / sqrt(sum(x^2, na.rm = TRUE))
  })
  return(cube_norm)
}

# --- Sum of Squares ---
sum_squares <- function(Y) {
  Y <- as.matrix(Y)   # force numeric matrix
  n <- nrow(Y)
  if (n < 2) return(list(ss = NA, sdiv = NA, fcsd = rep(NA, ncol(Y))))
  
  Y.cent <- scale(Y, center = TRUE, scale = FALSE)
  sij <- Y.cent^2
  SS.total <- sum(sij, na.rm = TRUE)
  SS.col <- colSums(sij, na.rm = TRUE)
  fcsd <- SS.col / SS.total
  sdiv <- SS.total / (n - 1)
  list(ss = SS.total, sdiv = sdiv, fcsd = fcsd)
}

# --- Gamma Diversity ---
gamma_diversity <- function(cube) {
  require(raster)
  require(tidyverse)
  
  cube_points <- rasterToPoints(cube, spatial = FALSE) %>% as_tibble()
  value_columns <- colnames(cube_points)[3:ncol(cube_points)]
  cube_points_sel <- cube_points %>% dplyr::select(all_of(value_columns))
  
  # Coerce to numeric matrix
  cube_mat <- as.matrix(cube_points_sel)
  
  sdiv_gamma <- sum_squares(cube_mat)
  
  output <- list(
    gamma_ss = sdiv_gamma$ss,
    gamma_sdiv = sdiv_gamma$sdiv,
    gamma_fcsd = sdiv_gamma$fcsd
  )
  return(output)
}
