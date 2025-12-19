## 01_SpecDiv
# Compute spectral gamma diversity

# --- Source functions ---
source("scripts/gamma_div_functions.R")

# --- Load cube ---
cube <- brick("working/cubes/YELL_hyperspectral_cube.tif")

# --- Normalize brightness (function assumed in gamma_div_functions.R) ---
cube_norm <- bright_norm(cube)

# --- Compute gamma diversity ---
result <- gamma_diversity(cube_norm, target_res = 60)

# --- Print results ---
cat("Gamma Sum of Squares:", result$gamma_ss, "\n")
cat("Spectral Gamma Diversity:", result$gamma_sdiv, "\n")