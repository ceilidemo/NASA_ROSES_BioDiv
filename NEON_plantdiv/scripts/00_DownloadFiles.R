## 00_DownloadFiles: Download data needed for plant div and RTM##

library(neonUtilities)
library(dplyr)

# DP1.10058.001 = plant presence percent cover
# DP1.10098.001 = vegetation structure

# Load datasets 
cover_all <- loadByProduct(
  dpID = "DP1.10058.001", 
  site = "all", 
  enddate = "2023-12", 
  package = "basic", 
  check.size = FALSE,
  include.provisional = TRUE
)

structure_all <- loadByProduct(
  dpID = "DP1.10098.001",
  site = "all",
  enddate = "2023-12",
  check.size = FALSE,
  include.provisional = TRUE
)

#save it all
save(structure_all, cover_all, file = "data_in/veg_data.RData")

