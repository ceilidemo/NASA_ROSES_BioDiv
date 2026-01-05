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

# filter to have the latest for each site
latest_cover_data <- cover_all$div_1m2Data %>%
  group_by(siteID, plotID, subplotID) %>%
  slice_max(endDate, n = 1, with_ties = FALSE) %>%
  ungroup()

latest_structure_data <- structure_all$vst_apparentindividual %>%
  group_by(siteID, plotID, individualID) %>%
  slice_max(date, n = 1, with_ties = FALSE) %>%
  ungroup()
