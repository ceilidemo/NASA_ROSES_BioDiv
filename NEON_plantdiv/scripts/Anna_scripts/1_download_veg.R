library(httr)
library(jsonlite)
library(dplyr)
library(downloader)
library(neonUtilities)

#### Currently API interface downloads released and provisional data
### Download from web and stack programmatically

# req <- GET("http://data.neonscience.org/api/v0/products/DP1.10058.001")
# req.content <- content(req, as="parsed")
# names(req.content$data)
# req.content$data$releases
# 
# req.text <- content(req, as="text")
# avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)

# tok <- "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbm5hLmsuc2Nod2VpZ2VyQGdtYWlsLmNvbSIsInNjb3BlIjoicmF0ZTpwdWJsaWMiLCJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImV4cCI6MTc2OTY5MDExMiwiaWF0IjoxNjEyMDEwMTEyLCJlbWFpbCI6ImFubmEuay5zY2h3ZWlnZXJAZ21haWwuY29tIn0.zSwI3pydhBnxXjZNP4S3QgFODsepLB-UPg9QYJ_bqE1OdIBuUga7Q2bQbFiDmLOg152GtsO0cZQmIxDhbmZBiA"
# zipsByProduct(dpID="DP1.10058.001",site = "BART",token = tok, 
#               package="expanded", check.size=T,savepath = "./orig_data/")

stackByTable(filepath="./orig_data/NEON_presence-cover-plant/", 
             saveUnzippedFiles = F)


#### Woody
zipsByProduct(dpID="DP1.10098.001",
              package="expanded", check.size=T,savepath = "./R_input/plants/")
stackByTable(filepath="./R_input/plants/filesToStack10098/", 
             folder=T,saveUnzippedFiles = T)


#### Chemistry
zipsByProduct(dpID="DP1.10026.001",
              package="expanded", check.size=T,savepath = "./R_input/plants/")
stackByTable(filepath="./R_input/plants/filesToStack10026/", folder=T,saveUnzippedFiles = T)
