# Prepare environmental variables
# Lei Song, lsong@ucsb.edu

# Load libraries
library(terra)
library(dplyr)
library(stringr)

dir_chelsa <- '/Volumes/gazelle/CHELSA_V2/GLOBAL/climatologies'
dst_dir <- "data/env"
if (!dir.exists(dst_dir)) dir.create(dst_dir)
years <- c("1981-2010", "2011-2040", "2041-2070", "2071-2100")
mod <- "GFDL-ESM4"
ssp <- 'ssp370'

# BIO1, BIO2, BIO12, BIO15, BIO10/(BIO10 + BIO11)
vars <- c("bio1", "bio2", "bio12", "bio15", "bio10", "bio11")
for (yr in years){
    if (yr == "1981-2010"){
        data_dir <- file.path(dir_chelsa, yr, "bio")
        lyrs <- rast(file.path(data_dir, sprintf("CHELSA_%s_%s_V.2.1.tif", vars, yr)))
    } else {
        data_dir <- file.path(dir_chelsa, yr, mod, ssp, "bio")
        lyrs <- rast(file.path(data_dir, sprintf("CHELSA_%s_%s_%s_%s_V.2.1.tif", 
                                                 vars, yr, tolower(mod), ssp)))
    }
    
    names(lyrs) <- vars
    lyrs$bio120 <- lyrs$bio10 / (lyrs$bio10 + lyrs$bio11)
    lyrs <- subset(lyrs, c(1:4, 7))
    writeRaster(lyrs, file.path(dst_dir, sprintf("chelsa_%s.tif", yr)))
}

# land cover
template <- rast("data/env/chelsa_1981-2010.tif", lyr = 1)
lc_dir <- "/Volumes/gazelle/future_land_projection_chen/Global\ 7-land-types\ LULC\ projection\ dataset\ under\ SSPs-RCPs"
ssp <- "SSP3_RCP70"
lc_2015 <- rast(file.path(lc_dir, "global_LULC_2015.tif"))
lc_2015 <- project(lc_2015, template, method = "near")
writeRaster(lc_2015, file.path(dst_dir, "LULC_2015.tif"), datatype = "INT1U")

for (yr in c(2040, 2070, 2100)){
    lc <- rast(file.path(lc_dir, "SSP3_RCP70", sprintf("global_SSP3_RCP70_%s.tif", yr)))
    lc <- project(lc, template, method = "near")
    writeRaster(lc, file.path(dst_dir, sprintf("LULC_%s.tif", yr)), datatype = "INT1U")
}
