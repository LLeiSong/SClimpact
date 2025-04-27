library(stringr)
library(terra)
library(sf)
library(dplyr)
library(parallel)
library(optparse)
select <- dplyr::select
filter <- dplyr::filter
source("/home/lsong/SCImpact/R/climate_change.R")

# Get model evaluation
root_dir <- "/home/lsong/SCImpact"
sdm_dir <- "/bigscratch/lsong/results/sdm"
dst_dir <- "/bigscratch/lsong/climate_change"
if (!dir.exists(dst_dir)) dir.create(dst_dir)

option_list <- list(
    make_option(c("-f", "--feature"),
                action = "store", type = 'character',
                help = "The climatic feature to process."))
opt <- parse_args(OptionParser(option_list = option_list))
feature <- opt$feature

species_list <- read.csv(file.path(
    root_dir, "results", "species_reliable_final.csv"))
species_list <- species_list$species

# Climate change analysis
template <- rast(
    file.path(root_dir,"data/variables/Env", "AllEnv.tif"), lyrs = 1)
values(template) <- NA

scenarios <- c(
    "ssp126_2011-2040", "ssp126_2041-2070", "ssp126_2071-2100",
    "ssp370_2011-2040", "ssp370_2041-2070", "ssp370_2071-2100",
    "ssp585_2011-2040", "ssp585_2041-2070", "ssp585_2071-2100")

for (scenario in scenarios){
    climate_change(
        feature, scenario, species_list, 
        template, root_dir, sdm_dir, dst_dir)
}
