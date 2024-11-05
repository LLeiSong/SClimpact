library(checkmate)
library(terra)
library(sf)
library(dplyr)
library(spatialEco)
library(pbmcapply)
library(optparse)
select <- dplyr::select
source("/home/lsong/SCImpact/R/utils.R")
source("/home/lsong/SCImpact/R/env_sampling.R")

option_list <- list(
    make_option(c("-s", "--species"), 
                action = "store", type = 'character',
                help = "The species to process."))
opt <- parse_args(OptionParser(option_list = option_list))
species <- opt$species
sp_list <- strsplit(species, ",")[[1]]

# Start the run
root_dir <- "/home/lsong/SCImpact"
src_dir <- file.path(root_dir, "data/occurrences/CSVs")
var_dir <- file.path(root_dir, "data/variables")
range_dir <- file.path(root_dir, "data/IUCN/Expert_Maps")
region_dir <- file.path(root_dir, "data/terr-ecoregions-TNC")
occ_dir <- file.path(root_dir, "data/occurrences/CSVs_thin")
bg_dir <- file.path(root_dir, "data/occurrences/bg")

env_sampling(sp_list, src_dir, var_dir, range_dir, 
             region_dir, occ_dir, bg_dir, 60000, 123, 30)

