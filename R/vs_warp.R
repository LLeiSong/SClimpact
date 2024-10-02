library(dbarts)
library(checkmate)
library(terra)
library(sf)
library(dplyr)
library(optparse)
select <- dplyr::select
source("R/utils.R")
source("R/variable_selection.R")

option_list <- list(
    make_option(c("-s", "--sp"), 
                action = "store", type = 'character',
                help = "The species to process."))
opt <- parse_args(OptionParser(option_list = option_list))
sp <- opt$sp

root_dir <- "data"
var_path <- file.path(root_dir, "variables/Env/AllEnv.tif")
occ_dir <- file.path(root_dir, "occurrences")
aoh_dir <- file.path(root_dir, "IUCN_AOH_100m/Mammals")
dst_dir <- file.path(root_dir, "variables/variable_list")
tmp_dir <- file.path(root_dir, "tmp")
seed <- as.integer(123)

var_selection(sp, var_path, occ_dir, aoh_dir, dst_dir, tmp_dir, seed)