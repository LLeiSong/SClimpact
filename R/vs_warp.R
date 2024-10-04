library(dbarts)
library(checkmate)
library(terra)
library(sf)
library(dplyr)
library(optparse)
select <- dplyr::select
source("/home/lsong/SCImpact/R/utils.R")
source("/home/lsong/SCImpact/R/var_selection.R")

option_list <- list(
    make_option(c("-s", "--sp"), 
                action = "store", type = 'character',
                help = "The species to process."))
opt <- parse_args(OptionParser(option_list = option_list))
sp <- opt$sp
sp_list <- strsplit(sp, ",")[[1]]

root_dir <- "/home/lsong/SCImpact/data"
var_path <- file.path(root_dir, "variables/Env/AllEnv.tif")
occ_dir <- file.path(root_dir, "occurrences")
range_dir <- file.path(root_dir, "IUCN/Expert_Maps")
aoh_dir <- file.path(root_dir, "IUCN_AOH_100m/Mammals")
dst_dir <- file.path(root_dir, "variables/variable_list")
tmp_dir <- file.path(root_dir, "tmp")
seed <- as.integer(123)

for (sp in sp_list){
    var_selection(sp, var_path, occ_dir, range_dir, 
                  aoh_dir, dst_dir, tmp_dir, seed)
}
