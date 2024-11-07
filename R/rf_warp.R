library(blockCV)
library(precrec)
library(randomForest)
library(checkmate)
library(terra)
library(sf)
library(dplyr)
library(optparse)
select <- dplyr::select
source("/home/lsong/SCImpact/R/utils.R")
source("/home/lsong/SCImpact/R/rf_dws.R")

option_list <- list(
    make_option(c("-s", "--sp"), 
                action = "store", type = 'character',
                help = "The species to process."))
opt <- parse_args(OptionParser(option_list = option_list))
sp <- opt$sp
sp_list <- strsplit(sp, ",")[[1]]

root_dir <- "/home/lsong/SCImpact"
var_dir <- file.path(root_dir, "data/variables")
occ_dir <- file.path(root_dir, "data/occurrences")
range_dir <- file.path(root_dir, "data/IUCN/Expert_Maps")
dst_dir <- file.path(root_dir, "results/sdm")
seed <- as.integer(123)

for (sp in sp_list){
    tryCatch({
        rf_dws(sp, occ_dir, var_dir, range_dir, dst_dir)
    }, error = function(e){print(e)})
}
