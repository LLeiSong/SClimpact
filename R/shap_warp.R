library(randomForest)
library(checkmate)
library(terra)
library(sf)
library(dplyr)
library(intervals)
library(blockCV)
library(fastshap)
library(doParallel)
library(optparse)
select <- dplyr::select
source("/home/lsong/SCImpact/R/utils.R")
source("/home/lsong/SCImpact/R/shap.R")

option_list <- list(
    make_option(c("-s", "--sp"), 
                action = "store", type = 'character',
                help = "The species to process."))
opt <- parse_args(OptionParser(option_list = option_list))
sp <- opt$sp
sp_list <- strsplit(sp, ",")[[1]]

root_dir <- "/home/lsong/SCImpact"
var_dir <- file.path(root_dir, "data/variables")
range_dir <- file.path(root_dir, "data/IUCN/Expert_Maps")
work_dir <- file.path(root_dir, "results/sdm")

for (sp in sp_list){
    tryCatch({
        shap(sp, work_dir, var_dir, range_dir)
    }, error = function(e){print(e)})
}
