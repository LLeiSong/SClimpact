library(stringr)
library(terra)
library(sf)
library(spatstat)
library(dplyr)
select <- dplyr::select
filter <- dplyr::filter
mask <- terra::mask
source("/home/lsong/SCImpact/R/climate_change.R")

root_dir <- "/home/lsong/SCImpact"
species_list <- read.csv(
    file.path(root_dir, "data/occurrences", "species_qualified_sdm.csv"))
species_list <- species_list$species

evals <- lapply(species_list, function(sp){
    dr <- file.path(root_dir, "results/sdm", sp)
    eval <- read.csv(file.path(dr, sprintf("cv_eval_%s.csv", sp)))
    eval <- eval %>% 
        mutate(accuracy = (tp + tn) / (tp + tn + fp + fn))
    eval <- eval %>% select(-fold) %>% group_by(type) %>% 
        summarise_all("mean") %>% 
        mutate(species = sp)
    if (all(eval$auc >= 0.7)){
        eval
    } else NULL
}) %>% bind_rows()

species_list <- unique(evals$species)

species_catalog <- read.csv(
    file.path(root_dir, "data/occurrences", "species_qualified_sdm.csv")) %>% 
    filter(species %in% species_list)

fname <- file.path(root_dir, "results", "species_reliable_final.csv")
if (!file.exists(fname)){
    write.csv(species_catalog, fname, row.names = FALSE)
}; rm(evals, species_catalog, fname)

scenarios <- c(
    "ssp126_2011-2040", "ssp126_2041-2070", "ssp126_2071-2100",
    "ssp370_2011-2040", "ssp370_2041-2070", "ssp370_2071-2100",
    "ssp585_2011-2040", "ssp585_2041-2070", "ssp585_2071-2100")

sdm_dir <- file.path(root_dir, "results/sdm")
dst_dir <- file.path(root_dir, "results/species_analysis")
if(!dir.exists(dst_dir)) dir.create(dst_dir)

for (sp in species_list){
    climate_change_sp(sp, scenarios, root_dir, sdm_dir, dst_dir)
}
