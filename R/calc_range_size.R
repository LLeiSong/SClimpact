# Load some markdown required packages
library(sf)
library(terra)
library(dplyr)
library(parallel)

root_dir <- "/scratch/ls1686/SClimpact"
sdm_dir <- file.path(root_dir, "results/sdm")
range_dir <- file.path(root_dir, "data/IUCN/Expert_Maps")
ncores <- 10

# Get a list of species
species_list <- read.csv(
    file.path(root_dir, "results/species_reliable_final.csv"))

b_brd <- mclapply(species_list$species, function(sp){
    # Load range polygon
    range <- st_read(
        file.path(range_dir, sprintf("%s.geojson", sp)), 
        quiet = TRUE)
    
    # Calculate range size and put together
    data.frame(
        species = sp, range_size = as.numeric(sum(st_area(range))) / 1000000)
}, mc.cores = ncores) %>% bind_rows()

write.csv(b_brd, file.path(root_dir, "results", "range_size_species.csv"),
          row.names = FALSE)