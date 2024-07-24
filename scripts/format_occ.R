#' @title Format occurrences.
#' @description Clean the species occurrences using 
#' package \code{\link{BIENWorkflow}}.
#' @param occurrences (`sf`) The simple feature of occurrences.
#' @param range_map (`sf`) The range_map of this species.
#' @param vars (`SpatRaster`) The raster stack of variables.
#' @param stats (`list`) The stats setting for function spatialStratify.
#' @param dirs (`list`) The list of directories
#' @return The occurrence records.
#'
#' @seealso
#' \code{\link{cleanOcc}}, \code{\link{spatialStratify}}
#'
#' @details
#' Convert the terrestrial continents into hex grids. Note that the map will be
#' shifted to keep Russia together. So it will be problematic if the user 
#' directly transform to geographic coordinate system.
#'
#' @importFrom dplyr select rename filter
#' @importFrom BIENWorkflow cleanOcc spatialStratify
#' @importFrom sf st_drop_geometry st_join
#' @importFrom terra extract
#' @export
#' @examples
#' # Query occurrence for species Lerista ameles from GBIF.
#' occ <- query_gbif(sci_name = "Lerista ameles")
#' 
# Function to format the occurrences based on the grid layer template
format_occ <- function(occurrences, range_map, vars, stats, dirs){
    csv <- occurrences %>% dplyr::select() %>% 
        st_coordinates() %>% data.frame() %>% 
        dplyr::rename(lat = Y, lon = X)
    csv_nm <- file.path(tempdir(), sprintf("%s.csv", species))
    write.csv(csv, csv_nm, row.names = FALSE)
    
    thin_flag <- ifelse(nrow(csv) >= 20, TRUE, FALSE)
    prog.points <- cleanOcc(speciesCSV = csv_nm, env = vars, dirs = dirs, 
                            doThin = thin_flag, doClean = TRUE, 
                            writeResult = FALSE)
    file.remove(csv_nm)
    
    pres <- prog.points$pres
    # Clean the points a bit
    pres <- cbind(pres, terra::extract(vars, pres, ID = FALSE))
    keep <- complete.cases(sf::st_drop_geometry(pres))
    pres <- pres[keep, ]
    
    # Filter by expert range map
    opt_old <- sf_use_s2()
    sf_use_s2(FALSE)
    pres <- st_join(pres, range_map %>% mutate(presence = 1) %>% 
                        st_make_valid()) %>% 
        dplyr::filter(!is.na(presence)) %>% dplyr::select(-presence)
    
    sf_use_s2(opt_old)
    if (nrow(pres) >= 10){
        # Cluster to folds
        pres <- spatialStratify(pres, stats)
    }
    return(pres)
}

# Load extra libraries
library(BIENWorkflow)

# Directories
data_dir <- "data"
dst_dir <- file.path(data_dir, "gbif_clean")
if (!dir.exists(dst_dir)) dir.create(dst_dir)
misc_dir <- file.path(data_dir, "misc")
if (!dir.exists(misc_dir)) dir.create(misc_dir)

# Other settings
sf_use_s2(FALSE)

# Get species info: fname, species name, number of occurrences
## Only keep the species with more than 10 occurrences to create meaningful
## model

catalog <- read.csv(file.path(data_dir, "catalog_species.csv")) %>% 
    filter(class == "MAMMALIA")
occ_num <- list.files(file.path(data_dir, "gbif"), full.names = TRUE) %>% 
    data.frame(fname = .) %>% 
    mutate(species = gsub(".geojson", "", basename(fname))) %>% 
    mutate(species = gsub("_", " ", species)) %>% 
    filter(species %in% catalog$sci_name)
nums <- sapply(occ_num$fname, function(fname){
    st_read(fname) %>% nrow()
})
occ_num <- occ_num %>% mutate(num = nums) %>% 
    filter(num >= 10); rm(nums)

# Clean occurrences of filtered species
vars <- rast(file.path("data/env", "chelsa_1981-2010.tif"))
load("data/IUCN/mammals_endangered.rda")

dirs <- NULL
dirs$miscDir <- misc_dir

coarse.grid.file <- paste0(dirs$miscDir, '/coarse_grid_for_thinning.tif')
if(!file.exists(coarse.grid.file)){
    env.template <- vars[[1]]
    coarse.grid <- terra::aggregate(env.template, fact = 2)
    terra::writeRaster(coarse.grid, file = coarse.grid.file, overwrite = TRUE)
}

# Start to format the occurrences
for (i in 1:nrow(occ_num)){
    # Get the record
    record <- occ_num %>% slice(i)
    message(i)
    
    # Extract variables
    species <- record$species
    range_map <- range_maps %>% filter(sci_name == species) %>% 
        dplyr::select(presence)
    occ <- st_read(record$fname)
    
    # Format the occurrences
    stats <- .statsStorage(species, type = "PPM", toDo = NULL)
    occ_new <- format_occ(occ, range_map, vars, stats, dirs)
    
    if (nrow(occ_new) > 10){
        st_write(occ_new, file.path(dst_dir, sprintf("%s.geojson", gsub(" ", "_", species))))
    }
}

