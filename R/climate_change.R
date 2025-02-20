#' @title climate_change
#' @description Function to calculate climate change impact.
#' @param feature (`character`) The target environmental feature.
#' @param scenario (`character`) The climate model scenario.
#' @param species_list (`vector`) A vector of species.
#' @param template (`spatRaster`) The template to burn the result in.
#' @param root_dir (`character`) The directory for variable list.
#' @param sdm_dir (`character`) The directory for sdm of species.
#' @param dst_dir (`character`) The destination directory to save out the results.
#' @return All files will be saved out.
#' 
#' @details
#' Make sure GDAL library is installed on your machine to make this work.
#' 
#' @import checkmate
#' @import terra
#' @import sf
#' @import dplyr
#' @import stringr
#' @export
#' @examples
#' \donttest{
#' feature <- "bio1"
#' scenario <- "ssp370_2041-2070"
#' species_list <- c("Mustela_erminea", "Mustela_nivalis")
#' root_dir <- "/home/lsong/SCImpact"
#' sdm_dir <- "/bigscratch/lsong/results/sdm"
#' dst_dir <- "/bigscratch/lsong/climate change"
#' climate_change(feature, scenario, species_list, template, 
#' root_dir, sdm_dir, dst_dir)
#' }
#' 

climate_change <- function(feature, 
                           scenario, 
                           species_list, 
                           template,
                           root_dir,
                           sdm_dir,
                           dst_dir){
    # Species list
    species_use <- lapply(species_list, function(sp){
        var_list <- read.csv(
            file.path(root_dir, "data/variables/variable_list",
                      sprintf("%s.csv", sp))) %>% 
            pull(var_uncorrelated) %>% na.omit()
        if (feature %in% var_list){
            sp
        } else NULL
    }) %>% unlist() %>% na.omit()
    
    message(sprintf("%s species for %s.", length(species_use), feature))
    
    fname <- file.path(
        dst_dir, sprintf("species_use_%s_%s.txt", feature, scenario))
    write.table(species_use, fname,
                row.names = FALSE, col.names = "species")
    
    ## Get file names for these species
    fnames <- lapply(species_use, function(sp){
        list.files(file.path(sdm_dir, sp), pattern = "shap_base", 
                   full.names = TRUE)
    }) %>% unlist()
    
    # Load dispersal rate and continents
    dispersal_rates <- read.csv(
        file.path(root_dir, "data/variables/species_dispersal_rate.csv"))
    continents <- st_read(
        file.path(root_dir, "data/variables/continents.geojson"), quiet = TRUE)
    
    ## Initialize layers
    dir_change <- c(rep(template, 4))
    names(dir_change) <- c("N to N", "N to P", "P to N", "P to P")
    num_species <- c(rep(template, 3))
    names(num_species) <- c("All", "P", "N")
    val_change <- c(rep(template, 2))
    names(val_change) <- c("P", "N")
    
    for (fname in fnames){
        # Load layers
        cur <- rast(fname) %>% subset(feature)
        fut <- rast(gsub("base", scenario, fname)) %>% subset(feature)
        
        # Load the range for clipping layers
        sp <- rev(strsplit(dirname(fname), "/")[[1]])[1]
        range <- st_read(
            file.path(root_dir, "data/IUCN/Expert_Maps", 
                      sprintf("%s.geojson", sp)), quiet = TRUE) %>% 
            st_transform(crs(cur))
        dispersal_rate <- dispersal_rates %>% 
            filter(species == gsub("_", " ", sp)) %>% pull(dispersal_rate)
        continent <- continents %>% 
            slice(unique(unlist(st_intersects(range, .))))
        
        # Make the mask
        end_year <- as.integer(strsplit(scenario, "_|-")[[1]][3])
        dispersal_distance <- dispersal_rate * (end_year - 2011 + 1) * 1000
        msk <- range %>% st_buffer(dispersal_distance) %>% 
            st_intersection(continent) %>% st_union() %>% st_as_sf()
        
        # Clip
        cur <- cur %>% mask(msk)
        fut <- fut %>% mask(msk)
        
        # Direction change
        lyr <- ((cur >= 0) + 1) * ((fut >= 0) + 3)
        lyr <- lyr %>% terra::extend(template, fill = NA)
        
        # Summarize the changes
        lyrs <- do.call(c, lapply(c(3, 4, 6, 8), function(x){
            lyr == x
        }))
        
        dir_change <- sum(dir_change, lyrs, na.rm = TRUE)
        
        # Value change
        pos_msk <- cur >= 0
        pos_msk[pos_msk == 0] <- NA
        neg_msk <- cur < 0
        neg_msk[neg_msk == 0] <- NA
        
        chg <- (fut - cur)
        lyrs <- c(mask(chg, pos_msk), mask(chg, neg_msk))
        names(lyrs) <- c("P", "N")
        
        lyrs <- lyrs %>% terra::extend(template, fill = NA)
        
        val_change <- sum(val_change, lyrs, na.rm = TRUE)
        
        # Count species
        cur <- cur %>% terra::extend(template, fill = NA)
        num_species[[2]] <- sum(num_species[[2]], cur >= 0, 
                                na.rm = TRUE, par = TRUE)
        num_species[[3]] <- sum(num_species[[3]], cur < 0, 
                                na.rm = TRUE, par = TRUE)
        num_species[[1]] <- sum(num_species[[1]], sum(cur >= 0, cur < 0), 
                                na.rm = TRUE, par = TRUE)
    }
    
    names(dir_change) <- c("N to N", "N to P", "P to N", "P to P")
    names(num_species) <- c("All", "P", "N")
    names(val_change) <- c("P", "N")
    
    # Save out
    fname <- file.path(
        dst_dir, sprintf("dir_change_%s_%s.tif", feature, scenario))
    writeRaster(dir_change, fname)
    
    fname <- file.path(
        dst_dir, sprintf("num_species_%s_%s.tif", feature, scenario))
    writeRaster(num_species, fname)
    
    fname <- file.path(
        dst_dir, sprintf("val_change_%s_%s.tif", feature, scenario))
    writeRaster(val_change, fname)
    
    message(sprintf(
        "Finish %s for %s %s.", length(species_use), feature, scenario))
}
