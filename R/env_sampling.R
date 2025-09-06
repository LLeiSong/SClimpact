#' @title env_sampling
#' @description Environmental subsampling occurrences, sample background points
#' and environmental subsampling the points as well.
#' @param src_dir (`character`) The directory of original CSVs of occurrences.
#' @param var_dir (`character`) The directory for environmental variables.
#' @param range_dir (`character`) The directory for IUCN ranges.
#' @param occ_dir (`character`) The directory to save thinned occurrences.
#' @param bg_dir (`character`) The directory to save background samples.
#' @param maxBGFit (`integer`) The maximum number of background samples. 
#' @param seed (`integer`) The seed for randomization. 
#' @param ncores (`integer`) The number of cores to use for parallel. 
#' @return CSVs of environmental thinned occurrences are saved to `occ_dir`. 
#' CSVs of background samples are saved to `bg_dir`.
#' 
#' @import checkmate
#' @import pbmcapply
#' @import terra
#' @import sf
#' @import dplyr
#' @export
#' @examples
#' \donttest{
#' src_dir <- "data/occurrences/CSVs"
#' var_dir <- "data/variables"
#' range_dir <- "data/IUCN/Expert_Maps"
#' region_dir <- "data/terr-ecoregions-TNC"
#' occ_dir <- "data/occurrences/CSVs_thin"
#' bg_dir <- "data/occurrences/bg"
#' env_sampling(src_dir, var_dir, range_dir, occ_dir, bg_dir)
#' }
#' 
env_sampling <- function(species = NULL,
                         src_dir = here("data/occurrences/CSVs"),
                         var_dir = here("data/variables"),
                         range_dir = here("data/IUCN/Expert_Maps"),
                         region_dir = here("data/terr-ecoregions-TNC"),
                         occ_dir = here("data/occurrences/CSVs_thin"),
                         bg_dir = here("data/occurrences/bg"),
                         maxBGFit = 60000,
                         seed = 123,
                         ncores = detectCores() - 2){
    # Check inputs
    assert_class(src_dir, "character")
    assert_class(var_dir, "character")
    assert_class(range_dir, "character")
    assert_class(region_dir, "character")
    assert_class(occ_dir, "character")
    assert_class(bg_dir, "character")
    
    # Create the dst_dir
    if (!dir.exists(occ_dir)) dir.create(occ_dir)
    if (!dir.exists(bg_dir)) dir.create(bg_dir)
    # terraOptions(tempdir = "/home/lsong/temp")
    
    # Pull the list of species
    if (is.null(species)){
        species <- list.files(src_dir)
        species <- gsub(".csv", "", species)
        
        species2 <- list.files(file.path(var_dir, "variable_list"))
        species2 <- gsub(".csv", "", species2)
        
        species <- intersect(species, species2)
    }
    
    message(sprintf("%s species to process.", length(species)))
    
    # Start to thin
    pbmclapply(species, function(sp){
        occ_path <- file.path(occ_dir, sprintf("%s.csv", sp))
        bg_path <- file.path(bg_dir, sprintf("%s.csv", sp))
        
        if (!(file.exists(occ_path) & file.exists(bg_path))){
            # Load raster
            vars <- rast(file.path(var_dir, "Env/AllEnv.tif"))
            var_list <- read.csv(
                file.path(var_dir, "variable_list", sprintf("%s.csv", sp)))
            if (length(na.omit(var_list$var_uncorrelated)) == 1){
                var_list <- na.omit(var_list$var)
            } else var_list <- na.omit(var_list$var_uncorrelated)
            vars <- subset(vars, var_list)
            
            # presences
            occ <- read.csv(
                file.path(src_dir, sprintf("%s.csv", sp)))
            occ <- st_as_sf(occ, coords = c("lon", "lat"), crs = 4326) %>% 
                st_transform(crs(vars))
            
            # Remove duplicates
            occ <- terra::extract(vars, occ, xy = TRUE, ID = FALSE)
            occ <- occ %>% unique() %>% na.omit()
            
            if (nrow(occ) == 0) stop("No occurrence left.")
            
            occ <- st_as_sf(occ, coords = c("x", "y"), crs = crs(vars))
            
            # Generate random background points from the occupied ecoregions
            ## The justification is background samples should cover the whole
            ## area, at least the area that can be reached by the species
            ## So here the area for background points is the ecoregions
            ## overlapped by range and the neighbors.
            range <- st_read(
                file.path(range_dir, sprintf("%s.geojson", sp)), 
                quiet = TRUE) %>% st_transform(crs(vars))
            ecoregions <- st_read(
                file.path(region_dir, "tnc_terr_ecoregions.shp"), 
                quiet = TRUE) %>% st_transform(crs(vars))
            
            # Get the touched ecoregions and the neighbors
            regions_touch <- ecoregions %>% 
                slice(unique(unlist(st_intersects(range, ecoregions)))) %>% 
                st_union() %>% st_buffer(1)
            regions_touch <- ecoregions %>% 
                slice(unique(unlist(st_intersects(regions_touch, ecoregions))))
            
            # Mask variables
            vars_buf <- mask(crop(vars, regions_touch), regions_touch)
            
            # get the largest extent of the region
            yext <- (st_bbox(regions_touch)[4] - st_bbox(regions_touch)[2]) / 5 # y range
            xext <- (st_bbox(regions_touch)[3] - st_bbox(regions_touch)[1]) / 5 # x range
            max_ext <- max(c(yext, xext))
            
            fnames <- list.files(src_dir, full.names = TRUE)
            occs <- lapply(fnames, function(fname){
                read.csv(fname) %>% 
                    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
                    st_transform(crs(vars))
            }) %>% bind_rows()
            occs <- occs %>% st_intersection(regions_touch %>% st_union())
            
            tgb_kde <- suppressMessages(
                spatialEco::sf.kde(
                x = occs,
                bw = round(max_ext, 2),
                ref = vars_buf[[1]],
                standardize = TRUE,
                mask = TRUE))
            
            # Get background
            set.seed(seed)
            # generate random samples from the KDE raster file
            kde_samples <- spatSample(
                tgb_kde, min(maxBGFit, global(vars_buf, fun = "notNA")[[1]]), 
                method = "weights", na.rm = TRUE, values = FALSE, cells = TRUE)
            bg <- extract(vars_buf, kde_samples$cell, xy = TRUE)
            
            # Merge variables
            occ <- st_coordinates(occ) %>% 
                cbind(st_drop_geometry(occ)) %>% 
                rename(x = X, y = Y)
            
            # Only if the deduction won't change the number significance
            if (length(var_list) > 1){
                if (nrow(occ) >= 20){
                    dat <- varela_sample_occ(occ, 100)
                    if (nrow(dat) >= 20) occ <- dat; rm(dat)
                }
                
                # Environmental subsampling bg points
                dat <- varela_sample_occ(bg, 100)
                if (nrow(dat) >= 20) bg <- dat; rm(dat)
            }
            
            # Save out
            occ <- occ %>% mutate(sp = sp, id = 1:nrow(.)) %>% 
                select(id, sp, x, y)
            write.csv(occ, occ_path, row.names = FALSE)
            bg <- bg %>% mutate(sp = sp) %>% select(sp, x, y)
            write.csv(bg, bg_path, row.names = FALSE)
        }
    }, mc.cores = ncores, ignore.interactive = TRUE)
}
