#' @title filter_gbif
#' @description Apply extra meta, quality, and spatial filters to the 
#' downloaded GBIF records to increase their quality.
#' #' @param x (`occ_download_get`) The \link{occ_download_get} object 
#' from `rgbif`. The output of function \link{query_gbif} can directly pass
#' into the function as `x`. Users can also run \link{occ_download_get} 
#' function to re-build it.
#' @param occ_dir (`character`) The destination directory to save out the 
#' queried occurrences.
#' @param aoh_dir (`character`) The directory to read AOH of species for
#' spatial filtering.
#' @param cores (`integer`) The number of cores to use for parallel.
#' The default is all available cores.
#' @param devision (`logical`) To save out csv for each species or not.
#' The default is TRUE to save out.
#' @param verbose (`logical`) If TRUE reports the name of the test and 
#' the number of records flagged.
#' @return (`data.frame`) The filtered `data.frame` of occurrences.
#' 
#' @import CoordinateCleaner
#' @importFrom checkmate assert_class
#' @importFrom rgbif occ_download_import
#' @importFrom dplyr filter distinct select %>% bind_rows
#' @importFrom here here
#' @importFrom parallel detectCores mclapply
#' @importFrom pbmcapply pbmclapply
#' @importFrom sf st_as_sf
#' @importFrom terra rast sprc merge
#' @importFrom purrr walk
#' @export
#' @examples
#' \donttest{
#' d <- query_gbif("Mammals")
#' filter_gbif(d)
#' }
#' 

filter_gbif <- function(
        x,
        occ_dir = here("data/occurrences"),
        aoh_dir = here("data/IUCN_AOH_100m/Mammals"),
        cores = detectCores(),
        division = TRUE,
        verbose = TRUE){
    
    # Check inputs
    assert_class(x, "occ_download_get")
    assert_class(occ_dir, "character")
    assert_class(aoh_dir, "character")
    
    # Define directory for occurrences
    if(!dir.exists(occ_dir)) stop("No directory to read occurrences.")
    if(!dir.exists(aoh_dir)) stop("No directory to read AOH.")
    
    # Create the quality flag meta for the occurrences
    quality_flag <- data.frame(
        flag = c("Raw", "blank (gbif)", "age (gbif)", "precision (gbif)", 
                 "validity (CC)", "duplicate (CC)", "equal (CC)", "zero (CC)", 
                 "capital (CC)", "centroid (CC)", "gbif (CC)", 
                 "institution (CC)", "insufficient", "AOH", "insufficient"),
        clean = NA, species = NA)
    
    # Import GBIF occurrences
    occ <- occ_download_import(x) # 25830904
    
    # Update flag meta
    quality_flag[1, 2] <- nrow(occ)
    quality_flag[1, 3] <- length(unique(occ$speciesKey))
    
    # Step 1: Basic meta filter (remove obvious noises/errors)
    ## 1.1 Using GBIF meta-data
    ### remove records without species info
    occ <- occ %>% filter(!is.na(species))
    
    ### Update flag meta
    quality_flag[2, 2] <- nrow(occ)
    quality_flag[2, 3] <- length(unique(occ$speciesKey))
    
    ### before the end of world war II
    occ <- occ %>% filter(year >= 1945 | is.na(year))
    
    ### Update flag meta
    quality_flag[3, 2] <- nrow(occ)
    quality_flag[3, 3] <- length(unique(occ$speciesKey))
    
    ### Remove low precise records
    occ <- occ %>% 
        filter(coordinatePrecision < 0.01 | is.na(coordinatePrecision)) %>% 
        filter(coordinateUncertaintyInMeters / 1000 <= 100 | 
                   is.na(coordinateUncertaintyInMeters)) %>%
        filter(!coordinateUncertaintyInMeters %in% c(301, 3036, 999, 9999))
    
    ### Update flag meta
    quality_flag[4, 2] <- nrow(occ)
    quality_flag[4, 3] <- length(unique(occ$speciesKey))
    
    ## 1.2 Use CoordinateCleaner to flag problematic records
    ## reorder the checks to increase efficiency.
    ### validity
    occ <- occ %>% cc_val(verbose = verbose)
    
    ### Update flag meta
    quality_flag[5, 2] <- nrow(occ)
    quality_flag[5, 3] <- length(unique(occ$speciesKey))
    
    ### duplicates
    occ <- occ %>% cc_dupl(verbose = verbose)
    
    ### Update flag meta
    quality_flag[6, 2] <- nrow(occ)
    quality_flag[6, 3] <- length(unique(occ$speciesKey))
    
    ### Other flags in CoordinateCleaner
    ### Do not use sea and country filter because: first, it is too slow. 
    ### second, the downstream AOH filter will do this as well.
    ### Use conservative radius (100m) for buffering filters. It is enough to 
    ### tolerate the mismatch with coordinates but also keep the real records.
    occ <- clean_coordinates(
        x = occ, capitals_rad = 100, centroids_rad = 100, aohi_rad = 100,
        tests = c("capitals", "centroids", "equal", 
                  "gbif", "institutions", "zeros"),
        verbose = verbose)
    
    ### Update flag meta in order
    quality_flag[7, 2] <- sum(occ$.equ)
    quality_flag[7, 3] <- length(unique(occ$speciesKey[occ$.equ]))
    quality_flag[8, 2] <- sum(occ$.zer & occ$.equ)
    quality_flag[8, 3] <- length(unique(occ$speciesKey[occ$.zer & occ$.equ]))
    quality_flag[9, 2] <- sum(occ$.cap & occ$.zer & occ$.equ)
    quality_flag[9, 3] <- length(unique(occ$speciesKey[
        occ$.cap & occ$.zer & occ$.equ]))
    quality_flag[10, 2] <- sum(occ$.cen & occ$.cap & occ$.zer & occ$.equ)
    quality_flag[10, 3] <- length(unique(occ$speciesKey[
        occ$.cen & occ$.cap & occ$.zer & occ$.equ]))
    quality_flag[11, 2] <- sum(occ$.gbf & occ$.cen & occ$.cap & 
                                   occ$.zer & occ$.equ)
    quality_flag[11, 3] <- length(unique(occ$speciesKey[
        occ$.gbf & occ$.cen & occ$.cap & occ$.zer & occ$.equ]))
    quality_flag[12, 2] <- sum(occ$.inst & occ$.gbf & occ$.cen & 
                                   occ$.cap & occ$.zer & occ$.equ)
    quality_flag[12, 3] <- length(unique(occ$speciesKey[
        occ$.inst & occ$.gbf & occ$.cen & occ$.cap & occ$.zer & occ$.equ]))
    
    # Clean the occ
    occ <- occ[occ$.summary, ] %>% 
        select(-names(occ)[grepl("\\.", names(occ))])
    
    ### Step 1.3: Remove species with less than 4 records
    species_to_move <- occ %>% group_by(species) %>% summarise(n = n()) %>% 
        filter(n < 4) %>% pull(species)
    occ <- occ %>% filter(!species %in% species_to_move)
    
    #### Update flag meta
    quality_flag[13, 2] <- nrow(occ)
    quality_flag[13, 3] <- length(unique(occ$speciesKey))
    
    # Save cleaned occ
    fname <- file.path(
        occ_dir, sprintf("occurrences_noaoh_%s", 
                         gsub(".zip", ".csv", basename(x))))
    write.csv(occ, fname, row.names = FALSE)
    
    # Step 2: AOH filter (remove spatial outliers)
    ## We use a more precise AOH instead of IUCN ranges.
    ## Because some species in AOH may be the same species in GBIF. So make
    ## the catalog for the transition with the reference of GBIF because we
    ## use the occurrence records for modeling.
    ### Prepare catalog
    fname <- file.path(
        occ_dir, sprintf("species_catalog_%s", 
                         gsub(".zip", ".csv", basename(x))))
    catalog <- read.csv(fname) %>% 
        filter(species %in% unique(occ$species))
    
    ### Filtering
    occ <- pbmclapply(unique(occ$species), function(sp){
        occ_tx <- occ %>% filter(species == sp)
        # Note that this could have multiple
        aoh_tx <- catalog %>% filter(species == sp)
        
        # Load AOH
        aoh <- sprc(lapply(aoh_tx$AOH_name, function(fn){
            rast(file.path(aoh_dir, sprintf("%s.tif", fn)))
        })); aoh <- merge(aoh); names(aoh) <- "range"
        
        # They have the same coordinate system, so it is good
        # Filter the occurrences
        pts <- occ_tx %>% select(decimalLatitude, decimalLongitude, gbifID) %>% 
            st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
                     crs = 4326)
        pts <- extract(aoh, pts, bind = TRUE)
        ids <- pts %>% st_as_sf() %>% filter(!is.na(range)) %>% pull(gbifID)
        
        occ_tx %>% filter(gbifID %in% ids)
    }, mc.cores = cores, ignore.interactive = TRUE)
    occ <- occ %>% bind_rows()
    
    ### Update flag meta
    quality_flag[14, 2] <- nrow(occ)
    quality_flag[14, 3] <- length(unique(occ$speciesKey))
    
    # Save cleaned occ
    fname <- file.path(
        occ_dir, sprintf("occurrences_afteraoh_%s", 
                         gsub(".zip", ".csv", basename(x))))
    write.csv(occ, fname, row.names = FALSE)
    
    # Step 3: Remove species with less than 4 records
    species_to_move <- occ %>% group_by(species) %>% summarise(n = n()) %>% 
        filter(n < 4) %>% pull(species)
    occ <- occ %>% filter(!species %in% species_to_move)
    
    ## Update flag meta
    quality_flag[15, 2] <- nrow(occ)
    quality_flag[15, 3] <- length(unique(occ$speciesKey))
    
    # Step 4: Remove environmental outliers and spatial thinning will be
    # applied within the modeling workflow
    
    # Save out
    fname <- file.path(
        occ_dir, sprintf("occurrences_clean_%s", 
                         gsub(".zip", ".csv", basename(x))))
    write.csv(occ, fname, row.names = FALSE)
    
    fname <- file.path(
        occ_dir, sprintf("quality_flag_%s", 
                         gsub(".zip", ".csv", basename(x))))
    write.csv(quality_flag, fname, row.names = FALSE)
    
    if (verbose){
        message(sprintf("%.1f%s of the occurrences are selected.", 
                        quality_flag[14, 2] / quality_flag[1, 2] * 100, "%"))
        message(sprintf("The cleaned data is saved to %s.", fname))
    }
    
    if (division){
        csv_path <- file.path(occ_dir, "CSVs")
        if (!dir.exists(csv_path)) dir.create(csv_path)
        
        walk(sort(unique(occ$species)), function(sp){
            occ_sp <- occ %>% filter(species == sp) %>% 
                mutate(gbif_id = gbifID, sp = species, 
                       lat = decimalLatitude, lon = decimalLongitude) %>% 
                select(gbif_id, sp, lat, lon)
            
            fname <- file.path(
                csv_path, sprintf("%s.csv", gsub(" ", "_", sp)))
            write.csv(occ_sp, fname, row.names = FALSE)
        })
        
        if (verbose){
            message(sprintf(paste0("The cleaned data is separated ",
                                   "by species and saved to %s."), csv_path))
        }
    }
}
