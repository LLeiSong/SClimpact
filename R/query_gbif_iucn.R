#' @title query_gbif.
#' @description Interactive function to download the raw occurrences 
#' for all species from GBIF website. Use a faster and more elegant way to
#' request with function `occ_download` in `rgbif`. And more importantly, can
#' ensure a specific citation for the data.
#' #' @param taxon (`character`) The taxon to process. Must be either Mammals or 
#' Birds.
#' @param occ_dir (`character`) The destination directory to save out the 
#' queried occurrences.
#' @return (`occ_download_get`) The \link{occ_download_get} object from `rgbif` 
#' or `NULL` if the request fails. The `occ_download_get` can directly be passed
#' to function \link{occ_download_import} for further processing.
#' 
#' @import checkmate
#' @importFrom rgbif name_backbone_checklist occ_download occ_download_wait 
#' occ_download_get pred_in pred pred_not
#' @importFrom dplyr pull select filter mutate %>%
#' @importFrom here here
#' @importFrom httr GET
#' @export
#' @examples
#' \donttest{
#' query_gbif("Mammals")
#' }
#' 

query_gbif <- function(
        range_path = here("data/IUCN/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp"),
        occ_dir = here("data/occurrences")){
    
    # Check the input
    assert_class(occ_dir, "character")
    
    # Define directory for occurrences
    if(!dir.exists(occ_dir)) dir.create(occ_dir)
    
    # Read the list of species from the downloaaded shapefile file
    species_catalog <- st_read(range_path) %>% 
        filter(terrestria == "true") %>% # just in case
        filter(freshwater == "false") %>% 
        filter(marine == "false") %>% 
        ## Subset the interested ones
        ## presence: extant (1), Possibly extant (2 or 3)
        ##           others: Possibly Extinct (4), and Extinct (5)
        ## origin: native (1) and reintroduced (2), and Assisted Colonisation (6)
        ##         others: introduced (3), Vagrant (4), and Origin Uncertain (5)
        ## season: resident (1)
        ##         others: breeding season (2) and non-breeding season (3)
        ##         Passage(4), and Seasonal Occurrence Uncertain (5)
        ## category: exclude EW and EX
        filter(presence %in% 1:3) %>% 
        filter(origin %in% c(1:2, 6)) %>% 
        filter(seasonal == 1) %>% 
        filter(!category %in% c("EW", "EX"))
    
    # Simplify the table
    species_catalog <- species_catalog %>% st_drop_geometry() %>% 
        select(sci_name, category) %>% unique()
    
    # Get the taxon keys in GBIF
    species_catalog <- 
        species_catalog %>%
        pull(sci_name) %>% 
        name_backbone_checklist() %>% 
        cbind(species_catalog, .)
    
    species_catalog <- species_catalog %>% filter(!is.na(species))
    
    # Diagnose some results
    species_to_check <- species_catalog %>% filter(matchType != "EXACT") %>%
        select(sci_name, scientificName, canonicalName, matchType)
    
    fname <- file.path(occ_dir, "species_to_check.csv")
    write.csv(species_to_check, fname, row.names = FALSE)
    message(sprintf(
        paste0("A CSV file has been saved to %s. ",
               "Please check the file manually, ",
               "and provide values of matchType to remove."), fname))
    
    matchType_to_remove <- readline(
        "Please enter matchTypes (use comma to seperate) to remove: ")
    matchType_to_remove <- trimws(strsplit(matchType_to_remove, ",")[[1]])
    file.remove(fname); rm(fname)
    
    # Remove some problematic species
    species_catalog <- species_catalog %>% 
        # means cannot find any match
        filter(!matchType %in% matchType_to_remove)
    
    # Update meta for download:
    # 1. Use acceptedUsageKey if any for search
    # 2. Update the scientific name for that
    species_catalog <- species_catalog %>% 
        mutate(usageKey_query = ifelse(
            is.na(acceptedUsageKey), usageKey, acceptedUsageKey))
    
    # Query the download and save out the meta
    ## Add a few internal filters for something that we definitely not want.
    ## 1. Zero coordinate : Coordinates are exactly (0,0). null island
    ## 2. Country coordinate mismatch : The coordinates fall outside of the 
    ##    given country’s polygon.
    ## 3. Coordinate invalid : GBIF is unable to interpret the coordinates.
    ## 4. Coordinate out of range : The coordinates are outside of the range 
    ##    for decimal lat/lon values ((-90,90), (-180,180)).
    ## 5. Only keep PRESENT records
    ## 6. Remove FOSSIL_SPECIMEN and LIVING_SPECIMEN
    occ_meta <- occ_download(
        pred_in("taxonKey", unique(species_catalog$usageKey_query)),
        pred("occurrenceStatus", "PRESENT"), 
        pred("hasCoordinate", TRUE), pred("hasGeospatialIssue", FALSE),
        pred_not(pred_in("basisOfRecord",
                         c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
        format = "SIMPLE_CSV")
    
    # Save out the processed catalog
    write.csv(species_catalog, 
              file.path(occ_dir, sprintf("species_catalog_%s.csv", occ_meta)),
              row.names = FALSE)
    
    # Download the occurrences
    status <- occ_download_wait(occ_meta)
    if (status$status == "SUCCEEDED"){
        d <- occ_download_get(occ_meta, path = occ_dir, overwrite = TRUE)
    } else {
        warning(sprintf("The request %s.", status$status))
        d <- NULL
    }
    
    # return the meta
    return(d)
}
