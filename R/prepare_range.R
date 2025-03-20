#' @title prepare_range.
#' @description Prepare the IUCN range for species.
#' @param range_path (`character`) The path of the IUCN range file.
#' @param occ_path (`character`) The path of cleaned occurrences.
#' @param sp_catalog_path (`character`) The path of species names catalog to save.
#' @param occ_dir (`character`) The directory to save out occurrences.
#' @param range_dir (`character`) The directory to save out the ranges.
#' @return No return value. A file of occurrences and range for each species
#' are save out to defined folders.
#' 
#' @details
#' Make sure setting IUCN Redlist API token correctly before any use. See
#' document of package `rredlist` for how to set up the token.
#' 
#' @import checkmate
#' @importFrom rgbif name_backbone_checklist occ_download occ_download_wait 
#' occ_download_get pred_in pred pred_not
#' @importFrom dplyr pull select filter mutate %>%
#' @importFrom here here
#' @importFrom httr GET
#' @importFrom rredlist rl_species
#' @importFrom purrr walk
#' @export
#' @examples
#' \donttest{
#' range_path <- here("data/IUCN/MAMMALS/MAMMALS.shp")
#' occ_path <- here("data/occurrences/occurrences_clean_0018860-240906103802322.csv")
#' sp_catalog_path <- here("data/occurrences/species_catalog_0018860-240906103802322.csv")
#' occ_dir <- here("data/occurrences/CSVs")
#' range_dir <- here("data/IUCN/Expert_Maps")
#' prepare_range(range_path, occ_path, sp_catalog_path, occ_dir, range_dir)
#' }
#' 
prepare_range <- function(range_path,
                          occ_path,
                          sp_catalog_path,
                          occ_dir,
                          range_dir){
    
    # Check the input
    assert_class(range_path, "character")
    assert_class(occ_path, "character")
    assert_class(sp_catalog_path, "character")
    assert_class(occ_dir, "character")
    assert_class(range_dir, "character")
    
    # Check dirs
    if(!file.exists(range_path)){
        stop("No range file found.")}
    
    if(!file.exists(occ_path)){
        stop("No occurrence file found.")}
    
    if(!file.exists(sp_catalog_path)){
        stop("No species catalog file found.")}
    
    if (!dir.exists(occ_dir)){
        dir.create(occ_dir, recursive = TRUE)}
    
    if (!dir.exists(range_dir)){
        dir.create(range_dir, recursive = TRUE)}
    
    # Load species list
    sp_records <- read.csv(occ_path)
    numbers <- sp_records %>% group_by(species) %>% summarise(n = n())
    numbers <- numbers %>% filter(n >= 15)
    sp_records <- sp_records %>% filter(species %in% unique(numbers$species))
    
    species_catalog <- read.csv(sp_catalog_path) %>% 
        filter(species %in% unique(sp_records$species))
    
    # Read the list of species from the downloaded shapefile file
    ranges <- st_read(range_path) %>% 
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
    
    # Get the synonyms with ranges as reference, ugh, such a big headache
    name_catalog <- lapply(unique(ranges$sci_name), function(nm){
        genus <- strsplit(nm, " ")[[1]][1]
        species <- strsplit(nm, " ")[[1]][2]
        
        tryCatch(nm_sd <- rl_species(genus, species),
                 error = function(e){
                     # With delays for the API
                     Sys.sleep(120)
                     tryCatch(nm_sd <- rl_species(genus, species),
                              error = function(e){
                                  # With delays for the API
                                  Sys.sleep(240)
                                  nm_sd <- rl_species(genus, species)
                              })
                 })
        synonyms <- nm_sd$taxon$synonyms
        
        if (!is.null(nrow(synonyms))){
            synonyms <- unlist(lapply(1:nrow(synonyms), function(n){
                sprintf("%s %s", 
                        synonyms[n, "genus_name"], 
                        synonyms[n, "species_name"])
            }))
        } else synonyms <- NA
        
        data.frame(range_name = nm,
                   synonyms = synonyms)
    }) %>% bind_rows()
    
    # Match them together
    bnm_sciname <- lapply(unique(ranges$sci_name), function(nm){
        nm_list <- name_catalog %>% filter(range_name == nm)
        if (nm %in% species_catalog$BINOMIAL){
            bnm <- nm
        } else {
            if (!all(is.na(nm_list$synonyms))){
                bnm <- species_catalog$BINOMIAL[
                    species_catalog$BINOMIAL %in% nm_list$synonyms]
                if (length(bnm) == 0) bnm <- NA
            } else {
                bnm <- NA
            }
        }
        data.frame(BINOMIAL = bnm, sci_name = nm)
    }) %>% bind_rows()
    
    # Some missing names
    # - Tamiops macclellandii: typo, Tamiops mcclellandii
    # - Nyctimene minutus: missing synonym record in IUCN, 
    # this is the synonym for Nyctimene varius
    # Pipistrellus tenuis, no map for the recent update
    # Tadarida insignis, recent update on the map reduce the range significantly
    # Semnopithecus dussumieri: recently recognized as an invalid taxon.
    # Scotonycteris ophiodon, disappear from the recent update
    # - Triaenops rufus, missing synonym record in IUCN, 
    # this is the synonym for Triaenops menamena
    # - Eptesicus lobatus, is known as Eptesicus serotinus.
    # - Hipposideros pendlebury, missing synonym record in IUCN,
    # this is Hipposideros pendleburyi.
    # Spalax istricus, sadly possibly extinct now.
    bnm_sciname$BINOMIAL[bnm_sciname$sci_name == "Tamiops mcclellandii"] <- 
        "Tamiops macclellandii"
    bnm_sciname$BINOMIAL[bnm_sciname$sci_name == "Nyctimene varius"] <- 
        "Nyctimene minutus"
    bnm_sciname$BINOMIAL[bnm_sciname$sci_name == "Triaenops menamena"] <- 
        "Triaenops rufus"
    bnm_sciname$BINOMIAL[bnm_sciname$sci_name == "Eptesicus serotinus"] <- 
        "Eptesicus lobatus"
    bnm_sciname$BINOMIAL[bnm_sciname$sci_name == "Hipposideros pendleburyi"] <- 
        "Hipposideros pendlebury"
    bnm_sciname <- rbind(
        bnm_sciname, c("Eptesicus serotinus", "Eptesicus serotinus"))
    
    # All good
    bnm_sciname <- bnm_sciname %>% na.omit()
    
    # Save out as querying this is time consuming and full of unexpected risks
    fname <- file.path(dirname(occ_dir), "species_name_catalog.csv")
    write.csv(bnm_sciname, fname, row.names = FALSE)
    
    # Cook the range for each species
    # Clean the csv at the same time
    sf_use_s2(FALSE)
    walk(sort(unique(sp_records$species)), function(sp){
        message(sp)
        bnm <- species_catalog %>% filter(species == sp)
        range_nms <- bnm_sciname %>% filter(BINOMIAL %in% bnm$BINOMIAL)
        
        sp_range <- ranges %>% filter(sci_name %in% range_nms$sci_name) %>% 
            select(sci_name)
        
        # # In case to fix some too complex polygons
        # sp_range <- sp_range %>% ms_simplify() %>% 
        #     st_make_valid()
        
        if (nrow(sp_range) > 0){
            fname <- file.path(
                range_dir, sprintf("%s.geojson", gsub(" ", "_", sp)))
            st_write(sp_range, fname)
            
            # Subset the occurrences
            occ_sp <- sp_records %>% filter(species == sp) %>% 
                mutate(gbif_id = gbifID, sp = species, 
                       lat = decimalLatitude, lon = decimalLongitude) %>% 
                select(gbif_id, sp, lat, lon)
            
            fname <- file.path(
                occ_dir, sprintf("%s.csv", gsub(" ", "_", sp)))
            write.csv(occ_sp, fname, row.names = FALSE)
        }
    })
}
