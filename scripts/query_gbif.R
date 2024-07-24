#' @title Query occurrence from different sources for vertebrates.
#' @description Take the scientific name to download occurrences from different
#' sources using package \code{\link{spocc}}.
#' @param sci_namne (`character`) The scientific name.
#' @param start_year (`character`) The start year for the query.
#' @param dst_dir (`character`) The destination directory for the data.
#' @return The occurrence records.
#'
#' @seealso
#' \code{\link{occ_search}}, \code{\link{st_as_sf}}
#'
#' @details
#' Convert the terrestrial continents into hex grids. Note that the map will be
#' shifted to keep Russia together. So it will be problematic if the user 
#' directly transform to geographic coordinate system.
#'
#' @importFrom rgbif  occ_search
#' @importFrom scrubr dframe coord_impossible coord_incomplete coord_unlikely
#' @importFrom sf st_as_sf st_write
#' @export
#' @examples
#' # Query occurrence for species Lerista ameles from GBIF.
#' occ <- query_gbif(sci_name = "Lerista ameles")

query_gbif <- function(sci_namne, 
                       start_year = 1981,
                       dst_dir = "."){
    ## Set the time interval for querying on GBIF
    start_year <- 1981
    year <- sprintf('%s,%s',  start_year, year(Sys.Date()))
    
    # Search
    tryCatch({
        occ <- occ_search(
            scientificName = sci_namne,
            # continent = "south_america",
            hasCoordinate = TRUE,
            limit = 100000, # maximum of the gbif API
            year = year,
            hasGeospatialIssue = FALSE)
        
        if (!is.null(occ$data)){
            # Clean the dataset
            occ <- dframe(occ$data) %>%
                coord_impossible() %>%
                coord_incomplete() %>%
                coord_unlikely()
            
            # Convert to sf
            occ <- st_as_sf(
                occ[, c("decimalLongitude", "decimalLatitude", "datasetKey", "scientificName")], 
                coords = c("decimalLongitude", "decimalLatitude"),
                crs = 4326)
            rm(start_year, year)
            
            fn <- file.path(dst_dir, sprintf("%s.geojson", gsub(" ", "_", sci_namne)))
            st_write(occ, fn)
        }
    }, error = function(e) e)
}

# # Application
# # Set directories
# occ_dir <- file.path("/home/lsong/data/SCImpact", 'gbif')
# sdm_dir <- file.path("/home/lsong/data/SCImpact", "sdm")
# for (dir_to in c(occ_dir, sdm_dir)){
#     if (!dir.exists(dir_to)) dir.create(dir_to)}
# # Start
# catalog <- read.csv("/home/lsong/data/catalog_species.csv")
# species_done <- list.files("/home/lsong/data/SCImpact/gbif")
# species_done <- gsub(".geojson", "", species_done)
# species_done <- gsub("_", " ", species_done)
# species <- setdiff(unique(catalog$sci_name), species_done)
# 
# for (sci_name in species){
#     query_gbif(sci_name)
# }
