#' @title query_aoh
#' @description Query AOHs of mammals or birds from Lumbierres et al. (2022).
#' Lumbierres, Maria; Dahal, Prabhat Raj; Soria, Carmen D. et al. (2022). 
#' Area of Habitat maps for the world's terrestrial birds and mammals [Dataset].
#'  Dryad. https://doi.org/10.5061/dryad.02v6wwq48.
#' #' @param taxon (`character`) The taxon to process. Must be either Mammals or 
#' Birds.
#' @param aoh_dir (`character`) The destination directory to save out the 
#' queried AOHs.
#' 
#' @importFrom checkmate assert_choice assert_class
#' @importFrom dplyr pull select filter mutate %>%
#' @importFrom here here
#' @importFrom XML htmlParse xpathSApply xmlGetAttr
#' @importFrom httr GET
#' @importFrom stringr str_extract str_detect
#' @importFrom pbapply pbwalk
#' @export
#' @examples
#' \donttest{
#' query_aoh("Mammals")
#' }
#' 
query_aoh <- function(taxon = "Mammals",
                      aoh_dir = here("data/IUCN_AOH_100m"),
                      verbose = FALSE){

    # Check the input
    assert_choice(
        taxon, c('Mammals', 'Birds'))
    assert_class(aoh_dir, "character")
    
    # Define directory for occurrences
    aoh_dir <- file.path(aoh_dir, taxon)
    if(!dir.exists(aoh_dir)) {
        dir.create(aoh_dir, recursive = TRUE)
    } else {
        stop("AOHs have been downloaded.")
    }
    
    if (verbose) message("Get the list of files from Dyrad.")
    # This dataset is too large to download directly
    # so need to download each file one by one
    # Parse the list of files on Dryad
    link <- "https://doi.org/10.5061/dryad.02v6wwq48"
    resource <- GET(link)
    
    # Get the ids for the files
    parse <- htmlParse(resource)
    links <- xpathSApply(parse, path = "//a", xmlGetAttr, "href") 
    ids <- links[str_detect(links, "/stash/downloads/file_stream")]
    ids <- str_extract(ids, "[0-9]+")
    
    # Package rdryad
    ## Functions to read API results
    .dGETasync <- function(urls, query = list()) {
        if (length(urls) > 30){
            urls <- split(urls, ceiling(seq_along(urls) / 30))
        } else{
            urls <- list(urls)
        }
        do.call(c, lapply(urls, function(url_list){
            con <- crul::Async$new(
                urls = url_list,
                opts = list(),
                headers = list(
                    `Content-Type` = "application/json",
                    `Accept` = "application/json"))
            res <- con$get(query = query)
            rst <- lapply(res, function(z) z$parse("UTF-8"))
            Sys.sleep(60)
            rst
        }))
    }
    
    .v2_parse <- function(x) {
        jsonlite::fromJSON(x, flatten = TRUE)
    }
    
    .parse_each <- function(x, dois) {
        stats::setNames(lapply(x, .v2_parse), dois)
    }
    
    # Get the urls for the files
    urls <- sapply(ids, function(id){
        sprintf("https://datadryad.org/api/v2/files/%s", id)
    })
    
    tmp <- .dGETasync(urls)
    files <- .parse_each(tmp, ids)
    
    # Filter files for taxon
    files <- lapply(files, function(file_obj){
        if (str_detect(file_obj$path, taxon) & 
            str_detect(file_obj$path, ".zip") &
            !str_detect(file_obj$path, "Richness")){
            file_obj
        } else NULL
    }); files <- files[!sapply(files, is.null)]
    
    if (verbose) message("Start to download.")
    # Download the filter files
    pbwalk(files, function(file_obj){
        fname <- file_obj$path
        path <- file_obj$`_links`$`stash:download`$href
        url <- file.path("https://datadryad.org", path)
        con <- crul::HttpClient$new(url = url)
        res <- con$get()
        
        file <- file.path(aoh_dir, fname)
        file_con <- file(file, "wb")
        writeBin(res$content, con = file_con)
        on.exit(close(file_con))
        flush(file_con); rm(con, file_con) # completely close the con
        
        # unzip the file
        unzip(file, exdir = aoh_dir)
        
        # Remove the zip file
        file.remove(file)
    })
}
