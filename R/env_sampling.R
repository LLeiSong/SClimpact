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
#' occ_dir <- "data/occurrences/CSVs_thin"
#' bg_dir <- "data/occurrences/bg"
#' env_sampling(src_dir, var_dir, range_dir, occ_dir, bg_dir)
#' }
#' 
env_sampling <- function(src_dir = here("data/occurrences/CSVs"),
                         var_dir = here("data/variables"),
                         range_dir = here("data/IUCN/Expert_Maps"),
                         occ_dir = here("data/occurrences/CSVs_thin"),
                         bg_dir = here("data/occurrences/bg"),
                         maxBGFit = 20000,
                         seed = 123){
    # Check inputs
    assert_class(src_dir, "character")
    assert_class(var_dir, "character")
    assert_class(range_dir, "character")
    assert_class(occ_dir, "character")
    assert_class(bg_dir, "character")
    
    # Create the dst_dir
    if (!dir.exists(occ_dir)) dir.create(occ_dir)
    if (!dir.exists(bg_dir)) dir.create(bg_dir)
    
    # Pull the list of species
    species <- list.files(src_dir)
    species <- gsub(".csv", "", species)
    
    # Start to thin
    pbmclapply(speices, function(sp){
        # Load raster
        vars <- rast(file.path(var_dir, "Env/AllEnv.tif"))
        var_list <- read.csv(
            file.path(var_dir, "variable_list", sprintf("%s.csv", sp)))
        vars <- subset(vars, na.omit(var_list$var_uncorrelated))
        
        # presences
        occ <- read.csv(
            file.path(src_dir, sprintf("%s.csv", sp)))
        occ <- st_as_sf(occ, coords = c("lon", "lat"), crs = 4326) %>% 
            st_transform(crs(vars))
        
        # Generate random background points for variable checking
        range <- st_read(file.path(range_dir, sprintf("%s.geojson", sp))) %>% 
            st_transform(crs(vars))
        occ_buf <- make_domain(occ, range); rm(range)
        occ_buf <- terra::rasterize(occ_buf, vars[[1]]) %>% terra::trim()
        
        # Remove presence
        occ_exist <- rasterize(occ, occ_buf)
        occ_buf <- mask(occ_buf, occ_exist, inverse = TRUE)
        vars_buf <- mask(crop(vars, occ_buf), occ_buf)
        
        # Get background
        set.seed(seed)
        bg <- spatSample(
            vars_buf, min(maxBGFit, global(vars_buf, fun = "notNA")[[1]]), 
            na.rm = TRUE, xy = TRUE)
        
        # Merge variables
        occ <- terra::extract(vars, occ, bind = TRUE)
        coords <- geom(occ)
        occ <- data.frame(x = coords[, "x"], y = coords[, "y"]) %>% 
            cbind(as.data.frame(occ)) %>% na.omit()
        
        # Only if the deduction won't change the number significance
        if (nrow(occ) >= 20){
            dat <- varela_sample_occ(occ %>% dplyr::select(-c(gbif_id, sp)), 100)
            if (nrow(dat) >= 20) {
                occ <- merge(dat, occ %>% dplyr::select(c(x, y, gbif_id, sp)), 
                             by = c("x", "y"), all.x = TRUE)
            }; rm(dat)
        }
        occ <- occ %>% select(gbif_id, sp, x, y)
        
        # Save out
        write.csv(occ, file.path(occ_dir, sprintf("%s.csv", sp)),
                  row.names = FALSE)
        
        # Environmental subsampling bg points
        dat <- varela_sample_occ(bg, 100)
        if (nrow(dat) >= 20) bg <- dat; rm(dat)
        bg <- bg %>% select(x, y) %>% mutate(sp = sp)
        
        # Save out
        write.csv(bg, file.path(bg_dir, sprintf("%s.csv", sp)),
                  row.names = FALSE)
    }, mc.cores = detectCores() - 2, ignore.interactive = TRUE)
}
