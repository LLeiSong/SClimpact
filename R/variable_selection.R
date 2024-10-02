#' @title var_selection
#' @description Use Bayesian Additive Regression Trees (BARTs) to select the 
#' most relevant variables.
#' #' @param sp (`character`) The scientific name for the species.
#' @param var_path (`character`) The path for the environmental variables.
#' @param occ_dir (`character`) The directory for occurrences.
#' @param aoh_dir (`character`) The directory for species AOH.
#' @param dst_dir (`character`) The destination directory to save out the 
#' CSVs of selected variables for each species.
#' @param tmp_dir (`character`) The directory to save temporary files. 
#' @param seed (`integer`) The seed for randomization. 
#' @return A csv of selected variables is saved to `dst_dir` for this species.
#' 
#' @details
#' Make sure GDAL library is installed on your machine to make this work.
#' 
#' @import checkmate
#' @importFrom dbarts bart
#' @importFrom dplyr pull select filter mutate %>%
#' @import terra
#' @import sf
#' @import dplyr
#' @export
#' @examples
#' \donttest{
#' sp <- "Puma_yagouaroundi"
#' var_path <- "data/variables/Env/AllEnv.tif"
#' occ_dir <- "data/occurrences"
#' aoh_dir <- "data/IUCN_AOH_100m/Mammals"
#' dst_dir <- "data/variables/variable_list"
#' tmp_dir = "data/tmp"
#' var_selection(sp, var_path, occ_dir, aoh_dir, dst_dir, 123)
#' }
#' 
var_selection <- function(sp = "Puma_yagouaroundi",
                          var_path = "data/variables/Env/AllEnv.tif",
                          occ_dir = "data/occurrences",
                          aoh_dir = "data/IUCN_AOH_100m/Mammals", 
                          dst_dir = "data/variables/variable_list",
                          tmp_dir = "data/tmp",
                          seed = 123,
                          iter = 20){
    # Check the input
    assert_class(sp, "character")
    assert_class(var_path, "character")
    assert_class(occ_dir, "character")
    assert_class(dst_dir, "character")
    assert_class(seed, "integer")
    
    # Create the dst_dir
    if (!dir.exists(dst_dir)) dir.create(dst_dir)
    
    # Define temporary folder
    if (!dir.exists(tmp_dir)) dir.create(tmp_dir)
    
    # Load raster
    vars <- rast(var_path)
    
    # presences
    occ <- read.csv(
        file.path(occ_dir, "CSVs", sprintf("%s.csv", sp)))
    occ <- st_as_sf(occ, coords = c("lon", "lat"), crs = 4326) %>% 
        st_transform(crs(vars))
    
    # Step 1: Detect the most contributing variables
    # Generate random background points for variable checking
    occ_buf <- buf_polygon(occ)
    occ_buf <- terra::rasterize(occ_buf, vars[[1]]) %>% terra::trim()
    
    # AOH masking
    ## Load Catalog
    aoh_path <- read.csv(
        file.path(occ_dir,
                  "species_catalog_0018860-240906103802322.csv")) %>% 
        filter(species == gsub("_", " ", sp))
    
    ## Load AOH
    aoh <- sprc(lapply(aoh_path$AOH_name, function(fn){
        rast(file.path(aoh_dir, sprintf("%s.tif", fn)))
    })); aoh <- merge(aoh); names(aoh) <- "aoh"
    
    ## Resample AOH using mode
    writeRaster(
        aoh, file.path(tmp_dir, sprintf("aoh_%s.tif",sp)), 
        datatype = "INT1U")
    te <- ext(occ_buf); tr <- res(occ_buf)
    options_warp <- paste0(
        '-multi -wo NUM_THREADS=ALL_CPUS ',
        '-co TILED=YES -co compress=lzw -co BIGTIFF=IF_NEEDED')
    command <- sprintf(
        paste0("gdalwarp -t_srs EPSG:%s -te %s %s %s %s",
               " -tr %s %s %s %s %s -r mode"),
        crs(occ_buf, describe = TRUE)$code,
        te[1], te[3], te[2], te[4], tr[1], tr[2], 
        options_warp, file.path(tmp_dir, sprintf("aoh_%s.tif", sp)), 
        file.path(tmp_dir, sprintf("aoh_rsp_%s.tif", sp)))
    system(command)
    
    ## Mask the buffered occurrences
    occ_buf <- mask(
        occ_buf, rast(file.path(tmp_dir, sprintf("aoh_rsp_%s.tif", sp))), 
        inverse = TRUE)
    
    ## Remove the temporary files
    file.remove(file.path(tmp_dir, sprintf("aoh_%s.tif", sp)))
    file.remove(file.path(tmp_dir, sprintf("aoh_rsp_%s.tif", sp)))
    
    # Remove presence
    occ_exist <- rasterize(occ, occ_buf)
    occ_buf <- mask(occ_buf, occ_exist, inverse = TRUE)
    vars_buf <- mask(crop(vars, occ_buf), occ_buf)
    
    set.seed(seed)
    seeds <- sample(1:100, iter)
    
    vars_voted <- do.call(rbind, lapply(seeds, function(seed){
        set.seed(seed)
        bg <- spatSample(vars_buf, nrow(occ), na.rm = TRUE, xy = TRUE)
        
        # Environmental subsampling occurrences
        occ <- terra::extract(vars, occ, bind = TRUE)
        coords <- geom(occ)
        occ <- data.frame(x = coords[, "x"], y = coords[, "y"]) %>% 
            cbind(as.data.frame(occ) %>% 
                      dplyr::select(-c(gbif_id, sp))) %>% na.omit()
        
        # Only if the deduction won't change the number significance
        if (nrow(occ) >= 20){
            dat <- varela_sample_occ(occ, 100)
            if (nrow(dat) >= 20) occ <- dat; rm(dat)
            
            dat <- varela_sample_occ(bg, 100)
            if (nrow(dat) >= 20) bg <- dat; rm(dat)
        }
        
        occ <- occ %>% select(-c(x, y))
        bg <- bg %>% select(-c(x, y))
        dat <- rbind(occ %>% mutate(presence = 1),
                     bg %>% mutate(presence = 0))
        
        # Start to diagnose the variables using BARTs
        set.seed(seed)
        vars_selected <- var_diagnose(
            dat[, setdiff(names(dat), c("presence", "x", "y"))], 
            dat[, "presence"], n.trees = 10, iter = 50)
        
        data.frame(var = vars_selected, seed = seed)
    }))
    
    ## Rearrange variables based on the voting results
    vars_selected <- vars_voted %>% group_by(var) %>% 
        summarise(n = n()) %>% 
        filter(n >= ceiling(iter * 0.5)) %>% 
        arrange(-n) %>% pull(var)
    
    if (length(vars_selected) == 0){
        vars_selected <- vars_voted %>% group_by(var) %>% 
            summarise(n = n()) %>% 
            filter(n == max(n)) %>% 
            pull(var)
    }
    
    # Step 2: Remove highly correlated variables
    if (length(vars_selected) > 3){
        # Reload presences
        occ <- read.csv(
            file.path(occ_dir, "CSVs", sprintf("%s.csv", sp)))
        occ <- st_as_sf(occ, coords = c("lon", "lat"), crs = 4326) %>% 
            st_transform(crs(vars))
        
        # Remove highly correlated variables
        vars_occ <- terra::extract(vars[[vars_selected]], occ, ID = FALSE) %>% 
            na.omit()
        vars_thin <- var_remove_cor(vars_occ, c("bio1", "bio12"))
    } else {
        vars_thin <- vars_selected
    }
    
    sq <- seq(max(length(vars_selected), length(vars_thin)))
    data.frame(var = vars_selected[sq], 
               var_uncorrelated = vars_thin[sq])
    
    fname <- file.path(dst_dir, sprintf("%s.csv", sp))
    write.csv(vars, fname, row.names = FALSE)
}
