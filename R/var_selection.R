#' @title var_selection
#' @description Use Bayesian Additive Regression Trees (BARTs) to select the 
#' most relevant variables.
#' @param sp (`character`) The scientific name for the species.
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
#' range_dir <- "data/IUCN/Expert_Maps"
#' aoh_dir <- "data/IUCN_AOH_100m/Mammals"
#' dst_dir <- "data/variables/variable_list"
#' tmp_dir = "data/tmp"
#' var_selection(sp, var_path, occ_dir, range_dir, aoh_dir, dst_dir, tmp_dir)
#' }
#' 
var_selection <- function(sp = "Puma_yagouaroundi",
                          var_path = "data/variables/Env/AllEnv.tif",
                          occ_dir = "data/occurrences",
                          range_dir = "data/IUCN/Expert_Maps",
                          aoh_dir = "data/IUCN_AOH_100m/Mammals", 
                          dst_dir = "data/variables/variable_list",
                          tmp_dir = "data/tmp",
                          seed = 123,
                          iter = 30){
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
    
    fname <- file.path(dst_dir, sprintf("%s_allruns.csv", sp))
    if (!file.exists(fname)){
        # presences
        occ <- read.csv(
            file.path(occ_dir, "CSVs", sprintf("%s.csv", sp)))
        occ <- st_as_sf(occ, coords = c("lon", "lat"), crs = 4326) %>% 
            st_transform(crs(vars))
        
        # Remove duplicates
        occ <- terra::extract(vars, occ, xy = TRUE, ID = FALSE)
        occ <- occ %>% unique() %>% na.omit()
        
        if (nrow(occ) == 0) stop("No occurrence left.")
        
        occ <- st_as_sf(occ, coords = c("x", "y"), crs = crs(vars))
        
        # Step 1: Detect the most contributing variables
        # Generate random background points for variable checking
        range <- st_read(file.path(range_dir, sprintf("%s.geojson", sp))) %>% 
            st_transform(crs(vars))
        occ_buf <- make_domain(occ, range); rm(range)
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
        occ_buf_tmp <- mask(
            occ_buf, rast(file.path(tmp_dir, sprintf("aoh_rsp_%s.tif", sp))), 
            inverse = TRUE)
        
        if (global(occ_buf_tmp, fun = "notNA")[[1]] == 0){
            occ_buf <- mask(
                buffer(occ_buf, res(occ_buf)[1] * 2), 
                rast(file.path(tmp_dir, sprintf("aoh_rsp_%s.tif", sp))), 
                inverse = TRUE)}
        
        ## Remove the temporary files
        file.remove(file.path(tmp_dir, sprintf("aoh_%s.tif", sp)))
        file.remove(file.path(tmp_dir, sprintf("aoh_rsp_%s.tif", sp)))
        
        # Remove presence
        occ_exist <- rasterize(occ, occ_buf)
        occ_buf <- mask(occ_buf, occ_exist, inverse = TRUE)
        vars_buf <- mask(crop(vars, occ_buf), occ_buf)
        
        # Make a series of seeds to use
        set.seed(seed)
        seeds <- sample(1:100, iter)
        
        vars_voted <- do.call(rbind, lapply(seeds, function(seed){
            # If too many, sample 1000 each time to use
            if (nrow(occ) > 3000){
                set.seed(seed)
                occ_to_use <- occ %>% sample_n(size = 1000)
            } else occ_to_use <- occ
            
            # Environmental subsampling occurrences
            occ_to_use <- st_coordinates(occ_to_use) %>% 
                cbind(st_drop_geometry(occ_to_use)) %>% 
                rename(x = X, y = Y)
            
            # Sample background
            set.seed(seed)
            bg <- spatSample(
                vars_buf, 
                min(nrow(occ_to_use), global(occ_buf, fun = "notNA")[[1]]), 
                na.rm = TRUE, xy = TRUE)
            
            # Only if the deduction won't change the number significance
            if (nrow(occ_to_use) >= 20){
                dat <- varela_sample_occ(occ_to_use, 100)
                if (nrow(dat) >= 20) occ_to_use <- dat; rm(dat)
                
                dat <- varela_sample_occ(bg, 100)
                if (nrow(dat) >= 20) bg <- dat; rm(dat)
            }
            
            occ_to_use <- occ_to_use %>% select(-c(x, y))
            bg <- bg %>% select(-c(x, y))
            dat <- rbind(occ_to_use %>% mutate(presence = 1),
                         bg %>% mutate(presence = 0))
            
            # Start to diagnose the variables using BARTs
            set.seed(seed)
            vars_selected <- var_diagnose(
                dat[, setdiff(names(dat), c("presence", "x", "y"))], 
                dat[, "presence"], n.trees = 10, iter = 50)
            
            data.frame(var = vars_selected, seed = seed)
        }))
        
        write.csv(vars_voted, fname, row.names = FALSE)
    } else vars_voted <- read.csv(fname)
    
    ## Rearrange variables based on the voting results
    vars_selected <- vars_voted %>% group_by(var) %>% 
        summarise(n = n()) %>% 
        filter(n >= floor(iter * 0.5)) %>% # 50% agree
        arrange(-n) %>% pull(var)
    
    if (length(vars_selected) == 0){
        vars_selected <- vars_voted %>% group_by(var) %>% 
            summarise(n = n()) %>% 
            filter(n == max(n)) %>% 
            pull(var)
    }
    
    # Always consider bio1 and bio12
    vars_selected <- unique(c(c("bio1", "bio12"), vars_selected))
    
    fname <- file.path(dst_dir, sprintf("%s.csv", sp))
    if (!file.exists(fname)){
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
        vars_candicates <- data.frame(var = vars_selected[sq], 
                                      var_uncorrelated = vars_thin[sq])
        
        
        write.csv(vars_candicates, fname, row.names = FALSE)
    }
}
