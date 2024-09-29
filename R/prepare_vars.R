#' @title prepare_vars
#' @description Prepare the environmental variables of baseline and future scenarios.
#' #' @param clim_dir (`character`) The directory for climate variables.
#' @param lc_dir (`character`) The directory for land cover layers.
#' @param dst_dir (`character`) The destination directory to save out the 
#' processed variables.
#' @param verbose (`logical`) If TRUE reports the progress.
#' @return No return. All files will be saved out to the dst_dir.
#' 
#' @import checkmate
#' @importFrom rnaturalearth ne_countries
#' @importFrom dplyr %>% filter
#' @importFrom stringr str_detect str_extract
#' @importFrom terra rast trim project mask
#' @export
#' @examples
#' \donttest{
#' prepare_vars(clim_dir, lc_dir, dst_dir, FALSE)
#' }
#' 
prepare_vars <- function(clim_dir,
                         lc_dir,
                         dst_dir,
                         verbose = TRUE){
    orig_setting <- terraOptions()
    terraOptions(memfrac = 0.8, tempdir = "data/temp")
    if (!dir.exists("data/temp")){
        dir.create("data/temp")}
    
    # clim_dir <- '/Volumes/gazelle/CHELSA_V2/GLOBAL/climatologies'
    # lc_dir <- "/Volumes/gazelle/future_land_projection_chen/Global\
    #  7-land-types\ LULC\ projection\ dataset\ under\ SSPs-RCPs"
    
    # Check inputs
    assert_class(clim_dir, "character")
    assert_class(lc_dir, "character")
    assert_class(dst_dir, "character")
    
    # Further check the existence
    if (!dir.exists(clim_dir)) stop("Wrong directory for climatic factors.")
    if (!dir.exists(lc_dir)) stop("Wrong directory for land cover factors.")
    if (!dir.exists(dst_dir)) dir.create(dst_dir, recursive = TRUE)
    
    # Create sub directories
    env_dir <- file.path(dst_dir, "Env")
    otherenv_dir <- file.path(dst_dir, "OtherEnv")
    
    for (dr in c(env_dir, otherenv_dir)){
        if (!dir.exists(dr)) dir.create(dr)}
    
    # Define the template for global
    template <- rast(xmin = -20037508.34, xmax = 20037508.34, 
                     ymin = -20048966.1, ymax = 20048966.1,
                     crs = "EPSG:3857", resolution = 10000)

    # Download the terrestrial boundary
    bry <- ne_countries(scale = 50, type = "countries", returnclass = "sf") %>% 
        filter(continent != "Antarctica") %>% # remove Antarctica
        st_transform(crs(template))
    
    # Define meta parameters
    years <- c("1981-2010", "2011-2040", "2041-2070", "2071-2100")
    mods <- c("GFDL-ESM4", "UKESM1-0-LL", "MPI-ESM1-2-HR", 
              "IPSL-CM6A-LR", "MRI-ESM2-0")
    ssps <- c("ssp126", 'ssp370', "ssp585")
    vars <- sprintf("bio%s", 1:19)
    
    for (yr in years){
        # base condition
        if (yr == "1981-2010"){
            if (verbose) message("Start to process Env.")
            # Climate
            data_dir <- file.path(
                clim_dir, yr, "bio", 
                sprintf("CHELSA_%s_%s_V.2.1.tif", vars, yr))
            clim <- rast(file.path(data_dir))
            names(clim) <- vars
            clim <- project(clim, template, 
                            method = "bilinear", threads = TRUE)
            
            # Land cover
            lc <- rast(file.path(lc_dir, "global_LULC_2015.tif"))
            
            # Transit land cover
            # Select Forest (forest and shrub), Grassland, Cropland, Urban
            #        deforestation, land degradation, human impact (ag + urban)
            types <- c(2, 3, 5)
            lc <- do.call(c, lapply(types, function(ty){
                if (ty == 5){
                    lc == 5 | lc == 6
                } else lc == ty
            })); names(lc) <- c("forest", "grassland", "human_impact")
            
            area <- !is.na(lc[[1]])
            area <- project(area, template, method = "sum", threads = TRUE)
            lc <- project(lc, template, method = "sum", threads = TRUE)
            lc <- lc / area
            
            # Put together
            clim <- c(clim, lc)
            clim <- mask(crop(clim, bry), bry)
            clim <- trim(clim)
            
            # Save out
            fname <- file.path(env_dir, "AllEnv.tif")
            writeRaster(clim, fname)
            write.table(names(clim), file.path(env_dir, "layerNames.csv"),
                        row.names = FALSE, col.names = FALSE)
            
            # Collect garbage
            tmpFiles(current = TRUE, remove = TRUE); gc()
            
        # Future conditions
        } else {
            for (mod in mods){
                for (ssp in ssps){
                    if (verbose) {
                        msg <- sprintf(
                            paste0("Start to process OtherEnv",
                                   " for year: %s, model: %s, and SSP: %s."),
                            yr, mod, ssp)
                        message(msg)}
                    # Climate
                    data_dir <- file.path(
                        clim_dir, yr, mod, ssp, "bio", 
                        sprintf("CHELSA_%s_%s_%s_%s_V.2.1.tif", 
                                vars, yr, tolower(mod), ssp))
                    clim <- rast(file.path(data_dir))
                    names(clim) <- vars
                    clim <- project(clim, template, 
                                    method = "bilinear", threads = TRUE)
                    
                    # Land cover
                    yr_lc <- strsplit(yr, "-")[[1]][2]
                    ssp_lc <- sprintf(
                        "SSP%s_RCP%s", 
                        as.integer(str_extract(ssp, "[0-9]+")) %/% 100,
                        as.integer(str_extract(ssp, "[0-9]+")) %% 100)
                    
                    lc <- rast(file.path(
                        lc_dir, ssp_lc, 
                        sprintf("global_%s_%s.tif", ssp_lc, yr_lc)))
                    
                    # Transit land cover
                    # Select Forest, Grassland, Cropland, Urban
                    # deforestation, land degradation, human impact (ag + urban)
                    types <- c(2, 3, 5)
                    lc <- do.call(c, lapply(types, function(ty){
                        if (ty == 5){
                            lc == 5 | lc == 6
                        } else lc == ty
                    })); names(lc) <- c("forest", "grassland", "human_impact")
                    
                    area <- !is.na(lc[[1]])
                    area <- project(area, template, 
                                    method = "sum", threads = TRUE)
                    lc <- project(lc, template, 
                                  method = "sum", threads = TRUE)
                    lc <- lc / area
                    
                    # Put together
                    clim <- c(clim, lc)
                    clim <- mask(crop(clim, bry), bry)
                    clim <- trim(clim)
                    
                    # Save out
                    fname <- file.path(
                        otherenv_dir, 
                        sprintf("%s_%s_%s.tif", 
                                str_extract(mod, "^[A-Z|a-z]+"), ssp, yr))
                    writeRaster(clim, fname)
                    
                    # Collect garbage
                    tmpFiles(current = TRUE, remove = TRUE); gc()
                }
            }
        }
    }
    
    # Turn the options back
    unlink("data/temp", recursive = TRUE)
    terraOptions(memfrac = orig_setting$memfrac, tempdir = orig_setting$tempdir)
}
