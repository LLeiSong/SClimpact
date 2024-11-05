library(rnaturalearth)
library(BIENWorkflow)
# to stop unneeded aux.xml being written with rasters
terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")

#====================================================================
#== 2. Create directories
#====================================================================
mc.cores <- 10
root_dir <- "/home/lsong/SCImpact/data"
myTopDir <- file.path(root_dir, "results")
myEnvDir <- file.path(root_dir, 'variables/Env/')
myOtherEnvDir <- file.path(root_dir, 'variables/OtherEnv/')

if (!file.exists(myTopDir))
    dir.create(myTopDir)
runName <- 'sci' # The unique name for a run
baseDir <- file.path(myTopDir, paste0(runName, '_inputs'))

# Species
mySpeciesDir <- file.path(root_dir, 'occurrences/CSVs_thin/')

# Set background dir based on how to collect background samples

myCustomBackgroundDir <- NULL
myCustomBackgroundTable <- NULL
myCustomDomainRegionRaster <- NULL

mySamplingModelDir <- NULL
mySamplingDataDir <- NULL
myOffsetDirs <- NULL
samplingDataDirNames <- NULL

myAlgorithmNames <- c('PPM')
sortDirNames <- c('PPM')
otherDirs <- list(rdistDir = paste0(baseDir, '/rDists'),
                  domains = paste0(baseDir, '/domains'))
offsetDataDirNames <- c('Expert_Maps')

allBaseDirs <- sdmDirectories(
    baseDir,
    # for inputs
    myAlgorithmNames = myAlgorithmNames,
    warn = FALSE,
    # just set FALSE
    mySpeciesDir = mySpeciesDir,
    # dir for species CSVs
    myEnvDir = myEnvDir,
    # dir for environmental variables
    myOffsetDirs = myOffsetDirs,
    # ??
    mySamplingDataDir = mySamplingDataDir,
    # dir for the samplings
    mySamplingModelDir = mySamplingModelDir,
    # dir for the sampling models
    myCustomBackgroundDir = myCustomBackgroundDir,
    # dir for family pres
    myCustomBackgroundTable = myCustomBackgroundTable,
    myCustomDomainRegionRaster = myCustomDomainRegionRaster,
    # domain for predict
    myOtherEnvDir = myOtherEnvDir,
    # dir for other environmental variables
    offsetDataDirNames = offsetDataDirNames,
    # dir for expert ranges
    samplingDataDirNames = samplingDataDirNames,
    overwrite = FALSE,
    sortDirNames = sortDirNames,
    otherDirs = otherDirs
)

#====================================================================
#== 3. Prep environmental layers
#====================================================================

# Pseudo Mercator
myprj <- st_crs(3857)
pointsProj <- st_crs(3857)

# Preprocess allenv layer
allEnvFile <- paste0(allBaseDirs$envDir, '/AllEnv.tif')
if (!file.exists(allEnvFile)) {
    original.files <- list.files(allBaseDirs$envDir, '.tif', full.names = T)
    env <- terra::rast(original.files)
    terra::crs(env) <- myprj
    terra::writeRaster(env, file = allEnvFile, overwrite = T)
    file.remove(original.files)
}

# Load it into the environment and set layer names
env <- terra::rast(allEnvFile)
layerNames <- try(read.csv(
    list.files(allBaseDirs$envDir, '.csv', full.names = T),
    header = FALSE,
    stringsAsFactors = FALSE
))
if (!class(layerNames) == 'try-error')
    names(env) <- layerNames[, 1]
message('Loaded env raster.')

# Make a mask of nonNA (land)
out.f <- file.path(allBaseDirs$miscDir, 'landMask.tif')
if (!file.exists(out.f)) {
    land <- mean(env)
    land[!is.na(land)] <- 1
    terra::writeRaster(
        land,
        file = out.f,
        overwrite = T,
        datatype = "INT2S",
        gdal = c("COMPRESS=LZW")
    )
}

## shapefile for plotting
world.shp <- ne_countries(scale = 50,
                          type = "countries",
                          returnclass = "sf") %>%
    filter(continent != "Antarctica") %>%
    st_transform(myprj)

## aggregate env layer for subsampling very large species (rather than running the spThin algorithm)
message('Prepped coarse grid for thinning')
coarse.grid.file <- paste0(allBaseDirs$miscDir, '/coarse_grid_for_thinning.tif')
if (!file.exists(coarse.grid.file)) {
    env.template <- terra::rast(allEnvFile)[[1]]
    coarse.grid <- terra::aggregate(env.template, fact = 2)
    terra::writeRaster(coarse.grid, file = coarse.grid.file, overwrite = TRUE)
}

save(allBaseDirs, file = paste0(baseDir, '/allBaseDirs.rdata'))

#====================================================================
#== 4. Sort species
#====================================================================

#--------------------------------------------------------------------
### 4.a sort species by algorithm
#--------------------------------------------------------------------
set.seed(1) # ensures thin algorithm is reproducible
sortDone <- checkSortDone(allBaseDirs, deleteSorted = FALSE)
if (length(sortDone$notSorted) >= 1 &
    !all(is.na(sortDone$notSorted))) {
    sortSpeciesBulk(
        allBaseDirs,
        pointsProj = pointsProj,
        myprj = myprj,
        doThin = TRUE,
        thinCutoff = 5,
        # doThin for more than 5 samples
        maxCells = 20000,
        # max number of samples to retain
        verbose = TRUE,
        overwrite = FALSE,
        doClean = TRUE,
        # retain only one samples for each grid
        mc.cores = mc.cores,
        # parallel
        longLatNames = c('x', 'y'),
        colsToKeep = 'id'
    )
}

#====================================================================
#== 5. Model settings
#====================================================================
cbs <- commonBIENSettings(myprj, world.shp, runName = runName)
toDo <- cbs$toDo
toSave <- cbs$toSave
toOverwrite <- cbs$toOverwrite

# these are changes to the defaults
toDo$misc$algorithm <- 'PPM'
toDo$sampling$samplingSettings <- c('noBias')
toDo$fitting$expertSettings <- list(
    prob = 1e-6,
    rate = 0,
    skew = 1e-6,
    shift = 0)
toDo$fitting$predictorSettings <- NULL
toDo$fitting$priorSettings <- NULL
toDo$fitting$algorithmSettings <- c('maxnet')
toDo$plotting$nPlotCol <- 4
toDo$plotting$nPlotRow <- 3
toDo$plotting$whichFuturesToPlot <- c('gf')

toDo$plotting$plotJustASampleOfMaps <- .02
toDo$plotting$responseCurves <- FALSE
toDo$transfer$otherEnvProjections <- TRUE
toDo$plotting$plotOtherEnvPred <- TRUE

toDo$background$makeBiasedBackground <- FALSE
toDo$background$customBG <- FALSE
toDo$fitting$evalFork <- FALSE
toDo$domain$standardizePredictors <- FALSE

# write out run setting used for metadata
runSettings <- list(toDo, toSave, toOverwrite)
save(runSettings, file = paste0(baseDir, '/runSettings.rdata'))
