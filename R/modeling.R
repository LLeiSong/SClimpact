#====================================================================
#== 1. Install packages
#====================================================================
rm(list=ls())
library(here)
library(BIENWorkflow)
# to stop unneeded aux.xml being written with rasters
terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")

#====================================================================
#== 2. Create directories
#====================================================================
myTopDir <- file.path("results/10kmCT")
mc.cores <- 10
myEnvDir <- 'results/experiments/10km_var_br_lc/AllEnv/'
myOtherEnvDir <- 'results/experiments/10km_var_br_lc/OtherEnv/'

if(!file.exists(myTopDir)) dir.create(myTopDir)
runName <- 'ele' # The unique name for a run
baseDir <- file.path(myTopDir, paste0(runName,'_inputs'))

# Species
mySpeciesDir <- here::here("results/experiments/CSVs")

# Set background dir based on how to collect background samples

myCustomBackgroundDir <- NULL
myCustomBackgroundTable <- NULL
myCustomDomainRegionRaster <- NULL

mySamplingModelDir <- NULL
mySamplingDataDir <- NULL
myOffsetDirs <- NULL
samplingDataDirNames <- NULL

myAlgorithmNames <- c('PPM','Points','RangeBag','PPM_RangeBag_As_Points','RangeBag_1530')
sortDirNames <- c('Points','RangeBag','PPM','PPM15_30','RangebagFromPPM', 'PointsFromRangeBag')
# CM 6/19/23: 'RangebagFromPPM', 'PointsFromRangeBag' are directories to move species with small sample sizes that are found during cleaning.
otherDirs <- list(rdistDir = paste0(baseDir,'/rDists'),
                  domains = paste0(baseDir,'/domains'))
offsetDataDirNames <- c('Expert_Maps')

allBaseDirs <- sdmDirectories(
    baseDir, # for inputs
    myAlgorithmNames = myAlgorithmNames,
    warn = FALSE, # just set FALSE
    mySpeciesDir = mySpeciesDir, # dir for species CSVs
    myEnvDir = myEnvDir, # dir for environmental variables
    myOffsetDirs = myOffsetDirs, # ??
    mySamplingDataDir = mySamplingDataDir, # dir for the samplings
    mySamplingModelDir = mySamplingModelDir, # dir for the sampling models
    myCustomBackgroundDir = myCustomBackgroundDir, # dir for family pres
    myCustomBackgroundTable = myCustomBackgroundTable,
    myCustomDomainRegionRaster = myCustomDomainRegionRaster, # domain for predict					 
    myOtherEnvDir = myOtherEnvDir, # dir for other environmental variables
    offsetDataDirNames = offsetDataDirNames, # dir for expert ranges
    samplingDataDirNames = samplingDataDirNames,
    overwrite = FALSE,
    sortDirNames = sortDirNames,
    otherDirs = otherDirs)

#====================================================================
#== 3. Prep environmental layers
#====================================================================

# Mollweide
# +proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs
myprj <- st_crs("EPSG:3857")
# Geographic coordinate:'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
pointsProj <- st_crs(4326)

# Preprocess allenv layer
allEnvFile <- paste0(allBaseDirs$envDir,'/AllEnv.tif')
if(!file.exists(allEnvFile)){
    original.files <- list.files(allBaseDirs$envDir,'.tif', full.names=T)
    env <- terra::rast(original.files)
    terra::crs(env) <- myprj
    terra::writeRaster(env,file = allEnvFile, overwrite = T)
    file.remove(original.files)
}

# Load it into the environment and set layer names
env <- terra::rast(allEnvFile)
layerNames <- try(
    read.csv(list.files(allBaseDirs$envDir,'.csv', full.names=T), 
             header=FALSE, stringsAsFactors=FALSE))
if(!class(layerNames) == 'try-error') names(env) <- layerNames[, 1]
message('Loaded env raster.')

# Make a mask of nonNA (land)
out.f <- file.path(allBaseDirs$miscDir, 'landMask.tif')
if(!file.exists(out.f)){
    land <- mean(env)
    land[!is.na(land)] <- 1 
    terra::writeRaster(
        land,file = out.f, overwrite = T, 
        datatype = "INT2S", gdal = c("COMPRESS=LZW"))
}

## shapefile for plotting 
world.shp <- ne_countries(scale = 50, type = "countries", returnclass = "sf") %>% 
    filter(continent != "Antarctica") %>% 
    st_transform(myprj)

## aggregate env layer for subsampling very large species (rather than running the spThin algorithm)
message('Prepped coarse grid for thinning')
coarse.grid.file <- paste0(allBaseDirs$miscDir, '/coarse_grid_for_thinning.tif')
if(!file.exists(coarse.grid.file)){
    env.template <- terra::rast(allEnvFile)[[1]]
    coarse.grid <- terra::aggregate(env.template, fact = 2)
    terra::writeRaster(coarse.grid, file = coarse.grid.file, overwrite = TRUE)
}

save(allBaseDirs, file = paste0(baseDir,'/allBaseDirs.rdata'))

#====================================================================
#== 4. Sort species
#====================================================================

#--------------------------------------------------------------------
### 4.a sort species by algorithm
#--------------------------------------------------------------------
# Note that to get the full tracking information on points used, your csvs must contain a column called ID, so that omitted records can be tracked. the name of this column is hard corded into the workflow.
set.seed(1) # ensures thin algorithm is reproducible
sortDone <- checkSortDone(allBaseDirs, deleteSorted = FALSE)
if(length(sortDone$notSorted) >= 1 & !all(is.na(sortDone$notSorted))) {
    sortSpeciesBulk(
        allBaseDirs, 
        pointsProj = pointsProj, 
        myprj = myprj,
        doThin = TRUE,
        thinCutoff = 5, # doThin for more than 5 samples
        maxCells = 20000, # max number of samples to retain
        verbose = TRUE,
        overwrite = FALSE,
        doClean = TRUE, # retain only one samples for each grid
        mc.cores = mc.cores, # parallel
        longLatNames = c('lon','lat'),
        colsToKeep = 'gbif_id')
}

### 4.b specify which species to run
#--------------------------------------------------------------------
speciesList <- c(
    list.files(file.path(allBaseDirs$speciesDir,'PPM15_30'), full.names = TRUE),
    list.files(file.path(allBaseDirs$speciesDir,'PPM'), full.names = TRUE))

# # Check some facts
# done <- checkSDMDone(allBaseDirs, speciesList, algorithm='PPM', checkFigs = T) 
# print('PPM status')
# str(done)
# speciesList=done$notRun
# speciesList=done$problems # for debugging

#====================================================================
#== 5. Model settings
#====================================================================
cbs <- commonBIENSettings(myprj, world.shp, runName = runName)
toDo <- cbs$toDo
toSave <- cbs$toSave
toOverwrite <- cbs$toOverwrite

# these are changes to the defaults
toDo$misc$algorithm <- 'PPM' # must match name specified in myAlgorithmNames
toDo$sampling$samplingSettings <- c('noBias')#,'targetBG') # note that the target group not tested
toDo$fitting$expertSettings <- list(prob=1e-6, rate=0, skew=1e-6, shift=0)
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

toDo$background$makeBiasedBackground=FALSE
toDo$background$customBG=FALSE
toDo$fitting$evalFork = FALSE
toDo$domain$standardizePredictors <- FALSE

# write out run setting used for metadata
runSettings <- list(toDo,toSave,toOverwrite)
save(runSettings, file = paste0(baseDir, '/runSettings.rdata'))

#--------------------------------------------------------------------
### 6. Run SDM Models
#--------------------------------------------------------------------
tt1 <- proc.time()
if(length(speciesList) > 0){
    lapply(1:length(speciesList), function(j){
        out <- sdmBIEN(
            speciesCSV = speciesList[j],
            allBaseDirs,
            toDo = toDo,
            toSave = toSave,
            toOverwrite = toOverwrite) 
        if(class(out) == 'try-error') print(out)
        gc()
    })
    (tt2 = proc.time() - tt1)
}

# # clean up problematic 
# done = checkSDMDone(allBaseDirs, speciesList, algorithm = 'PPM') 
# str(done)

#--------------------------------------------------------------------
### 7. Run Range bagging Models
#--------------------------------------------------------------------

speciesListRB <- c(
    list.files(file.path(allBaseDirs$speciesDir,'RangeBag'), full.names = T),
    list.files(file.path(allBaseDirs$speciesDir,'RangebagFromPPM'), full.names = T))
doneRB <- checkSDMDone(allBaseDirs, speciesList = speciesListRB, 
                       algorithm = 'RangeBag', checkFigs = FALSE) 
speciesListRB <- doneRB$notRun

#Settings
future <- TRUE
csrb <- commonRangeBaggingSettings(myprj, world.shp, runName = runName)
toDo <- csrb$toDo
toSave <- csrb$toSave
toOverwrite <- csrb$toOverwrite
toDo$misc$bienMetadata <- TRUE
if(future){
    toDo$transfer$otherEnvProjections <- TRUE
    toSave$shapeFileTransfer <- FALSE
    toDo$transfer$whichFutures <- 'all'
    toDo$plotting$whichFuturesToPlot <- c('gf','ip')
    toDo$plotting$plotJustASampleOfMaps <- .02
    toDo$plotting$plotOtherEnvPred <- TRUE
    toDo$plotting$nPlotCol <- 4
    toDo$plotting$nPlotRow <- 3
} else {
    toDo$plotting$nPlotCol <- 2
    toDo$plotting$plotJustASampleOfMaps <- 1
}

# write out run setting used for metadata
runSettingsRB <- list(toDo, toSave, toOverwrite)
save(runSettingsRB, file = paste0(baseDir, '/runSettingsRangeBag.rdata'))
# Run Models
if(length(speciesListRB) > 0){
    mclapply(1:length(speciesListRB), function(j){
        out <- sdmRangeBag(
            speciesCSV = speciesListRB[j],
            allBaseDirs,
            toDo = toDo,
            toSave = toSave,
            toOverwrite = toOverwrite) 
        if(class(out) == 'try-error') print(out)													 
        gc()
    }, mc.cores = mc.cores)
}

#--------------------------------------------------------------------
### 8. Run Models for species with too few points for a model
#------------------------------------------------------------
speciesListPt <- c(
    list.files(file.path(allBaseDirs$speciesDir, 'Points'), full.names = T),
    list.files(file.path(allBaseDirs$speciesDir, 'PointsFromRangebag'), full.names = T))

# for picking up after where you left off after a run has started
donePt <- checkSDMDone(
    allBaseDirs, speciesList = speciesListPt, 
    algorithm = 'Points',checkFigs = FALSE)
speciesListPt <- donePt$notRun

# Settings
csp <- commonPointsSettings(myprj, world.shp, runName=runName)
toDo <- csp$toDo
toSave <- csp$toSave
toOverwrite <- csp$toOverwrite

toDo$plotting$plotMaps <- FALSE  
toDo$plotting$plotJustASampleOfMaps <- 0

# write out run setting used for metadata
runSettingsPoints <- list(toDo, toSave, toOverwrite)
save(runSettingsPoints, file = paste0(baseDir,'/runSettingsPoints.rdata'))

# Run Models
if(length(speciesListPt > 0)){
    mclapply(1:length(speciesListPt), function(j){
        out <- pointsMap(
            speciesCSV = speciesListPt[j],
            allBaseDirs,
            toDo = toDo,
            toSave = toSave,
            toOverwrite = toOverwrite) 
        if(class(out) == 'try-error') print(out)													 
        gc()
    }, mc.cores = mc.cores)
}

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# 9. Make maps of just the raw presence records for all species
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# this corresponds to an extra step running the 'Points' algorithm for species that were modeled with Range bagging or PPMs
speciesListPPMRBAsPt <- c(
    list.files(file.path(allBaseDirs$speciesDir, 'PPM'), full.names = T),
    list.files(file.path(allBaseDirs$speciesDir, 'PPM15_30'), full.names = T),
    list.files(file.path(allBaseDirs$speciesDir, 'RangeBag'), full.names = T),
    list.files(file.path(allBaseDirs$speciesDir, 'RangebagFromPPM'), 
               full.names = T))

# check whether any already done
donePt <- checkSDMDone(
    allBaseDirs, speciesList = speciesListPPMRBAsPt, 
    algorithm = 'PPM_RangeBag_As_Points', checkFigs = FALSE)
speciesListPPMRBAsPt <- donePt$notRun

# Settings
csp <- commonPointsSettings(myprj,world.shp,runName=runName)
toDo <- csp$toDo
toSave <- csp$toSave
toOverwrite <- csp$toOverwrite
toDo$misc$algorithm <- 'PPM_RangeBag_As_Points'
toDo$plotting$plotMaps <- FALSE  
toDo$plotting$plotJustASampleOfMaps <- 0
# write out run setting used for metadata
runSettingsPoints <- list(toDo,toSave,toOverwrite)
save(runSettingsPoints, file = paste0(baseDir,'/runSettingsPPMRBAsPoints.rdata'))

# Run Models
if(length(speciesListPPMRBAsPt)>0){
    mclapply(1:length(speciesListPPMRBAsPt), function(j) {
        out <- pointsMap(
            speciesCSV = speciesListPPMRBAsPt[j],
            allBaseDirs,
            toDo = toDo,
            toSave = toSave,
            toOverwrite = toOverwrite)
        if (class(out) == 'try-error') print(out)
        gc()
    }, mc.cores = mc.cores)
}

#--------------------------------------------------------------------
### 10. Do Range Bagging for PPM15-30 in case we want to avoid PPMs for them
#--------------------------------------------------------------------
speciesListRB1530 <- list.files(
    file.path(allBaseDirs$speciesDir,'PPM15_30'), full.names = T)
doneRB1530 <- checkSDMDone(
    allBaseDirs, speciesList = speciesListRB1530, 
    algorithm='RangeBag_1530', checkFigs = FALSE) 
speciesListRB1530 <- doneRB1530$notRun

#Settings
future <- TRUE
csrb <- commonRangeBaggingSettings(myprj,world.shp,runName=runName)
toDo <- csrb$toDo
toSave <- csrb$toSave
toOverwrite <- csrb$toOverwrite
toDo$misc$bienMetadata <- TRUE
toDo$misc$algorithm <- 'RangeBag_1530'
toDo$transfer$whichFutures <- 'all'	

if(future) {
    toDo$transfer$otherEnvProjections <- TRUE
    toSave$shapeFileTransfer <- TRUE
    toDo$transfer$whichFutures <- 'all'
    toDo$plotting$whichFuturesToPlot <- c('gf', 'ip')
    toDo$plotting$plotJustASampleOfMaps <- .02
    toDo$plotting$plotOtherEnvPred <- TRUE
    toDo$plotting$nPlotCol <- 4
    toDo$plotting$nPlotRow <- 3
} else {
    toDo$plotting$nPlotCol <- 2
    toDo$plotting$plotJustASampleOfMaps <- 1
}

# write out run setting used for metadata
runSettingsRB <- list(toDo, toSave, toOverwrite)
save(runSettingsRB, file = paste0(baseDir,'/runSettingsRangeBag1530.rdata'))

# Run Models
if(length(speciesListRB1530) > 0){
    mclapply(1:length(speciesListRB1530), function(j){
        out <- sdmRangeBag(
            speciesCSV = speciesListRB1530[j],
            allBaseDirs,
            toDo = toDo,
            toSave = toSave,
            toOverwrite = toOverwrite) 
        if(class(out) == 'try-error') print(out)													 
        gc()
    }, mc.cores = mc.cores)
}
