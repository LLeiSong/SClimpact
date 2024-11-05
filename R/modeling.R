# Load the parameters
load(paste0(baseDir, '/allBaseDirs.rdata'))
load(paste0(baseDir, '/runSettings.rdata'))
toDo <- runSettings$toDo
toSave <- runSettings$toSave
toOverwrite <- runSettings$toOverwrite

### 4.b specify which species to run
#--------------------------------------------------------------------
speciesList <- list.files(
    file.path(allBaseDirs$speciesDir,'PPM'), full.names = TRUE)

#--------------------------------------------------------------------
### 6. Run SDM Models
#--------------------------------------------------------------------
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
}
