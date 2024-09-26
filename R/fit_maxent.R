# Configure HPC
# First, install conda on HPC
# wget https://repo.anaconda.com/archive/Anaconda3-2023.07-2-Linux-x86_64.sh
# sh Anaconda3.*.sh
# Then: 
# 1. source ~/anaconda3/bin/activate
# 2. conda create --name sdm
# 3. conda activate sdm
## Install R and relevant packages:
## conda install conda-forge::r-essentials
## conda install conda-forge::r-terra
## conda install conda-forge::r-sf
## conda install r-fastshap
## conda install conda-forge::r-optparse
## conda install conda-forge::r-biomod2
## and other possible ones

# Load libraries
library(biomod2)
library(dplyr)
library(terra)
library(sf)
library(stringr)
library(fastshap)
library(optparse)

# Define the function
## Fit SDM with different parameters
fit_maxent <- function(species, occ, domain, vars, vars_future, dst_dir){
    message(species)
    
    # Mask vars by domain
    cls <- levels(vars$lc)[[1]]
    vars <- mask(crop(vars, domain), domain)
    
    # Make a non-NA mask for background sampling
    msk <- lapply(1:nlyr(vars), function(n) !is.na(vars[[n]]))
    msk <- sum(do.call(c, msk)) == nlyr(vars)
    msk[msk == 0] <- NA
    
    # Update vars
    vars <- vars * msk
    levels(vars$lc) <- cls
    
    # Update vars_future
    vars_future <- mask(crop(vars_future, domain), domain)
    
    # Make a non-NA mask for background sampling
    msk <- lapply(1:nlyr(vars_future), function(n) !is.na(vars_future[[n]]))
    msk <- sum(do.call(c, msk)) == nlyr(vars_future)
    msk[msk == 0] <- NA
    
    # Update vars
    vars_future <- vars_future * msk
    levels(vars_future$lc) <- cls
    
    # Dynamically make background samples of 1000
    set.seed(123) # make sure variability for each parameter pair
    
    occ <- occ %>% mutate(observation = 1) %>% dplyr::select(observation, folds)
    
    obs_data <- BIOMOD_FormatingData(
        resp.var = occ$observation,
        expl.var = vars,
        dir.name = dst_dir,
        resp.xy = st_coordinates(occ),
        resp.name = species,
        PA.nb.rep = max(occ$folds),
        PA.strategy = "random",
        PA.nb.absences = min(1000, nrow(occ)),
        na.rm = TRUE)
    
    CVtable <- do.call(cbind, lapply(sort(unique(occ$folds)), function(i){
        set.seed(i)
        real <- occ$folds != i
        bg <- obs_data@PA.table[, i]
        bg[bg == FALSE] <- NA
        bg <- bg[c(nrow(occ) + 1):length(bg)]
        bg[!is.na(bg)][sample(length(bg[!is.na(bg)]), sum(!real))] <- FALSE
        c(real, bg)
    })) %>% data.frame()
    names(CVtable) <- sprintf("_PA%s_RUN%s", 
                              sort(unique(occ$folds)), sort(unique(occ$folds)))
    CVtable$`_allData_allRun` <- obs_data@PA.table[, 1]
    CVtable$`_allData_allRun`[CVtable$`_allData_allRun` == FALSE] <- NA
    
    # Set parameters
    ## BIOMOD2 4.2-5-2
    user.maxent <- rep(list(list(memory_allocated = 1024 * 10, 
                                 max_heap_size = '10G')), max(occ$folds) + 1)
    names(user.maxent) <- c('_allData_allRun', 
                            sprintf("_PA%s_allRun", sort(unique(occ$folds))))
    
    ## BIOMOD2 4.2-6
    # user.maxent <- list('for_all_datasets' = list(memory_allocated = 1024 * 10, 
    #                                               max_heap_size = '10G'))
    
    user.val <- list(MAXENT.binary.MAXENT.MAXENT = user.maxent)
    opt.u <- bm_ModelingOptions(data.type = 'binary',
                                models = c('MAXENT'),
                                strategy = 'user.defined',
                                user.val = user.val,
                                bm.format = obs_data)
    for (nm in names(CVtable)){
        opt.u@options$MAXENT.binary.MAXENT.MAXENT@args.values[[nm]] <- 
            opt.u@options$MAXENT.binary.MAXENT.MAXENT@args.values$for_all_datasets
    }
    
    # Run the model
    models <- BIOMOD_Modeling(
        bm.format = obs_data,
        models = c("MAXENT"),
        modeling.id = "maxent",
        CV.strategy = "user.defined",
        CV.user.table = CVtable,
        var.import = 3, 
        metric.eval = c("ACCURACY", 'TSS', 'ROC'),
        OPT.user = opt.u, 
        CV.do.full.models = TRUE,
        seed.val = 123,
        nb.cpu = 10)
    
    # Only use full name for results
    mod_nm <- sprintf("%s_allData_allRun_MAXENT", gsub(" ", ".", species))
    model_proj <- BIOMOD_Projection(
        bm.mod = models,
        proj.name = "Current",
        new.env = vars,
        models.chosen = mod_nm,
        metric.binary = "all",
        metric.filter = "all",
        nb.cpu = 10)
    
    model_proj_future <- BIOMOD_Projection(
        bm.mod = models,
        proj.name = "Future",
        new.env = vars_future,
        models.chosen = mod_nm,
        metric.binary = "all",
        metric.filter = "all",
        nb.cpu = 10)
    
    pred <- get_predictions(model_proj)
    pred_future <- get_predictions(model_proj_future)
    
    fn <- file.path(dst_dir, sprintf("current_suit_%s.tif", gsub(" ", "_", species)))
    writeRaster(pred, fn)
    
    fn <- file.path(dst_dir, sprintf("future_suit_%s.tif", gsub(" ", "_", species)))
    writeRaster(pred_future, fn)
    
    # Variable importance
    var_imp <- get_variables_importance(models)
    fn <- file.path(dst_dir, sprintf("variable_importance_%s.csv", gsub(" ", "_", species)))
    write.csv(var_imp, fn, row.names = FALSE)
    
    # Get spatial SHAP
    pfun <- function(X.model, newdata) {
        # predict(X.model, newdata, type = "prob")[, "1"] # RF
        nm <- X.model@models.computed
        nm <- nm[str_detect(nm, "allData_allRun_MAXENT")]
        pred <- BIOMOD_Projection(
            bm.mod = X.model,
            proj.name = "Future",
            new.env = newdata,
            models.chosen = nm,
            nb.cpu = 10)
        
        # results
        get_predictions(pred)$pred / 1000
    }
    
    # Calculate SHAP values
    new_data <- values(vars_future) %>% na.omit() %>% data.frame() %>% 
        mutate(lc = factor(lc, levels = cls$id, labels = cls$cover))
    shap_explain <- explain(
        models,
        X = obs_data@data.env.var,
        newdata = new_data,
        nsim = 50,
        adjust = TRUE,
        pred_wrapper = pfun)
    
    fut_shap <- vars_future
    vals <- values(fut_shap)
    vals[complete.cases(vals), ] <- shap_explain
    values(fut_shap) <- vals
    
    fn <- file.path(dst_dir, sprintf("future_shap_%s.tif", gsub(" ", "_", species)))
    writeRaster(fut_shap, fn)
    
    new_data <- values(vars) %>% na.omit() %>% data.frame() %>% 
        mutate(lc = factor(lc, levels = cls$id, labels = cls$cover))
    shap_explain <- explain(
        models,
        X = obs_data@data.env.var,
        newdata = new_data,
        nsim = 50,
        adjust = TRUE,
        pred_wrapper = pfun)
    
    cur_shap <- vars
    vals <- values(cur_shap)
    vals[complete.cases(vals), ] <- shap_explain
    values(cur_shap) <- vals
    
    fn <- file.path(dst_dir, sprintf("current_shap_%s.tif", gsub(" ", "_", species)))
    writeRaster(cur_shap, fn)
}

# Directories
data_dir <- "data"
occ_dir <- file.path(data_dir, "gbif_clean")
dst_dir <- file.path(data_dir, "sdm")
if (!dir.exists(dst_dir)) dir.create(dst_dir)

# Get input
# Parse inline arguments
option_list <- list(make_option(c("-i", "--index"), 
                                action = "store", type = 'integer',
                                help = "The row of species to process."))
opt <- parse_args(OptionParser(option_list = option_list))
i <- opt$index

# Collect occurrences
occ_list <- list.files(occ_dir, full.names = TRUE) %>% 
    data.frame(fname = .) %>% 
    mutate(species = gsub(".geojson", "", basename(fname))) %>% 
    mutate(species = gsub("_", " ", species))

# Load environmental variables
fut_year <- 2070
cls <- data.frame(
    id = 1:7,
    cover = c("Water", "Forest", 
              "Grassland",
              "Barren",
              "Cropland",
              "Urban",
              "Permanent snow and ice"))

vars <- rast(file.path("data/env", "chelsa_1981-2010.tif"))
vars$lc <- rast(file.path("data/env", "LULC_2015.tif"))
levels(vars$lc) <- cls
vars_future <- rast(
    file.path("data/env", sprintf("chelsa_%s-%s.tif", 
                                  fut_year - 30 + 1, fut_year)))
vars_future$lc <- rast(file.path("data/env", sprintf("LULC_%s.tif", fut_year)))
levels(vars_future$lc) <- cls

# Load other data
load("data/IUCN/mammals_endangered.rda")
ecoregions <- st_read(
    file.path(data_dir, "terr-ecoregions-TNC", "tnc_terr_ecoregions.shp"))

# Other settings
sf_use_s2(FALSE)

# Start process
species <- occ_list[i, "species"]
occ <- st_read(occ_list[i, "fname"])

# Get domain
range_map <- range_maps %>% filter(sci_name == species)
ecoregion <- ecoregions %>% 
    slice(unique(unlist(st_intersects(range_map, .))))
domain <- ecoregions %>% 
    slice(unique(unlist(st_intersects(ecoregion %>% st_buffer(0.0083), .))))
domain <- domain %>% st_union() %>% st_simplify() %>% 
    st_as_sf()
rm(range_map, ecoregion)

# Fit maxent model and calculate SHAP
fit_maxent(species, occ, domain, vars, vars_future, dst_dir)
    