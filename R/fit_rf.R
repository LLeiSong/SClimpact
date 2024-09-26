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
fit_rf <- function(species, occ, domain, vars, vars_future, n_core, dst_dir){
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
    
    # Update the cores
    # ram_need <- mem_info(vars, n = 1)
    # n_core <- min(floor(ram_need[2] * ram_need[3] / ram_need[1]), n_core)
    
    # Dynamically make background samples of 1000
    set.seed(123) # make sure variability for each parameter pair
    
    occ <- occ %>% mutate(observation = 1) %>% dplyr::select(observation, folds)
    
    if (nrow(occ) > 5000){
        set.seed(123)
        occ <- occ %>% sample_n(size = 5000)
    }
    
    obs_data <- BIOMOD_FormatingData(
        resp.var = occ$observation,
        expl.var = vars,
        dir.name = dst_dir,
        resp.xy = st_coordinates(occ),
        resp.name = species,
        PA.nb.rep = max(occ$folds),
        PA.strategy = "random",
        PA.nb.absences = nrow(occ),
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

    # Run the model
    models <- BIOMOD_Modeling(
        bm.format = obs_data,
        models = c("RF"),
        modeling.id = "RF",
        CV.strategy = "user.defined",
        CV.user.table = CVtable,
        var.import = 3, 
        OPT.strategy = 'bigboss',
        metric.eval = c("ACCURACY", 'TSS', 'ROC'),
        CV.do.full.models = TRUE,
        seed.val = 123,
        nb.cpu = n_core)
    
    # Variable importance
    fn <- file.path(dst_dir, sprintf("evaluation_%s.csv", gsub(" ", "_", species)))
    if (!file.exists(fn)){
        eval <- get_evaluations(models)
        write.csv(eval, fn, row.names = FALSE)
    }
    
    # Only use full name for results
    mod_nm <- sprintf("%s_allData_allRun_RF", gsub(" ", ".", species))
    
    # Predict
    fn <- file.path(dst_dir, sprintf("current_suit_%s.tif", gsub(" ", "_", species)))
    if (!file.exists(fn)){
        model_proj <- BIOMOD_Projection(
            bm.mod = models,
            proj.name = "Current",
            new.env = vars,
            models.chosen = mod_nm,
            metric.binary = "all",
            metric.filter = "all",
            nb.cpu = n_core)
        pred <- get_predictions(model_proj)
        writeRaster(pred, fn)
    }
    
    fn <- file.path(dst_dir, sprintf("future_suit_%s.tif", gsub(" ", "_", species)))
    if (!file.exists(fn)){
        model_proj_future <- BIOMOD_Projection(
            bm.mod = models,
            proj.name = "Future",
            new.env = vars_future,
            models.chosen = mod_nm,
            metric.binary = "all",
            metric.filter = "all",
            nb.cpu = n_core)
        pred_future <- get_predictions(model_proj_future)
        writeRaster(pred_future, fn)
    }

    # Variable importance
    fn <- file.path(dst_dir, sprintf("variable_importance_%s.csv", gsub(" ", "_", species)))
    if (!file.exists(fn)){
        var_imp <- get_variables_importance(models)
        write.csv(var_imp, fn, row.names = FALSE)
    }
    
    # Get spatial SHAP
    pfun <- function(X.model, newdata) {
        predict(X.model, newdata, type = "prob")[, "1"] # RF
    }
    
    # Calculate SHAP values
    BIOMOD_LoadModels(models, mod_nm)
    
    fn <- file.path(dst_dir, sprintf("future_shap_%s.tif", gsub(" ", "_", species)))
    if (!file.exists(fn)){
        new_data <- values(vars_future) %>% na.omit() %>% data.frame() %>% 
            mutate(lc = factor(lc, levels = cls$id, labels = cls$cover))
        
        shap_explain <- fastshap::explain(
            get(mod_nm)@model,
            X = obs_data@data.env.var,
            newdata = new_data,
            nsim = 200,
            adjust = TRUE,
            pred_wrapper = pfun)
        
        fut_shap <- vars_future
        vals <- values(fut_shap)
        vals[complete.cases(vals), ] <- shap_explain
        values(fut_shap) <- vals
        writeRaster(fut_shap, fn)
    }
    
    fn <- file.path(dst_dir, sprintf("current_shap_%s.tif", gsub(" ", "_", species)))
    if (!file.exists(fn)){
        new_data <- values(vars) %>% na.omit() %>% data.frame() %>% 
            mutate(lc = factor(lc, levels = cls$id, labels = cls$cover))
        
        shap_explain <- fastshap::explain(
            get(mod_nm)@model,
            X = obs_data@data.env.var,
            newdata = new_data,
            nsim = 200,
            adjust = TRUE,
            pred_wrapper = pfun)
        
        cur_shap <- vars
        vals <- values(cur_shap)
        vals[complete.cases(vals), ] <- shap_explain
        values(cur_shap) <- vals
        
        
        writeRaster(cur_shap, fn)
    }
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
                                help = "The row of species to process."),
                    make_option(c("-c", "--cores"), 
                                action = "store", type = 'integer',
                                help = "The number of cores to use."))
opt <- parse_args(OptionParser(option_list = option_list))
i <- opt$index
n_core <- opt$cores
terraOptions(memfrac = 0.9)

# Collect occurrences
occ_list <- list.files(occ_dir, full.names = TRUE) %>% 
    data.frame(fname = .) %>% 
    mutate(species = gsub(".geojson", "", basename(fname))) %>% 
    mutate(species = gsub("_", " ", species))

# # Remove species have been done
# done <- list.files(file.path(dst_dir), pattern = "current_shap")
# done <- gsub("current_shap_", "", done)
# done <- gsub(".tif", "", done)
# done <- gsub("_", " ", done)
# occ_list <- occ_list %>% filter(!species %in% done); rm(done)

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
vars <- aggregate(vars, 6)
vars$lc <- aggregate(rast(file.path("data/env", "LULC_2015.tif")), 6, fun = "modal")
levels(vars$lc) <- cls

vars_future <- rast(
    file.path("data/env", sprintf("chelsa_%s-%s.tif", 
                                  fut_year - 30 + 1, fut_year)))
vars_future <- aggregate(vars_future, 6)
vars_future$lc <- aggregate(rast(file.path("data/env", sprintf("LULC_%s.tif", fut_year))), 6,
                            fun = "modal")

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

# Fit RF model and calculate SHAP
fit_rf(species, occ, domain, vars, vars_future, n_core, dst_dir)
