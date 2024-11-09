#' @title rf_dws
#' @description Fit a down sampled RF SDM and do prediction.
#' @param sp (`character`) The species name with format Name_name.
#' @param occ_dir (`character`) The directory to save thinned occurrences.
#' @param var_dir (`character`) The directory for environmental variables.
#' @param range_dir (`character`) The directory for IUCN ranges.
#' @param dst_dir (`character`) The directory to save the result.
#' @param seed (`integer`) The seed for randomization. 
#' @return All results are saved to `dst_dir` under the folder with same name
#' as sp.
#' 
#' @import blockCV
#' @import precrec
#' @import randomForest
#' @import terra
#' @import sf
#' @import dplyr
#' @export
#' @examples
#' \donttest{
#' sp <- "Akodon_albiventer"
#' occ_dir <- "data/occurrences"
#' var_dir <- "data/variables"
#' range_dir <- "data/IUCN/Expert_Maps"
#' dst_dir <- "results/sdm"
#' rf_dws(sp, occ_dir, var_dir, range_dir, dst_dir)
#' }
#' 

rf_dws <- function(sp, 
                   occ_dir = "data/occurrences",
                   var_dir = "data/variables",
                   range_dir = "data/IUCN/Expert_Maps",
                   dst_dir = "results/sdm",
                   seed = 123){
    ########### Create folder for results ###########
    if (!dir.exists(dst_dir)) dir.create(dst_dir, recursive = TRUE)
    sp_dir <- file.path(dst_dir, sp)
    if (dir.exists(sp_dir)) unlink(sp_dir, recursive = TRUE)
    dir.create(sp_dir)
    
    ########### Load data ###########
    ## Points
    occ <- read.csv(file.path(occ_dir, "CSVs_thin", sprintf("%s.csv", sp)))
    bg <- read.csv(file.path(occ_dir, "bg", sprintf("%s.csv", sp)))
    
    ## Environmental variables
    vars <- rast(file.path(var_dir, "Env/AllEnv.tif"))
    var_list <- read.csv(
        file.path(var_dir, "variable_list", sprintf("%s.csv", sp)))
    if (length(na.omit(var_list$var_uncorrelated)) == 1){
        var_list <- na.omit(var_list$var)
    } else var_list <- na.omit(var_list$var_uncorrelated)
    vars <- subset(vars, var_list)
    
    # Load range polygon
    range <- st_read(
        file.path(range_dir, sprintf("%s.geojson", sp)), 
        quiet = TRUE) %>% st_transform(crs(vars))
    
    ########### Split for cross validation ###########
    ## Use NbClust to automatically decide the optimal number of clusters
    ## Min is 5 and max is 15 to save time.
    set.seed(seed)
    sac <- cv_spatial_autocor(
        x = st_as_sf(occ, coords = c("x", "y"), 
                     crs = 3857) %>% mutate(occ = 1),
        column = "occ", plot = FALSE, progress = FALSE)
    
    scv <- tryCatch({
        cv_spatial(
            x = st_as_sf(occ, coords = c("x", "y"), 
                         crs = 3857) %>% mutate(occ = 1),
            column = "occ",
            # in case for super sparse distributed species
            k = min(5, nrow(sac$plots$data)),
            size = sac$range,
            selection = "random",
            biomod2 = FALSE,
            iteration = 200, seed = seed,
            progress = FALSE, report = FALSE, plot = FALSE) 
    }, error = function(e){
        cv_spatial(
            x = st_as_sf(occ, coords = c("x", "y"), 
                         crs = 3857) %>% mutate(occ = 1),
            column = "occ",
            # in case for super sparse distributed species
            k = min(5, nrow(sac$plots$data)),
            size = sac$range,
            selection = "random",
            biomod2 = FALSE, hexagon = FALSE, # For one super spread species
            iteration = 200, seed = seed,
            progress = FALSE, report = FALSE, plot = FALSE) 
    })
    
    # Merge the folds
    occ <- occ %>% mutate(fold = scv$folds_ids)
    
    # Get pseudo absence, assuming the points outside of range are most likely
    # be true absence for a solid evaluation.
    ## Notice: species living in relatively islands may failed to extract
    ## any bg_eval using this method. If this is the case, randomly select
    ## the required number of bg_eval.
    bg_eval <- tryCatch({
        set.seed(seed)
        bg %>% st_as_sf(coords = c("x", "y"), crs = crs(vars)) %>% 
            vect() %>% mask(., vect(range %>% st_buffer(50000)), 
                            inverse = TRUE) %>% st_as_sf() %>% 
            sample_n(size = nrow(occ), replace = TRUE) %>% 
            st_coordinates() %>% as.data.frame() %>% 
            mutate(sp = unique(occ$sp), id = 1:nrow(.), fold = occ$fold) %>% 
            rename(x = X, y = Y) %>% 
            select(id, sp, x, y, fold)
    }, error = function(e){
        set.seed(seed)
        bg %>% sample_n(size = nrow(occ), replace = FALSE) %>% 
            mutate(id = 1:nrow(.), fold = occ$fold) %>% 
            select(id, sp, x, y, fold)
    })
    
    # Save out
    write.csv(occ, file.path(sp_dir, sprintf("p_%s.csv", sp)), 
              row.names = FALSE)
    write.csv(bg_eval, file.path(sp_dir, sprintf("bg_eval_%s.csv", sp)), 
              row.names = FALSE)
    write.csv(bg, file.path(sp_dir, sprintf("bg_%s.csv", sp)), 
              row.names = FALSE)
    
    ########### Extract environmental values ###########
    occ <- extract(vars, occ %>% select(x, y), ID = FALSE) %>% 
        cbind(occ %>% select(fold)) %>% mutate(occ = 1)
    bg <- extract(vars, bg %>% select(x, y), ID = FALSE) %>% 
        mutate(occ = 0)
    bg_eval <- extract(vars, bg_eval %>% select(x, y), ID = FALSE) %>% 
        cbind(bg_eval %>% select(fold)) %>% mutate(occ = 0)
    
    ########### SDM ###########
    # Cross validation
    cv_sdm <- lapply(sort(unique(occ$fold)), function(fld){
        # Subset data
        training <- occ %>% filter(fold != fld) %>% select(-fold) %>% 
            rbind(bg)
        # Evaluation
        testing <- occ %>% filter(fold == fld) %>% select(-fold) %>% 
            rbind(bg_eval %>% filter(fold == fld) %>% select(-fold))
        
        # convert response to factor for classification
        training$occ <- as.factor(training$occ)
        testing$occ <- as.factor(testing$occ)
        
        # calculate sub-samples
        prNum <- sum(training$occ == 1) # number of presence records
        spsize <- c("0" = prNum, "1" = prNum) # sample size for both classes
        
        # RF with down-sampling
        set.seed(seed)
        mod <- randomForest(
            occ ~ .,
            data = training,
            ntree = 1000,
            sampsize = spsize,
            replace = TRUE)
        
        # Evaluation
        pred_training <- predict(mod, training, type = "prob", index = 2)[, 2]
        pred_testing <- predict(mod, testing, type = "prob", index = 2)[, 2]
        pred_br_testing <- predict(mod, testing)
        
        eval_training <- mod_eval(pred_training, mod$predicted, training$occ)
        eval_testing <- mod_eval(pred_testing, pred_br_testing, testing$occ)
        
        rbind(eval_training %>% mutate(type = "training"),
              eval_testing %>% mutate(type = "testing")) %>% 
            mutate(fold = fld)
    }) %>% bind_rows()
    
    write.csv(cv_sdm, file.path(sp_dir, sprintf("cv_eval_%s.csv", sp)),
              row.names = FALSE)
    
    # Full model
    training <- occ %>% select(-fold) %>% rbind(bg)
    training$occ <- as.factor(training$occ)
    # calculate sub-samples
    prNum <- sum(training$occ == 1) # number of presence records
    spsize <- c("0" = prNum, "1" = prNum) # sample size for both classes
    
    # RF with down-sampling
    set.seed(seed)
    mod <- randomForest(
        occ ~ .,
        data = training,
        ntree = 1000,
        sampsize = spsize,
        importance = TRUE, 
        localImp = TRUE,
        replace = TRUE)
    save(mod, file = file.path(sp_dir, sprintf("rf_dws_%s.rda", sp)))
    
    ########### Prediction ###########
    # Baseline
    # Crop the variables
    msk <- range %>% st_buffer(90000) # 1km/yr until 2100
    vars <- mask(crop(vars, msk), msk)
    
    pred_current <- predict(vars, mod, type = "prob", index = 2)
    names(pred_current) <- "baseline"
    writeRaster(pred_current, file.path(sp_dir, sprintf("pred_base_%s.tif", sp)))
    
    # Future
    fnames <- list.files(file.path(var_dir, "OtherEnv"), full.names = TRUE)
    preds_future <- do.call(c, lapply(fnames, function(fname){
        lyr <- rast(fname) %>% crop(., msk) %>% mask(., msk)
        lyr <- subset(lyr, var_list)
        predict(lyr, mod, type = "prob", index = 2)
    }))
    names(preds_future) <- gsub(".tif", "", basename(fnames))
    
    writeRaster(preds_future, file.path(sp_dir, sprintf("pred_scenarios_%s.tif", sp)))
}