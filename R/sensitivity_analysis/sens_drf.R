# Load some markdown required packages
library(sf)
library(terra)
library(dplyr)
library(fastshap)
library(randomForest)
library(pROC)
library(virtualspecies)
library(doParallel)
library(optparse)
source("R/utils.R")

root_dir <- "/scratch/ls1686/SClimpact"
occ_dir <- file.path(root_dir, "data/occurrences")
var_dir <- file.path(root_dir, "data/variables")
range_dir <- file.path(root_dir, "data/IUCN/Expert_Maps")
ncores <- 64

# Get a list of virtual species
sps <- read.csv(
    file.path(root_dir, "results/validate/evals_raw.csv")) %>% 
    filter(item == "dat") %>% 
    filter(r >= 0.95 & iou >= 0.95)

# Omniscient modeling
dst_dir <- file.path(root_dir, "results/validate")

# Random forest modeling
rf_evals <- lapply(sps$sp, function(sp){
    message(sp)
    seed <- sum(utf8ToInt(sp))
    
    ########### Create folder for results ###########
    if (!dir.exists(dst_dir)) dir.create(dst_dir, recursive = TRUE)
    sp_dir <- file.path(dst_dir, sp)
    if (dir.exists(sp_dir)) unlink(sp_dir, recursive = TRUE)
    dir.create(sp_dir)
    
    ########### Load data ###########
    ## background, bigger
    bg <- read.csv(file.path(occ_dir, "bg", sprintf("%s.csv", sp)))
    
    ## Environmental variables
    vars <- rast(file.path(var_dir, "Env/AllEnv.tif"))
    var_list <- read.csv(
        file.path(var_dir, "variable_list", sprintf("%s.csv", sp)))
    if (length(na.omit(var_list$var_uncorrelated)) == 1){
        var_list <- na.omit(var_list$var)
    } else var_list <- na.omit(var_list$var_uncorrelated)
    vars <- subset(vars, var_list)
    
    bg_vals <- extract(vars, bg[, c("x", "y")], ID = FALSE) %>% na.omit()
    
    # Load range polygon
    range <- st_read(
        file.path(range_dir, sprintf("%s.geojson", sp)), 
        quiet = TRUE) %>% st_transform(crs(vars))
    
    vars_sp <- mask(crop(vars, range), range) %>% trim()
    
    # background, smaller
    buf <- min(st_bbox(range)[3] - st_bbox(range)[1], 
               st_bbox(range)[4] - st_bbox(range)[2]) / 8
    range_smaller <- st_buffer(range, -buf)
    
    while (all(st_is_empty(range_smaller))){
        buf <- buf / 2
        range_smaller <- st_buffer(range, -buf)}
    
    set.seed(seed)
    bg_smaller_vals <- spatSample(
        mask(crop(vars, range_smaller), range_smaller) %>% trim(),
        size     = nrow(bg_vals),
        method   = "random",
        as.points = FALSE,
        na.rm  = TRUE)
    
    nms <- names(vars_sp)
    
    # Save out the importance for linear
    formula <- paste(nms, collapse = " + ")
    
    # Step 5: Randomly generate a virtual species with automatic
    v_sp <- generateRandomSp(vars_sp, formula = formula, seed = seed)
    
    var_scale_prs <- var_rescale_parameters(v_sp, vars_sp)
    v_sp$details$rescale.each.response.parameters <- var_scale_prs
    
    prev <- prevalence(v_sp$pa.raster)
    
    # randomly select between 20 to 200 points for presences
    ## assume it is a perfect detection
    np <- sample(20:200, 1, prob = dpois(x = 1:181, lambda = round(prev * 181)))
    thred <- mean(values(v_sp$suitab.raster), na.rm = TRUE)
    v_sp$pa.raster <- v_sp$suitab.raster >= thred
    v_sp$pa.raster[v_sp$pa.raster == 0] <- NA
    
    set.seed(seed)
    presence.points <- spatSample(
        v_sp$pa.raster, size = np, method = "stratified", 
        na.rm = TRUE, xy = TRUE,
        weights = v_sp$suitab.raster)
    
    pres <- extract(vars, presence.points[, c("x", "y")], ID = FALSE)
    
    # Calculate Shapley values as a omniscient case
    ## Define the wrapper function for random forest
    pfun <- function(X.model, newdata) {
        predict(X.model, newdata, type = "prob")[, 2]}
    
    # Get new data
    dat <- values(vars_sp) %>% na.omit() %>% data.frame()
    
    set.seed(seed)
    dat_less <- dat[sample(1:nrow(dat), nrow(dat) / 2), ]
    
    items <- list("dat" = dat, "dat_less" = dat_less,
                  "bg_bigger" = bg_vals, "bg_smaller" = bg_smaller_vals)
    
    results <- lapply(c("dat", "dat_less", "bg_bigger", "bg_smaller"), function(x){
        bg <- items[[x]]
        
        training <- rbind(pres %>% mutate(occ = 1),
                          bg %>% mutate(occ = 0))
        training$occ <- as.factor(training$occ)
        
        #### Training ##########
        # calculate sub-samples
        prNum <- sum(training$occ == 1) # number of presence records
        bgNum <- sum(training$occ == 0)
        if (bgNum < prNum) prNum <- bgNum
        spsize <- c("0" = prNum, "1" = prNum) # sample size for both classes
        
        # RF with down-sampling
        set.seed(seed)
        mod <- randomForest(
            occ ~ .,
            data = training,
            ntree = 500,
            sampsize = spsize,
            importance = TRUE, 
            localImp = TRUE,
            replace = TRUE)
        
        save(mod, file = file.path(sp_dir, sprintf("rf_dws_%s.rda", sp)))
        
        pred <- predict(mod, dat, type = "prob", index = 2)[, 2]
        reals <- values(v_sp$suitab.raster) %>% na.omit() %>% as.vector()
        mod_cor <- cor(pred, reals)
        
        #### SHAP ##########
        # Calculate Shapley values
        registerDoParallel(cores = ncores)
        set.seed(seed, "L'Ecuyer-CMRG")
        
        shap_explain <- fastshap::explain(
            mod, X = training %>% select(-occ),
            newdata = dat,
            nsim = 1000,
            adjust = TRUE, 
            parallel = TRUE,
            pred_wrapper = pfun)
        
        stopImplicitCluster()
        
        shaps <- vars_sp
        vals <- values(shaps)
        vals[complete.cases(vals), ] <- shap_explain
        values(shaps) <- vals
        
        # Calculate suitability of each variables
        parameters <- v_sp$details$parameters
        minmaxs <- v_sp$details$rescale.each.response.parameters
        
        suit_vars <- lapply(1:length(parameters), function(i) {
            nm <- names(parameters)[i]
            cur.seq <- dat[[nm]]
            values <- do.call(match.fun(parameters[[i]]$fun), 
                              args = c(list(cur.seq), parameters[[i]]$args))
            
            # Rescale
            minmax <- minmaxs %>% filter(var == nm)
            data.frame(x =  (values - minmax$min) / (minmax$max - minmax$min)) %>% 
                rename("{nm}" := x)
        }) %>% bind_cols()
        
        suits <- vars_sp
        vals <- values(suits)
        vals[complete.cases(vals), ] <- as.matrix(suit_vars %>% select(names(shaps)))
        values(suits) <- vals
        
        # Check their correlation
        r_ss <- mean(diag(cor(shap_explain, suit_vars %>% select(names(shaps)))))
        
        # Compare Shapley values and suitability of each variables
        ## If they match, that means Shapley values are effective to detect variable contribution
        suits_bg <- lapply(1:length(parameters), function(i) {
            nm <- names(parameters)[i]
            cur.seq <- train[[nm]]
            values <- do.call(match.fun(parameters[[i]]$fun), 
                              args = c(list(cur.seq), parameters[[i]]$args))
            
            # Rescale
            minmax <- minmaxs %>% filter(var == nm)
            data.frame(x =  (values - minmax$min) / (minmax$max - minmax$min)) %>% 
                rename("{nm}" := x)
        }) %>% bind_cols()
        
        avgs_bg <- suits_bg %>% select(names(shaps)) %>% colMeans()
        
        # 1) intersection and union per layer
        r1 <- shaps >= 0
        r2 <- suits >= avgs_bg
        r <- r1 + r2
        
        # 2) sum up per layer (assuming equal-area pixels)
        intersect_sum <- global(r == 2, "sum", na.rm = TRUE)[, 1]
        union_sum <- global(r >= 1, "sum", na.rm = TRUE)[, 1]
        
        # 3) band-wise IoU
        iou <- mean(intersect_sum / union_sum)
        
        data.frame(item = x, r = r_ss, iou = iou)
    }) %>% bind_rows() %>% mutate(sp = sp)
}) %>% bind_rows()

write.csv(rf_evals, file.path(dst_dir, "evals_rf_raw.csv"), row.names = FALSE)
