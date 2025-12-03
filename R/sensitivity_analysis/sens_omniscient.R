# Load some markdown required packages
library(sf)
library(terra)
library(dplyr)
library(fastshap)
library(randomForest)
library(pROC)
library(virtualspecies)
library(doParallel)

root_dir <- "/scratch/ls1686/SClimpact"
occ_dir <- file.path(root_dir, "data/occurrences")
var_dir <- file.path(root_dir, "data/variables")
range_dir <- file.path(root_dir, "data/IUCN/Expert_Maps")
ncores <- 64

source("/home/ls1686/SClimpact/R/utils.R")

# Get a list of virtual species
species_list <- read.csv(
    file.path(root_dir, "results/species_reliable_final.csv"))

species_list <- lapply(1:nrow(species_list), function(i){
    r <- species_list[i, ]
    
    var_list <- read.csv(
        file.path(var_dir, "variable_list", sprintf("%s.csv", r$species)))
    if (length(na.omit(var_list$var_uncorrelated)) == 1){
        var_list <- na.omit(var_list$var)
    } else var_list <- na.omit(var_list$var_uncorrelated)
    
    r %>% mutate(var_n = length(var_list))
}) %>% bind_rows()

# Just a random middle size range species as the template
set.seed(123)
sps <- species_list %>%
    filter(var_n <= 3) %>% 
    sample_n(size = 100)

# Omniscient modeling
dst_dir <- file.path(root_dir, "results/validate")

evals <- lapply(sps$species, function(sp){
    message(sp)
    seed <- sum(utf8ToInt(sp))
    
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
    
    # Calculate Shapley values as a omniscient case
    ## Define the wrapper function
    pfun <- function(X.model, newdata) {
        pred_equation(X.model, newdata)}
    
    # Get new data
    dat <- values(vars_sp) %>% na.omit() %>% data.frame()
    
    set.seed(seed)
    dat_less <- dat[sample(1:nrow(dat), nrow(bg_smaller_vals)), ]
    
    items <- list("dat" = dat, "dat_less" = dat_less,
                  "bg_bigger" = bg_vals, "bg_smaller" = bg_smaller_vals)
    
    results <- lapply(c("dat", "dat_less", "bg_bigger", "bg_smaller"), function(x){
        train <- items[[x]]
        
        # Calculate Shapley values
        registerDoParallel(cores = ncores)
        set.seed(seed, "L'Ecuyer-CMRG")
        
        shap_explain <- explain(
            v_sp,
            X = train,
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
        ## Note that the suitability will be rescaled to match with species
        ## But this does not mean 0.5 is necessarily the actual threshold
        ## This is why background samples are important.
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
        r1 <- shaps >= 0 # Using Shapley values
        r2 <- suits >= avgs_bg # Using the defined domain average
        r <- r1 + r2
        
        # 2) sum up per layer (assuming equal-area pixels)
        intersect_sum <- global(r == 2, "sum", na.rm = TRUE)[, 1]
        union_sum <- global(r >= 1, "sum", na.rm = TRUE)[, 1]
        
        # 3) band-wise IoU
        iou <- mean(intersect_sum / union_sum)
        
        data.frame(item = x, r = r_ss, iou = iou)
    }) %>% bind_rows() %>% mutate(sp = sp)
}) %>% bind_rows()

write.csv(evals, file.path(dst_dir, "evals_raw.csv"), row.names = FALSE)
