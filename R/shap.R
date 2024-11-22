#' @title shap
#' @description Calculate Shapley values for baseline and future 
#' based on scenarios.
#' @param sp (`character`) The species name with format Name_name.
#' @param work_dir (`character`) The directory to save load the modeling 
#' results and save new results.
#' @param var_dir (`character`) The directory for environmental variables.
#' @param range_dir (`character`) The directory for IUCN ranges.
#' @param max_nshap (`integer`) The maximum nsim of Monte Carlo simulation
#' for `fastshap` to calculate Shapely values. 
#' @param seed (`integer`) The seed for randomization. 
#' @return All results are saved to `work_dir` under the folder with same name
#' as sp.
#' 
#' @import randomForest
#' @import terra
#' @import sf
#' @import dplyr
#' @import intervals
#' @import blockCV
#' @import fastshap
#' @import doParallel
#' @export
#' @examples
#' \donttest{
#' sp <- "Akodon_albiventer"
#' work_dir <- "results/sdm"
#' var_dir <- "data/variables"
#' range_dir <- "data/IUCN/Expert_Maps"
#' shap(sp, work_dir, var_dir, range_dir)
#' }
#' 

shap <- function(sp,
                 work_dir = "results/sdm",
                 var_dir = "data/variables",
                 range_dir = "data/IUCN/Expert_Maps",
                 max_nshap = 10000, # Set 1e*
                 seed = 123,
                 ncores = 10){ 
    if (!dir.exists(work_dir)) stop("No such folder to work with.")
    sp_dir <- file.path(work_dir, sp)
    if (!dir.exists(sp_dir)) stop("Run model for this species first.")
    
    # Load the model
    load(file.path(sp_dir, sprintf("rf_dws_%s.rda", sp)))
    
    # Load training
    occ <- read.csv(file.path(sp_dir, sprintf("p_%s.csv", sp)))
    bg <- read.csv(file.path(sp_dir, sprintf("bg_%s.csv", sp)))
    bg_eval <- read.csv(file.path(sp_dir, sprintf("bg_eval_%s.csv", sp)))
        
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
    
    # Load dispersal rate and continents
    dispersal_rate <- read.csv(
        file.path(var_dir, "species_dispersal_rate.csv")) %>% 
        filter(species == gsub("_", " ", sp)) %>% pull(dispersal_rate)
    continents <- st_read(
        file.path(var_dir, "continents.geojson"), quiet = TRUE) %>% 
        slice(unique(unlist(st_intersects(range, .))))
    
    # Reconstruct training
    occ <- extract(vars, occ %>% select(x, y), ID = FALSE) %>% 
        cbind(occ %>% select(fold)) %>% mutate(occ = 1) %>% na.omit()
    bg <- extract(vars, bg %>% select(x, y), ID = FALSE) %>% 
        mutate(occ = 0) %>% na.omit()
    bg_eval <- extract(vars, bg_eval %>% select(x, y), ID = FALSE) %>% 
        mutate(occ = 0) %>% na.omit()
    training <- occ %>% select(-fold) %>% rbind(bg) %>% 
        select(-occ)
    testing <- occ %>% select(-fold) %>% rbind(bg_eval) %>% 
        select(-occ)
    if (nrow(testing) > 3000){
        set.seed(seed)
        testing <- testing %>% sample_n(size = 3000)
    }
    rm(occ, bg, bg_eval)
    
    ##############
    # 1. Define pfun for shap
    pfun <- function(X.model, newdata) {
        predict(X.model, newdata, type = "prob")[, 2]}
    
    # 2. Diagnose the enough max_nshap
    ## Make a reasonable list of nsim to test
    nshap_list <- 10^(1:(nchar(as.integer(max_nshap)) - 1))
    if (any(nshap_list > 1000) & max(nshap_list) > 1000){
        nshap_list <- c(nshap_list, seq(1000, max(nshap_list), 3000)) %>% 
            unique() %>% sort()
    }
    shap_list <- lapply(nshap_list, function(nsim){
        registerDoParallel(cores = min(ncores, ncol(testing)))
        set.seed(seed, "L'Ecuyer-CMRG")
        
        shap_explain <- fastshap::explain(
            mod,
            X = training,
            newdata = testing,
            nsim = nsim,
            adjust = TRUE, parallel = TRUE,
            pred_wrapper = pfun)
        
        stopImplicitCluster()
        
        shap_explain
        })
    
    shap_cors <- lapply(1:(length(shap_list) - 1), function(i){
        data.frame(nshap = nshap_list[i],
            cor = mean(diag(cor(shap_list[[i]], 
                                shap_list[[length(shap_list)]]))))
    }) %>% bind_rows()
    shap_cors <- rbind(
        shap_cors, data.frame(nshap = nshap_list[length(nshap_list)],
                              cor = 1))
    
    write.csv(shap_cors, file.path(sp_dir, sprintf("shap_cor_%s.csv", sp)))
    
    # If the increase is less than 0.01, then stop
    nshap <- shap_cors %>% 
        mutate(gain = c(shap_cors$cor[2:nrow(.)] - 
                            shap_cors$cor[1:(nrow(.) - 1)], 0)) %>% 
        filter(gain < 0.001) %>% slice(1) %>% pull(nshap)
    
    # Now start to calculate the shapely values
    # Mask out
    dispersal_distance <- dispersal_rate * (2100 - 2011 + 1) * 1000
    msk <- range %>% st_buffer(dispersal_distance) %>% 
        st_intersection(continents)
    
    # Baseline
    fn <- file.path(sp_dir, sprintf("shap_base_%s_%s.tif", sp, nshap))
    if (!file.exists(fn)){
        vars <- mask(crop(vars, msk), msk)
        new_data <- values(vars) %>% na.omit() %>% data.frame()
        
        # The real calculation
        registerDoParallel(cores = min(ncores, ncol(new_data)))
        set.seed(seed, "L'Ecuyer-CMRG")
        
        shap_explain <- fastshap::explain(
            mod, X = training,
            newdata = new_data,
            nsim = nshap,
            adjust = TRUE, parallel = TRUE,
            pred_wrapper = pfun)
        
        stopImplicitCluster()
        
        shap_base <- vars
        vals <- values(shap_base)
        vals[complete.cases(vals), ] <- shap_explain
        values(shap_base) <- vals
        writeRaster(shap_base, fn, overwrite = TRUE)
    }
    
    # Scenarios
    fnames <- list.files(file.path(var_dir, "OtherEnvMean"), full.names = TRUE)
    for (fname in fnames){
        scenario <- gsub(".tif", "", basename(fname))
        fn <- file.path(sp_dir, sprintf("shap_%s_%s_%s.tif", scenario, sp, nshap))
        
        if (!file.exists(fn)){
            vars <- rast(fname) %>% crop(., msk) %>% mask(., msk)
            vars <- subset(vars, var_list)
            new_data <- values(vars) %>% na.omit() %>% data.frame()
            
            # The real calculation
            registerDoParallel(cores = min(ncores, ncol(new_data)))
            set.seed(seed, "L'Ecuyer-CMRG")
            
            shap_explain <- fastshap::explain(
                mod, X = training,
                newdata = new_data,
                nsim = nshap,
                adjust = TRUE, parallel = TRUE,
                pred_wrapper = pfun)
            
            stopImplicitCluster()
            
            shap_scn <- vars
            vals <- values(shap_scn)
            vals[complete.cases(vals), ] <- shap_explain
            values(shap_scn) <- vals
            writeRaster(shap_scn, fn, overwrite = TRUE)
        }
    }
}
