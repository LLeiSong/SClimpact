# Get the variable importance for a BART model
var_imp <- function(model){
    # Extract importance from model result
    if(class(model) == 'rbart') {
        basenames <- unlist(attr(model$fit[[1]]$data@x,"drop"))
        names <- names(which(basenames == FALSE))
        importance <- rowMeans(model$varcount/colSums(model$varcount))
        fitobj <- model$fit[[1]]
    }
    if(class(model) == 'bart') {
        basenames <- unlist(attr(model$fit$data@x,"drop"))
        names <- names(which(basenames == FALSE))
        importance <- colMeans(model$varcount/rowSums(model$varcount))
        fitobj <- model$fit
    }
    
    var.df <- data.frame(names, importance)
    
    # Set importance to 0 for dropped variables
    missing <- attr(fitobj$data@x, "term.labels")[
        !(attr(fitobj$data@x, "term.labels") %in% 
              names(unlist(attr(fitobj$data@x, "drop"))))]
    
    if(length(missing) > 0) {
        missing.df <- data.frame(names = missing, importance = 0)
        var.df <- rbind(var.df, missing.df)
    }
    
    var.df$names <- factor(var.df$names)
    var.df <- transform(var.df, names = reorder(names, -importance))
    
    return(var.df)
}

quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}

# The main function to diagnose the variable importance
var_diagnose <- function(bart.x, 
                         bart.y,
                         n.trees = 10, 
                         iter = 200){
    # Check data structure
    if (nrow(bart.x) != length(bart.y)){
        stop("Not the same length of dependent and independent variables.")}
    
    # Clean the data
    comp <- complete.cases(bart.x) & complete.cases(bart.y)
    bart.x <- bart.x[comp, ]
    bart.y <- bart.y[comp]
    
    quietly(mod_base <- bart(y.train = bart.y, x.train = bart.x,
                             ntree = 200, keeptrees = TRUE))
    
    # Get the keeped variable names
    if(class(mod_base) == 'rbart') {
        mod_obj <- mod_base$fit[[1]]
    } else if(class(mod_base) == 'bart') {
        mod_obj <- mod_base$fit
    }
    
    nms_keep <- names(which(unlist(attr(mod_obj$data@x,"drop")) == FALSE))
    bart.x <- bart.x %>% select(all_of(nms_keep))
    
    # Start the diagnose the variables
    nvars <- ncol(bart.x)
    varnums <- c(1:nvars)
    varlist.orig <- varlist <- colnames(bart.x)
    
    rmses <- data.frame(Variable.number = c(), RMSE = c())
    dropped.varlist <- c()
    
    for(j in c(nvars:3)) {
        rmse.list <- c()
        for(i in 1:iter) {
            quietly(
                mod_j <- bart(y.train = bart.y, x.train = bart.x[, varnums],
                              ntree = n.trees, keeptrees = TRUE))
            
            vi.j <- var_imp(mod_j)
            if(i == 1) {
                vi.j.df <- vi.j
            } else {
                vi.j.df[, i + 1] <- vi.j[, 2]
            }
            
            pred.p <- colMeans(pnorm(mod_j$yhat.train))[bart.y == 1]
            pred.a <- colMeans(pnorm(mod_j$yhat.train))[bart.y == 0]
            rmsej.i <- Metrics::rmse(
                c(pred.p, pred.a), # pred
                c(rep(1, length(pred.p)), rep(0, length(pred.a)))) # true
            rmse.list <- c(rmse.list, rmsej.i)
        }
        
        vi.j <- data.frame(
            names = vi.j.df[, 1],
            importance = rowMeans(vi.j.df[, -1])) %>% 
            arrange(importance)
        
        drop.var <- vi.j[1, 1]
        dropped.varlist <- c(dropped.varlist, as.character(drop.var))
        
        rmsej <- mean(rmse.list)
        
        rmses <- rbind(rmses, c(nvars - j, rmsej))
        colnames(rmses) <- c('VarsDropped','RMSE')
        
        varnums <- varnums[!(varnums == which(varlist.orig == drop.var))]
        varlist <- varlist.orig[varnums]
    }
    
    varlist.final <- varlist.orig[
        !(varlist.orig %in% dropped.varlist[
            0:(which(rmses$RMSE == min(rmses$RMSE)) - 1)])]
    return(varlist.final)
}

# Add buffer into occurrences
buf_polygon <- function(occ, buffer_dist = NA){
    # Convert to terra
    occ <- vect(occ)
    
    # Check
    if(is.na(terra::linearUnits(occ))){
        stop("The spatial units of the data cannot be found!")}
    
    if (is.na(buffer_dist)) {
        # Uses the 95% quantile of the minimum distance between each point
        Distance <- terra::distance(occ)
        Distance <- as.matrix(Distance)
        mindist <- c()
        for (q in 1:ncol(Distance)) {
            DistanceZero <- Distance[which(Distance[, q] > 0), q]
            mindist <- c(mindist, min(DistanceZero))}
        BufferWidth <- 2 * stats::quantile(mindist, 0.95)
        rm(Distance, mindist); gc()
    } else {
        BufferWidth <- buffer_dist
    }
    
    # Makes a buffer (width dependent) around the first occurrence point
    combinedPolygon <- buffer(occ, BufferWidth)
    combinedPolygon <- aggregate(combinedPolygon)
    
    return(combinedPolygon)
}

# Combine both range and buffered occurrences
make_domain <- function(occ, range, buffer_dist = NA){
    # Convert to terra
    occ <- vect(occ)
    range <- vect(range)
    
    # Check
    if (crs(occ) != crs(range)){
        range <- project(range, crs(occ))
    }
    
    if(is.na(terra::linearUnits(occ))){
        stop("The spatial units of the data cannot be found!")}
    
    if (is.na(buffer_dist)) {
        # Uses the 95% quantile of the minimum distance between each point
        err_detect <- try({
            Distance <- terra::distance(occ)
            Distance <- as.matrix(Distance)
            mindist <- c()
            for (q in 1:ncol(Distance)) {
                DistanceZero <- Distance[which(Distance[, q] > 0), q]
                mindist <- c(mindist, min(DistanceZero))}
            BufferWidth <- 2 * stats::quantile(mindist, 0.95)
            rm(Distance, mindist); gc()
        }, silent = TRUE)
        
        if (inherits(err_detect, "try-error")) BufferWidth <- 10000 # one pixel
        
    } else {
        BufferWidth <- buffer_dist
    }
    
    # Makes a buffer (width dependent) around the first occurrence point
    occ_buf <- buffer(occ, BufferWidth)
    range_buf <- buffer(range, BufferWidth)
    combinedPolygon <- rbind(occ_buf, range_buf)
    combinedPolygon <- aggregate(combinedPolygon)
    
    return(combinedPolygon)
}

# Function to sample occurrences using environmental bins
varela_sample_occ <- function (dat, no_bin = 25) {
    out_ptz <- dat[, 1:2]
    for(i in 3:length(names(dat))) {
        # Cut into bins
        k <- dat[!is.na(dat[, i]), i]
        rg <- range(k)
        res <- (rg[2] - rg[1]) / no_bin
        
        # Calculate value
        d <- (dat[, i] - rg[1]) / res
        f <- ceiling(d)
        f[f == 0] <- 1
        
        # Collect result and set name
        names(f) <- names(dat)[i]
        out_ptz <- cbind(out_ptz, f)
        names(out_ptz)[length(names(out_ptz))] <- names(dat)[i]
    }
    
    # subsample the bin membership df to come up with the filled bins
    sub_ptz <- dplyr::distinct(out_ptz[, -1:-2])
    # count the number of filled bins
    no_grps <- nrow(sub_ptz)
    # add a column with the group membership number; this number is arbitrary
    sub_ptz$grp <- c(1:no_grps)
    out_ptz <- suppressMessages(left_join(out_ptz, sub_ptz))
    
    final_out <- lapply(1:no_grps, function(i){
        grp_mbrs <- out_ptz[out_ptz$grp == i, c(1, 2)]
        grp_mbrs[sample(1:nrow(grp_mbrs), 1), ]
    }) %>% bind_rows()
    
    final_out <- merge(final_out, dat, 
                       by = c("x", "y"), all.x = TRUE)
    
    # Return result
    return(final_out)
}

# Make the buffer polygon for background sampling
# Assume the projection is in meter
background_sampling <- function(sp_name, 
                                occ_sf, 
                                vars,
                                dist = NA, # dist in meter
                                nbg = 5000, 
                                spatial_weight = 0.5){
    # Make template
    template <- vars[[1]]
    combinedPolygon <- buf_polygon(occ_sf)
    
    # Load range map and AOH
    range <- vect(file.path(
        "data/IUCN/Expert_Maps",
        sprintf("%s.geojson", gsub(" ", "_", sp_name)))) %>% 
        project(crs(template))
    
    # Because MaxEnt work with background not pseudo-samples
    # so the objective here is to provide a possible range, rather than the 
    # possible set of absence.
    range <- terra::union(buffer(range, BufferWidth), convHull(combinedPolygon))
    
    # Get the numbers
    nbgBuff <- round(nbg * spatial_weight)
    nbgFull <- round(nbg - nbgBuff)
    
    # Sampling from the full training area--------------------
    full_map <- mask(crop(template, range), range)
    
    if (nbgFull * 10 > terra::ncell(full_map)) {
        RastFullBG <- terra::ncell(full_map)
    } else {
        RastFullBG <- nbgFull * 10
    }
    
    full_vars <- mask(crop(vars, full_map), full_map)
    RandomTrain <- spatSample(full_vars, RastFullBG, na.rm = TRUE, xy = TRUE)
    nbinscount <- nbins
    
    NUniqueClim <- nrow(unique(RandomTrain[, 3:ncol(RandomTrain)]))
    
    if (NUniqueClim < nbgFull) {
        nbgFull <- NUniqueClim
        nbinscount <- 98
    }
    
    FullPointsEnv <- varela_sample_occ(RandomTrain, nbinscount)
    
    if (nrow(FullPointsEnv) >= nbgFull) {
        while(nrow(FullPointsEnv) > nbgFull && nbinscount > 2) {
            nbinscount <- nbinscount - 1
            FullPointsEnv <- varela_sample_occ(RandomTrain, nbinscount)
        }
        
        FullPointsEnv2 <- varela_sample_occ(RandomTrain, (nbinscount + 1))
        diff1 <- abs(nrow(FullPointsEnv) - nbgFull)
        diff2 <- abs(nrow(FullPointsEnv2) - nbgFull)
        
        if (diff2 < diff1) {
            FullPointsEnv <- FullPointsEnv2
        }
    } else {
        while(nrow(FullPointsEnv) < nbgFull && nbinscount < max(nbins, 99)) {
            nbinscount <- nbinscount + 1
            FullPointsEnv <- varela_sample_occ(RandomTrain, nbinscount)
        }
        
        FullPointsEnv2 <- varela_sample_occ(RandomTrain, (nbinscount - 1))
        
        diff1 <- abs(nrow(FullPointsEnv) - nbgFull)
        diff2 <- abs(nrow(FullPointsEnv2) - nbgFull)
        
        if (diff2 < diff1) {
            FullPointsEnv <- FullPointsEnv2
        }
    }
    
    FullPointsEnv <- FullPointsEnv[stats::complete.cases(FullPointsEnv), ]
    
    # Sampling from the buffer area------------------
    buff_map <- mask(crop(template, combinedPolygon), combinedPolygon)
    
    if (nbgFull * 10 > terra::ncell(buff_map)) {
        RastFullBG <- terra::ncell(buff_map)
    } else {
        RastFullBG <- nbgFull * 10
    }
    
    #randomly samples nbg*5 number of points from the shapefile
    buff_vars <- mask(crop(vars, buff_map), buff_map)
    
    BuffClim <- spatSample(buff_vars, RastFullBG, na.rm = TRUE, xy = TRUE)
    nbinscount <- nbins
    
    #Prints a warning if the number of unique climates is less than the buffer background points wanted
    NUniqueClim <- nrow(unique(BuffClim))
    
    if (NUniqueClim < nbgFull) {
        nbgFull <- NUniqueClim
        nbinscount <- 98
    }
    
    if (NUniqueClim < nbgBuff) {
        nbgBuff <- NUniqueClim
        nbinscount <- 98
    }
    
    BuffPointsEnv <- varela_sample_occ(BuffClim, nbinscount)
    
    if (nrow(BuffPointsEnv) >= nbgBuff) {
        while(nrow(BuffPointsEnv) > nbgBuff && nbinscount > 2) {
            nbinscount <- nbinscount - 1
            BuffPointsEnv <- varela_sample_occ(BuffClim, nbinscount)
        }
        
        BuffPointsEnv2 <- varela_sample_occ(BuffClim, (nbinscount + 1))
        
        diff1 <- abs(nrow(BuffPointsEnv) - nbgBuff)
        diff2 <- abs(nrow(BuffPointsEnv2) - nbgBuff)
        
        if (diff2 < diff1) {
            BuffPointsEnv <- BuffPointsEnv2
        }
    } else {
        
        while(nrow(BuffPointsEnv) < nbgBuff && nbinscount < max(nbins, 99)) {
            nbinscount <- nbinscount + 1
            BuffPointsEnv <- varela_sample_occ(BuffClim, nbinscount)
        }
        
        BuffPointsEnv2 <- varela_sample_occ(BuffClim, (nbinscount - 1))
        
        diff1 <- abs(nrow(BuffPointsEnv) - nbgBuff)
        diff2 <- abs(nrow(BuffPointsEnv2) - nbgBuff)
        
        if (diff2 < diff1) {
            BuffPointsEnv <- BuffPointsEnv2
        }
    }
    
    BuffPointsEnv <- BuffPointsEnv2[stats::complete.cases(BuffPointsEnv2), ]
    
    if (spatial_weight == 0) {
        background <- FullPointsEnv
        BuffPointsEnv <- FullPointsEnv[c(), ]
    } else if (spatial_weight == 1) {
        background <- BuffPointsEnv
        FullPointsEnv <- BuffPointsEnv[c(), ]
    } else {
        background <- rbind(FullPointsEnv, BuffPointsEnv)
    }
    
    background <- data.frame("Species" = rep(spplist[s], nrow(background)), background)
}

# Remove correlated variables
var_remove_cor <- function(env,
                           preferedPredictors,
                           tieRule = 'last',
                           corThreshold = .7) {
    if (!is.null(preferedPredictors)) {
        keepers <- which(names(env) %in% preferedPredictors)
        if (length(keepers) > 0) {
            env <- cbind(env[keepers], env[-keepers])
        }
    }
    c1 <- suppressWarnings(cor(env, use = 'complete.obs'))
    tossed <- NULL
    bad <- which(apply(c1, 1, function(x) {
        all(is.na(x)) | all(is.nan(x)) | all(is.null(x)) | sum(x, na.rm = T) == 1
    }))
    
    if (length(bad) > 0) {
        tossed = c(tossed, names(bad))
        c1.index = which(colnames(c1) %in% tossed)
        c1 = c1[-c1.index, -c1.index]
    }
    # number of variables a variable is too correlated with
    too.cor <- apply(abs(c1) > corThreshold, 1, sum) - 1
    
    while (any(too.cor >= 1)) {
        most.cor = too.cor[too.cor == max(too.cor)]
        if (tieRule == 'last') {
            toss = tail(most.cor, 1)
        } else if (tieRule == 'random') {
            toss = sample(1:length(most.cor), 1)
        }
        tossed = c(tossed, names(toss))
        c1.index = which(colnames(c1) %in% names(toss))
        c1 = c1[-c1.index, -c1.index]
        if (length(c1) == 1){
            if (c1 == 1) break
        } else too.cor = apply(abs(c1) > corThreshold, 1, sum) - 1
    }
    if (!is.null(tossed)) env = env[-which(names(env) %in% tossed)]
    names(env)
}

mod_eval <- function(score, pred, label){
    # TSS
    cm <- table(pred, label)
    if (all(rownames(cm) == 'TRUE')) {
        cm <- rbind('FALSE' = c(0, 0), cm)
    } else if (all(rownames(cm) == 'FALSE')) {
        cm <- rbind(cm, 'TRUE' = c(0, 0))}
    
    # Traditional ones: sensitivity (TPR), specificity (TNR), and TSS
    # TPR  = TP / (TP + FN)
    # TNR = TN / (TN + FP)
    tn <- cm[1]; fn <- cm[3]; tp <- cm[4]; fp <- cm[2]
    tpr <- cm[4] / sum(cm[4], cm[3])
    tnr <- cm[1] / sum(cm[1], cm[2])
    tss <- tpr + tnr - 1
    
    aucs <- evalmod(scores = score, labels = label)
    aucs <- auc(aucs)
    
    data.frame(tn = tn, fn = fn, tp = tp, fp = fp,
               tpr = tpr, tnr = tnr, tss = tss,
               auc = aucs$aucs[1], prc = aucs$aucs[2])
}

# Author:   Lei Song
# contact:  lsong@ucsb.edu
# Date:     2024-12-21
# Description: Generating virtual species

# An updated version of generateRandomSp to take formula
# Reference: virtualspecies R package: https://github.com/Farewe/virtualspecies
generateRandomSp <- function(
        raster.stack, 
        approach = "automatic",
        rescale = TRUE,
        convert.to.PA = TRUE,
        relations = c("gaussian", "linear", "logistic", "quadratic"),
        rescale.each.response = TRUE,
        realistic.sp = TRUE,
        species.type = "multiplicative",
        formula = NULL,
        niche.breadth = "any",
        sample.points = FALSE, 
        nb.points = 10000,
        PA.method = "probability",
        adjust.alpha = TRUE,
        beta = "random", # random, 0.5
        alpha = -0.1,
        species.prevalence = NULL,
        seed = 123,
        plot = TRUE){
    
    set.seed(seed)
    if(inherits(raster.stack, "Raster")) {
        raster.stack <- rast(raster.stack)
    }
    if(!(inherits(raster.stack, "SpatRaster")))
    {
        stop("raster.stack must be a SpatRaster object")
    }
    
    
    if(approach == "automatic")
    {
        if(nlyr(raster.stack) <= 5)
        {
            approach <- "response"
        } else
        {
            approach <- "pca"
        }
    } else if (approach == "random")
    {
        approach <- sample(c("response", "pca"), 1)
    } else if(!(approach %in% c("response", "pca")))
    {
        stop("Argument approach was misspecified. Either choose 'automatic', ",
             "'random', 'response' or 'pca'.")
    }
    
    var.names <- names(raster.stack)
    
    if(approach == "pca")
    {
        results <- generateSpFromPCA(raster.stack,
                                     rescale = rescale,
                                     niche.breadth = niche.breadth,
                                     sample.points = sample.points, 
                                     nb.points = nb.points,
                                     plot = FALSE)
    } else if (approach == "response")
    {
        parameters <- list()
        message(" - Determining species' response to predictor variables\n")
        if(any(!(relations %in% c("gaussian", "linear", "logistic", "quadratic"))))
        {
            stop(paste("Wrong relation type specified: pick among '", 
                       paste(c("gaussian", "linear", "logistic", "quadratic"), 
                             collapse = " "), "'",
                       collapse = " "))
        }
        valid.cells <- setValues(raster.stack[[1]], 1)
        var.order <- sample(var.names, length(var.names), replace = F)
        for (i in 1:length(var.order))
        {
            
            cur.var <- var.order[i]
            cur.rast <- raster.stack[[cur.var]]
            if(realistic.sp) cur.rast <- cur.rast * valid.cells # Cur.rast is 
            # here restricted to current suitable conds
            
            type <- sample(relations, 1)
            
            min_ <- global(cur.rast, "min", na.rm = TRUE)[1, 1]
            max_ <- global(cur.rast, "max", na.rm = TRUE)[1, 1]
            
            
            if (type == "gaussian")
            {
                parameters[[cur.var]] <- list(
                    fun = 'dnorm',
                    args = c(mean = sample(seq(min_,
                                               max_, 
                                               length = 100000), 1),
                             sd = sample(seq(0, 
                                             (max_ - min_), 
                                             length = 100000), 1))
                )
            } else if (type == "linear")
            { # At the moment this is not really useful because the rescale will 
                # transform the results in either 0:1 or 1:0, regardless of the slope
                # To be improved later
                parameters[[cur.var]] <- list(
                    fun = 'linearFun',
                    args = c(a = sample(seq(-1, 1, length = 100), 1),
                             b = sample(seq(min_, 
                                            max_, 
                                            length = 100000), 1))
                )
            } else if (type == "logistic")
            {
                beta.t <- sample(seq(min_,
                                     max_,
                                     length = 1000000), 1)
                alpha.t <-  sample(c(seq((max_ - min_)/1000,
                                         (max_ - min_)/100, length = 10),
                                     seq((max_ - min_)/100,
                                         (max_ - min_)/10, length = 100),
                                     seq((max_ - min_)/10,
                                         (max_ - min_)*10, length = 10)), size = 1)
                if(realistic.sp == TRUE)
                {
                    if(beta.t > max_)
                    {
                        alpha.t <- alpha.t
                    } else if (beta.t < min_)
                    {
                        alpha.t <- -alpha.t
                    } else
                    {
                        alpha.t <- sample(c(alpha.t, -alpha.t), 1)
                    }
                }
                
                parameters[[cur.var]] <- list(fun = 'logisticFun',
                                              args = c(alpha = alpha.t,
                                                       beta = beta.t)
                )
            } else if (type == "quadratic")
            {
                max.point <- sample(seq(min_,
                                        max_, 
                                        length = 1000), 1)
                a <- sample(seq(-.01, -20, length = 10000), 1)
                b <- - max.point * 2 * a
                parameters[[cur.var]] <- list(fun = 'quadraticFun',
                                              args = c(a = a,
                                                       b = b,
                                                       c = 0)
                )
                
            }
            
            # Restricting values to suitable conditions
            tmp.rast <- app(raster.stack[[cur.var]], fun = function(x)
            {
                do.call(match.fun(parameters[[cur.var]]$fun), 
                        args = c(list(x), parameters[[cur.var]]$args))
            }
            )
            tmp.rast <- (tmp.rast - global(tmp.rast, "min", na.rm = TRUE)[1, 1]) /
                (global(tmp.rast, "max", na.rm = TRUE)[1, 1] - 
                     global(tmp.rast, "min", na.rm = TRUE)[1, 1])
            valid.cells <- valid.cells * (tmp.rast > 0.05)
        }
        message(" - Calculating species suitability\n")
        results <- generateSpFromFun(raster.stack, 
                                     parameters, 
                                     rescale = rescale,
                                     formula = formula,
                                     species.type = species.type, 
                                     plot = FALSE,
                                     rescale.each.response = rescale.each.response)
    }
    
    if(convert.to.PA == TRUE) {
        message(" - Converting into Presence - Absence\n")
        
        # Need to adjust alpha to appropriate scale if rescale = FALSE
        if(rescale == FALSE) {
            if(adjust.alpha)
            {
                alpha <- diff(c(global(results$suitab.raster,
                                       min, na.rm = TRUE)[1, 1],
                                global(results$suitab.raster,
                                       max, na.rm = TRUE)[1, 1])) * alpha
            }
            results <- convertToPA(results, 
                                   PA.method = PA.method,
                                   alpha = alpha,
                                   beta = beta,
                                   species.prevalence = species.prevalence,
                                   plot = FALSE)
            
            if(plot) plot(results)
        } else {
            
            results <- convertToPA(results, 
                                   PA.method = PA.method,
                                   alpha = alpha,
                                   beta = beta,
                                   species.prevalence = species.prevalence,
                                   plot = FALSE)
            
            if(plot) plot(results)
        }
    } else {
        if(plot) plot(results)
    }
    
    return(results)
}

# Function to calculate the rescale parameters for variables
var_rescale_parameters <- function(sp, vars){
    # Reconstruct new data
    vars <- as.data.frame(vars) %>% na.omit()
    if (!is.null(sp$details$parameters)){
        parameters <- sp$details$parameters
        lapply(1:length(vars), function(i) {
            nm <- names(parameters)[i]
            cur.seq <- vars[[nm]]
            values <- do.call(match.fun(parameters[[i]]$fun), 
                              args = c(list(cur.seq), parameters[[i]]$args))
            
            # Rescale
            data.frame(var = nm, min = min(values), max = max(values))
        }) %>% bind_rows()
    } else {
        lapply(names(vars), function(nm) {
            values <- vars[[nm]]
            
            # Rescale
            data.frame(var = nm, min = min(values), max = max(values))
        }) %>% bind_rows()
    }
}

# Prediction function for virtual species for the case of being omniscient
pred_equation <- function(sp, newdata){
    # Reconstruct new data
    newdata <- newdata %>% select(sp$details$variables)
    parameters <- sp$details$parameters
    minmaxs <- sp$details$rescale.each.response.parameters
    
    newdata <- lapply(1:length(parameters), function(i) {
        nm <- names(parameters)[i]
        cur.seq <- newdata[[nm]]
        values <- do.call(match.fun(parameters[[i]]$fun), 
                          args = c(list(cur.seq), parameters[[i]]$args))
        
        # Rescale
        minmax <- minmaxs %>% filter(var == nm)
        data.frame(x =  (values - minmax$min) / (minmax$max - minmax$min)) %>% 
            rename("{nm}" := x)
    }) %>% bind_cols()
    
    
    # Apply formula, rescale, and convert to PA
    prediction <- rlang::eval_tidy(
        rlang::parse_expr(sp$details$formula), data = newdata)
    
    prediction <- (prediction - min(prediction)) / 
        (max(prediction) - min(prediction))
}

# Calculate the prevalence of a virtual species
prevalence <- function(x){
    y <- length(which(values(x) == 1)) / length(which(!is.na(values(x))))
    return(y)
}
