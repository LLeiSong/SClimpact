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
}

# Function to sample occurrences using environmental bins
varela_sample_occ <- function (dat, no_bin) {
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
    out_ptz <- left_join(out_ptz, sub_ptz)
    
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
                           predictorsToKeepNoMatterWhat,
                           tieRule = 'last',
                           corThreshold = .7){
    out=try({
        if(!is.null(predictorsToKeepNoMatterWhat)){
            keepers <- which(names(env) %in% predictorsToKeepNoMatterWhat)
            envToKeep <- env[keepers]
            env <- env[-keepers]
        }
        c1 <- suppressWarnings(cor(env, use = 'complete.obs'))		
        tossed <- NULL
        bad <- which(apply(c1, 1, function(x){
            all(is.na(x)) | all(is.nan(x)) | all(is.null(x)) | sum(x,na.rm=T)==1
        }))
        
        if(length(bad)>0){
            tossed=c(tossed,names(bad)) 
            c1.index=which(colnames(c1) %in% tossed)
            c1=c1[-c1.index,-c1.index]
        }
        # number of variables a variable is too correlated with
        too.cor <- apply(abs(c1) > corThreshold,1,sum) - 1 # since diagonal always 1
        
        while(any(too.cor>=1)){
            most.cor=too.cor[too.cor==max(too.cor)]
            if(tieRule=='last') {
                toss=tail(most.cor, 1)
            } else if(tieRule=='random'){
                toss=sample(1:length(most.cor),1)
            }
            tossed=c(tossed,names(toss)) 
            c1.index=which(colnames(c1) %in% names(toss))
            c1=c1[-c1.index,-c1.index]
            too.cor=apply(abs(c1)>corThreshold,1,sum)-1
        }
        if(!is.null(tossed)){
            env=env[-which(names(env) %in% tossed)]
        } else {tossed='none' }
        if(!is.null(predictorsToKeepNoMatterWhat)){
            env <- cbind(envToKeep, env)
        }	
        env
    })
    names(out)
}
