# Check species affected by NAs in bio14
# The truth is the analysis should be be affected that much by bio14

library(rnaturalearth)

fnames <- list.files("data/variables/OtherEnvMean", full.names = TRUE)

msk <- do.call(c, lapply(fnames, function(fn){
    rast(fn, lyrs = "bio14")
})) %>% mean()

msk[is.na(msk)] <- -999
msk[msk != -999] <- NA
msk[msk == -999] <- 1

msk <- as.polygons(msk)

afc <- ne_countries(scale = 50) %>% 
    filter(continent == "Africa") %>% 
    st_union() %>% st_transform(crs(msk))

msk <- st_as_sf(msk) %>% st_intersection(afc)

dispersal_rates <- read.csv("data/variables/species_dispersal_rate.csv")

fnames <- list.files(file.path("data/variables", "variable_list"),
                     full.names = TRUE)
fnames <- fnames[!grepl("_allruns", fnames)]
sp_check <- lapply(fnames, function(fn){
    if (gsub("_", " ", gsub(".csv", "", basename(fn))) %in% dispersal_rates$species){
        var_list <- read.csv(fn)
        if (length(na.omit(var_list$var_uncorrelated)) == 1){
            var_list <- na.omit(var_list$var)
        } else var_list <- na.omit(var_list$var_uncorrelated)
        
        if ("bio14" %in% var_list){
            dispersal_rate <- dispersal_rates %>% 
                filter(species == gsub("_", " ", gsub(".csv", "", basename(fn)))) %>% 
                pull(dispersal_rate)
            range <- st_read(
                file.path("data/IUCN/Expert_Maps",
                          sprintf("%s.geojson", gsub(".csv", "", basename(fn))))) %>% 
                st_transform(st_crs(msk)) %>% 
                st_buffer(dist = dispersal_rate * 90 * 1000)
            
            if (length(st_intersects(msk, range) %>% unlist()) == 0){
                affect_by_bio14 <- FALSE
                affect_area <- 0
            } else {
                affect_by_bio14 <- TRUE
                affect_area <- (sum(st_area(st_intersection(msk, range))) / 
                                    sum(st_area(range)) * 100) %>% 
                    units::drop_units() %>% round(0)
            }
        } else {
            affect_by_bio14 <- FALSE
            affect_area <- 0
        }
        
        data.frame(sp = gsub(".csv", "", basename(fn)),
                   affect_by_bio14 = affect_by_bio14, 
                   affect_area = affect_area)
    }
}) %>% bind_rows()
