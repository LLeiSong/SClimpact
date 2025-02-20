library(dplyr)
root_dir <- "/home/lsong/SCImpact"
species_list <- read.csv("/bigscratch/lsong/species_qualified_sdm.csv")
species_list <- species_list$species

var_list <- lapply(species_list, function(sp){
    var_list <- read.csv(
        file.path(root_dir, "data/variables/variable_list",
                  sprintf("%s.csv", sp))) %>% 
        pull(var_uncorrelated) %>% na.omit()
}) %>% unlist()

var_list <- table(var_list) / length(species_list) * 100
var_list <- sort(var_list[var_list > 10], decreasing = TRUE)
features <- names(var_list)

catalog <- lapply(features, function(feature){
    msg <- system(
        sprintf("sbatch schedulers/climate_change.sh %s", feature), 
        intern = TRUE)
    
    data.frame(feature = feature, slurm = msg)
}) %>% bind_rows()

write.csv(catalog, file.path(root_dir, "climate_change_slurm.csv"), row.names = FALSE)
