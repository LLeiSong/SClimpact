library(dplyr)
root_dir <- "/home/lsong/SCImpact"
species_list <- read.csv(
    file.path(root_dir, "data/occurrences", "species_qualified_sdm.csv")) %>% 
    arrange(num_occ)
species_list <- species_list$species

species_done <- lapply(species_list, function(sp){
    dr <- file.path(root_dir, "results/sdm", sp)
    if (dir.exists(dr)){
        if (length(list.files(dr)) == 54){
            sp
        } else NULL
    } else NULL
})

species_done <- species_done[!sapply(species_done, is.null)]
species_done <- unlist(species_done)

species_list <- setdiff(species_list, species_done)

message(sprintf("%s species to process.", length(species_list)))

chunk_length <- 20
species_list <- split(
    species_list, ceiling(seq_along(species_list) / chunk_length))

for (i in 1:length(species_list)){
    sp <- species_list[[i]]
    sp <- paste(sp, collapse = ",")
    system(sprintf("sbatch schedulers/shap.sh %s", sp))
}