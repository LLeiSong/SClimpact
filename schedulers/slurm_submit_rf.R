root_dir <- "/home/lsong/SCImpact/data"
species_list <- read.csv(
    file.path(root_dir, "occurrences", "species_qualified_sdm.csv")) %>% 
    arrange(num_occ)
species_list <- species_list$species

message(sprintf("%s species to process.", length(species_list)))

chunk_length <- 40
species_list <- split(
    species_list, ceiling(seq_along(species_list) / chunk_length))

for (i in 1:length(species_list)){
    sp <- species_list[[i]]
    sp <- paste(sp, collapse = ",")
    system(sprintf("sbatch schedulers/rf_dws.sh %s", sp))
}