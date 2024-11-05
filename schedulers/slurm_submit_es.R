root_dir <- "/home/lsong/SCImpact/data"

species <- list.files(file.path(root_dir, "occurrences/CSVs"))
species <- gsub(".csv", "", species)
species2 <- list.files(file.path(root_dir, "variables", "variable_list"))
species2 <- gsub(".csv", "", species2)
species <- intersect(species, species2)

species1 <- list.files(file.path(root_dir, "occurrences/CSVs_thin"))
species1 <- gsub(".csv", "", species1)
species2 <- list.files(file.path(root_dir, "occurrences/bg"))
species2 <- gsub(".csv", "", species2)
species2 <- intersect(species1, species2)
species <- setdiff(species, species2)

message(sprintf("%s species to process.", length(species)))

# # If use one node
# sp <- paste(species, collapse = ",")
# system(sprintf("sbatch schedulers/env_sampling.sh %s", sp))

chunk_length <- 100
species_list <- split(
    species, ceiling(seq_along(species) / chunk_length))

for (i in 1:length(species_list)){
    sp <- species_list[[i]]
    sp <- paste(sp, collapse = ",")
    system(sprintf("sbatch schedulers/env_sampling.sh %s", sp))
}