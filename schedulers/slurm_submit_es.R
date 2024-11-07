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

# When all species are done, run this:
fnames <- list.files(file.path(root_dir, "occurrences", "CSVs_thin"), 
                     full.names = TRUE)
## Get the numbers
nums <- lapply(fnames, function(fname){
    occ <- read.csv(fname)
    bg <- read.csv(gsub("CSVs_thin", "bg", fname))
    data.frame(species = gsub(".csv", "", basename(fname)),
               num_occ = nrow(occ),
               num_bg = nrow(bg))
}) %>% bind_rows()

## Only keep the species with >=20 occurrences
nums <- nums %>% filter(num_occ >= 20)
fname <- file.path(root_dir, "occurrences", "species_qualified_sdm.csv")
write.csv(nums, fname, row.names = FALSE)
