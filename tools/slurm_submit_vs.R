root_dir <- "/home/lsong/SCImpact/data"
species_list <- list.files(file.path(root_dir, "occurrences/CSVs"))
species_list <- gsub(".csv", "", species_list)

fnames <- list.files(file.path(root_dir, "variables/variable_list"))
fnames_to_move <- list.files(file.path(root_dir, "variables/variable_list"),
                             pattern = "allruns.csv")
fnames <- setdiff(fnames, fnames_to_move)
fnames <- gsub(".csv", "", fnames)
species_list <- setdiff(species_list, fnames)

message(sprintf("%s species to process.", length(species_list)))

# chunk_length <- 2
# species_list <- split(
# species_list, ceiling(seq_along(species_list) / chunk_length))

for (i in 1:length(species_list)){
    sp <- species_list[[i]]
    sp <- paste(sp, collapse = ",")
    print(sp)
    system(sprintf("sbatch schedulers/select_variable.sh %s", sp))
}