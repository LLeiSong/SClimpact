library(dplyr)
root_dir <- "/home/lsong/SCImpact"
species_list <- read.csv(
    file.path(root_dir, "data/occurrences", "species_qualified_sdm.csv"))
species_list <- species_list$species

species_smr <- lapply(species_list, function(sp){
    dr <- file.path(root_dir, "results/sdm", sp)
    data.frame(species = sp,
               fnum = length(list.files(dr)))
}) %>% bind_rows()

if (file.exists(file.path(root_dir, "species_smr.csv"))){
    species_smr_old <- read.csv(file.path(root_dir, "species_smr.csv"))
    species_smr$change <- species_smr$fnum - species_smr_old$fnum
} else {
    species_smr$change <- 0
}

write.csv(species_smr, file.path(root_dir, "species_smr.csv"), row.names = FALSE)

species_done <- species_smr %>% filter(fnum == 18) %>% pull(species)
species_to_fsh <- species_smr %>% filter(fnum >= 15 & fnum < 18) %>% pull(species)
species_in_process <- species_smr %>% filter(fnum > 7 & fnum < 18) %>% pull(species)
species_has_process <- species_smr %>% filter(change > 0) %>% pull(species)
species_to_do <- species_smr %>% filter(fnum == 7) %>% pull(species)

message(sprintf("%s species are done.", length(species_done)))
message(sprintf("%s species are in process.", length(species_in_process)))
message(sprintf("- %s species have progress.", length(species_has_process)))
message(sprintf("- %s species are to be finished soon.", length(species_to_fsh)))
message(sprintf("%s species are waiting.", length(species_to_do)))
