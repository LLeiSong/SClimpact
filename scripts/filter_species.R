# Utilities: filter the species
# Author: Lei Song (lsong@ucsb.edu)

# Get the list of species
library(here)
library(sf)
library(dplyr)
library(readxl)
library(stringr)

# Define directory
iucn_dir <- here("data/IUCN")

# Taxon groups
taxons <- c("MAMMALS", "REPTILES", "AMPHIBIANS", "BIRDS")

for (taxon in taxons) {
    message(sprintf("Process taxon: %s", taxon))
    if (taxon != "BIRDS"){
        fnames <- list.files(file.path(iucn_dir, taxon), pattern = ".shp$",
                             full.names = TRUE)
        plys <- do.call(rbind, lapply(fnames, st_read)) %>% 
            filter(marine == "false" & terrestria == "true") %>% 
            # Temporarily select these
            filter(category %in% c("EN", "CR")) %>% 
            filter(presence %in% c(1, 2, 3)) %>% 
            filter(origin %in% c(1, 2))
    } else {
        plys <- st_read(file.path(iucn_dir, "BOTW.gdb"), 
                        layer = "All_Species") %>% 
            filter(presence %in% c(1, 2, 3)) %>% 
            filter(origin %in% c(1, 2))
        
        # Read categories
        fname <- file.path(iucn_dir, "HBW-BirdLife_Checklist_v7_Dec22",
                           "HBW_BirdLife_List of birds_v7.xlsx")
        sheet <- excel_sheets(fname)
        categories <- read_excel(fname, sheet = sheet) %>% 
            filter(`2022 IUCN Red List category` %in% c("EN", "CR")) %>% 
            rename(sci_name = `Scientific name`,
                   category = `2022 IUCN Red List category`) %>% 
            select(sci_name, category)
        
        plys <- plys %>% filter(sci_name %in% categories$sci_name)
        plys <- left_join(plys, categories)
    }
    
    # Save out
    st_write(plys, file.path(iucn_dir, sprintf("%s_endangered.geojson", taxon)))
}

# Make a catalog
fnames <- list.files(iucn_dir, pattern = "endangered.geojson", full.names = TRUE)

catalog <- do.call(rbind, lapply(fnames, function(fname){
    if (str_detect(fname, "BIRDS")){
        Sys.setenv(OGR_GEOJSON_MAX_OBJ_SIZE = 0)
        st_read(fname) %>% st_drop_geometry() %>% 
            mutate(class = "BIRDS") %>% 
            select(sci_name, presence, origin, seasonal, class, category) %>% 
            unique()
    } else {
        st_read(fname) %>% st_drop_geometry() %>% 
            select(sci_name, presence, origin, seasonal, class, category) %>% 
            unique()
    }
}))
write.csv(catalog, file.path(iucn_dir, "catalog_species.csv"), row.names = FALSE)
