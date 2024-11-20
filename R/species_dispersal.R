library(sf)
library(dplyr)
library(readxl)
library(rnaturalearth)

# Load species list and catalog
species_list <- read.csv(
    file.path("data/occurrences", "species_qualified_sdm.csv")) %>% 
    mutate(species = gsub("_", " ", species)) %>% select(species)

species_catalog <- read.csv(
    file.path("data/occurrences", 
              "species_catalog_0018860-240906103802322.csv"))

# Load species trait
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.gd0m3
traits_mammals <- read_excel(
    file.path("data/generation_length", 
              "Generation\ Lenght\ for\ Mammals.xlsx"), sheet = 1)

# Calculate the natal dispersal distance
## https://www.pnas.org/doi/epdf/10.1073/pnas.1116791109
traits_mammals <- traits_mammals %>% 
    mutate(AdultBodyMass_kg = AdultBodyMass_g / 1000) %>% 
    select(Order, Scientific_name, AdultBodyMass_kg, GenerationLength_d) %>% 
    mutate(dispersal_distance = ifelse(Order == "Carnivora",
                                       3.45 * AdultBodyMass_kg^0.89,
                                       1.45 * AdultBodyMass_kg^0.54)) %>% 
    mutate(dispersal_rate = dispersal_distance / (GenerationLength_d / 365))

# Attach the values to species
species_dispersal <- left_join(
    species_list, traits_mammals, 
    by = c('species' = 'Scientific_name'))

# Deal with missing species for dispersal distance
## Using the median of the order
order_dispersal <- species_dispersal %>% 
    filter(!is.na(Order)) %>% 
    group_by(Order) %>% 
    summarise(dispersal_rate = median(dispersal_rate, na.rm = TRUE))

# One missing Order, using the most similar Order instead
order_dispersal <- order_dispersal %>% 
    rbind(data.frame(Order = "Artiodactyla",
        dispersal_rate = order_dispersal[
            order_dispersal$Order == "Perissodactyla", "dispersal_rate"]))

# Concatenate everything together
species_left <- species_dispersal %>% filter(is.na(Order))
species_catalog <- species_catalog %>% 
    filter(species %in% species_left$species) %>% 
    select(species, order) %>% unique()
species_left <- left_join(species_left, species_catalog, by = "species") %>% 
    select(species, order) %>% 
    left_join(order_dispersal, by = c("order" = "Order"))

species_dispersal <- species_dispersal %>% 
    filter(!is.na(Order)) %>% 
    select(species, Order, dispersal_rate) %>% 
    rename(order = Order) %>% 
    rbind(species_left)

write.csv(species_dispersal, 
          file.path("data/variables/species_dispersal_rate.csv"),
          row.names = FALSE)

# Because considering dispersal will affect landmass crossing, 
# so save a continent file as well
continents <- ne_countries(scale = 50) %>% 
    filter(continent != "Antarctica") %>% 
    st_transform(3857) %>% 
    summarise(geometry = st_union(geometry)) %>% 
    st_cast("POLYGON") %>% 
    mutate(id = 1:nrow(.)) %>% 
    select(id, geometry)
st_write(continents, file.path("data/variables/continents.geojson")) 
