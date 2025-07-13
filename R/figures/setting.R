#### Load libraries ####
library(stringr)
library(terra)
library(sf)
library(lwgeom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(magick)
library(ggsci)
library(ggpubr)
library(rnaturalearth)
library(rmapshaper)
library(tidyterra)
library(forcats)
library(tidytext)
library(biscale)
library(Hmisc)
library(spatstat)
library(here)
library(ggtext)
library(showtext)
library(flextable)
library(officer)
library(parallel)
library(colorspace)
library(cowplot)
sf_use_s2(FALSE)
mask <- terra::mask
get_legend <- ggpubr::get_legend

#### Model evaluation ####
root_dir <- here()
sdm_dir <- file.path(root_dir, "results/sdm")
sp_analysis_dir <- file.path(root_dir, "results/species_analysis")
cc_dir <- file.path(root_dir, "results/climate_change")
fig_dir <- file.path(root_dir, "results/figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir)
tbl_dir <- file.path(root_dir, "results/tables")
if (!dir.exists(tbl_dir)) dir.create(tbl_dir)

#### Global parameters ####
time_periods <- c("2011-2040", "2041-2070", "2071-2100")
ssps <- c("ssp126", "ssp370", "ssp585")
analysis_crs <- "EPSG:3857"
plot_crs <- "ESRI:54030"

#### Load qualified species list ####
species_list <- read.csv(
    file.path(root_dir, "results/species_reliable_final.csv")) %>% 
    pull(species)

#### Load selected environmental features ####
var_list <- lapply(species_list, function(sp){
    var_list <- read.csv(
        file.path(root_dir, "data/variables/variable_list",
                  sprintf("%s.csv", sp))) %>% 
        pull(var_uncorrelated) %>% na.omit()
}) %>% unlist()

var_list <- table(var_list) / length(species_list) * 100
var_list <- sort(var_list[var_list > 10], decreasing = TRUE)
features <- names(var_list)

#### Load IUCN status ####
sp_names <- read.csv(
    file.path(root_dir, "data/occurrences", 
              "species_catalog_0018860-240906103802322.csv"))
sp_names <- sp_names %>% select(species, Red_List_Category_2019) %>% 
    rename(category = Red_List_Category_2019) %>% 
    filter(species %in% gsub("_", " ", species_list)) %>% 
    mutate(category = factor(
        category, levels = c("CR", "EN", "VU", "NT", "LC", "DD"))) %>% 
    group_by(species) %>% arrange(category) %>% slice_head(n = 1) %>% 
    # Vulnerable, Endangered, or Critically Endangered is endangered
    mutate(status = ifelse(category %in% c("CR", "EN", "VU"), "EN", "NEN"))

#### Define layer for low, middle and high latitude areas ####
areas <- st_bbox(c(xmin = -180, ymin = -90, xmax = 180, ymax = 90)) %>% 
    st_as_sfc() %>% st_as_sf(crs = 4326)
mid_lines <- st_linestring(rbind(c(-180, 30), c(180, 30))) %>% 
    st_sfc(crs = 4326) %>% st_as_sf() %>% 
    rbind(st_linestring(rbind(c(-180, -30), c(180, -30))) %>% 
              st_sfc(crs = 4326) %>% st_as_sf())

high_lines <- st_linestring(rbind(c(-180, 60), c(180, 60))) %>% 
    st_sfc(crs = 4326) %>% st_as_sf() %>% 
    rbind(st_linestring(rbind(c(-180, -60), c(180, -60))) %>% 
              st_sfc(crs = 4326) %>% st_as_sf())

areas <- st_split(areas, st_union(mid_lines, high_lines)$x) %>% 
    st_collection_extract("POLYGON") %>% 
    mutate(area = c("High", "Middle", "Low", "Middle", "High")) %>% 
    rename(geometry = x) %>% 
    st_transform(crs = analysis_crs)

rm(mid_lines, high_lines)

#### World boundary without antarctic ####
bry <- ne_countries(scale = "medium") %>%
    filter(continent != "Antarctica") %>%
    st_union() %>% st_transform(plot_crs)

#### Variable groups ####
grps <- data.frame(var = paste0("BIO", c(1, 5, 9, 3, 4, 2, 8, 7)),
                   group = "Temperature", color = "#e66101") %>% 
    rbind(data.frame(var = paste0("BIO", c(14, 12, 15, 18, 19)),
                     group = "Precipitation", color = "#072ac8")) %>% 
    rbind(data.frame(var = c("FOR", "HLU", "GRA"),
                     group = "Land cover", color = "#1c2541")) %>% 
    mutate(group = factor(
        group, levels = c("Temperature", "Precipitation", "Land cover")))
