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
library(here)
sf_use_s2(FALSE)
library(ggtext)
library(showtext)
library(flextable)
font_add_google('Merriweather')
showtext_auto()
showtext_opts(dpi = 500)
mask <- terra::mask

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
