library(stringr)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(rnaturalearth)
library(tidyterra)
sf_use_s2(FALSE)

# Get model evaluation
dst_dir <- "results/sdm"

evals <- list.files(dst_dir, full.names = TRUE, pattern = "evaluation_")
evals <- do.call(rbind, lapply(evals, function(fname){
    numbers <- read.csv(fname) %>% 
        mutate(species = str_extract(full.name, "[a-z|A-Z]+.[a-z|A-Z]+")) %>% 
        mutate(species = gsub("\\.", " ", species)) %>% 
        filter(PA != "allData" & run != "allRun") %>% 
        group_by(species, metric.eval) %>% 
        dplyr::summarise(worst = min(validation),
                         best = max(validation),
                         validation = mean(validation))
    
    if (numbers[3, 5] > 0.4){
        if (numbers[2, 5] > 0.7){
            numbers
        }
    }
}))

# 110 species left, these models can well explain the species distribution
# Visualize the results
evals <- evals %>% 
    mutate(metric.eval = factor(metric.eval,
                                levels = c("ACCURACY", "ROC", "TSS"),
                                labels = c("Accuracy", "AUC", "TSS")))
ggplot(data = evals, 
       aes(x = validation, after_stat(density), fill = metric.eval)) +
    geom_density(alpha = 0.30) +
    xlab("Metric value(0 - 1)") +
    ylab("Density") +
    scale_fill_aaas(name = "Metric") +
    theme_pubclean(base_size = 14)

ggsave("figures/model_eval.png",
       width = 4.5, height = 4.5, dpi = 500)

# Subset the species
species_list <- unique(evals$species)

# No. of associated species
template <- rast(file.path("data/env", "chelsa_1981-2010.tif"), lyrs = 1)
template <- aggregate(template, 6)
values(template) <- 0

fnames <- list.files(dst_dir, full.names = TRUE, pattern = "current_suit") %>% 
    data.frame(fname = .) %>% 
    mutate(species = gsub("current_suit_", "", basename(fname))) %>% 
    mutate(species = gsub(".tif", "", species)) %>% 
    mutate(species = gsub("_", " ", species)) %>% 
    filter(species %in% species_list)

for (fname in fnames$fname){
    lyr <- rast(fname)
    lyr <- !is.na(lyr)
    template <- terra::extend(lyr, template, fill = 0) + template
}

bry <- rnaturalearth::ne_countries(type = "countries")
template <- template %>% crop(bry) %>% mask(bry)
template[template == 0] <- NA
# template <- terra::trim(template)

ggplot() +
    geom_spatraster(data = template) +
    geom_sf(data = st_union(bry), color = "black", fill = "transparent") +
    scale_fill_whitebox_c(
        name = "No. of\nspecies",
        palette = "muted",
        n.breaks = 12,
        guide = guide_legend(reverse = TRUE))+
    theme_void(base_size = 14)

ggsave("figures/no_species.png",
       width = 8, height = 4, dpi = 500)

# Climate change analysis
template <- rast(file.path("data/env", "chelsa_1981-2010.tif"), lyrs = 1)
template <- aggregate(template, 6)
values(template) <- 0

climate_change <- function(lyr_nm, template){
    # Current
    fnames <- list.files(dst_dir, full.names = TRUE,
                         pattern = "current_shap") %>% 
        data.frame(fname = .) %>% 
        mutate(species = gsub("current_suit_", "", basename(fname))) %>% 
        mutate(species = gsub(".tif", "", species)) %>% 
        mutate(species = gsub("_", " ", species)) %>% 
        filter(species %in% species_list) %>% pull(fnames)
    
    ## Direction change
    dir_change <- do.call(c, lapply(fnames, function(fname){
        cur <- rast(fname) %>% subset(lyr_nm)
        fut <- rast(gsub("current", "future", fname)) %>% subset(lyr_nm)
        
        lyr <- ((cur >= 0) + 1) * ((fut >= 0) + 3)
        lyr %>% terra::extend(template, fill = NA)
    }))
    
    dir_change <- do.call(c, lapply(c(3, 4, 6, 8), function(x){
        lyrs <- do.call(c, lapply(dir_change, function(lyr){
            lyr == x
        }))
        
        mean(lyrs, na.rm = TRUE) * 100
    }))
    names(dir_change) <- c("N to N", "N to P", "P to N", "P to P")
    
    writeRaster(dir_change, file.path(sprintf("results/dir_change_%s.tif", lyr_nm)))
    message("Finish direction change detection.")
    
    # Value change
    val_change <- lapply(fnames, function(fname){
        cur <- rast(fname) %>% subset(lyr_nm)
        fut <- rast(gsub("current", "future", fname)) %>% subset(lyr_nm)
        pos_msk <- cur >= 0
        pos_msk[pos_msk == 0] <- NA
        neg_msk <- cur < 0
        neg_msk[neg_msk == 0] <- NA
        
        chg <- (fut - cur)
        lyrs <- c(mask(chg, pos_msk), mask(chg, neg_msk))
        names(lyrs) <- c("P", "N")
        
        lyrs %>% terra::extend(template, fill = NA)
    })
    
    val_change <- do.call(c, lapply(c("P", "N"), function(x){
        lyrs <- do.call(c, lapply(val_change, function(lyr){
            lyr[[x]]
        }))
        
        mean(lyrs, na.rm = TRUE)
    }))
    names(val_change) <- c("Suitable", "Non-suitable")
    
    writeRaster(val_change, file.path(sprintf("results/val_change_%s.tif", lyr_nm)))
    message("Finish value change detection.")
}

nms <- c("bio1", "bio12", "lc")
for (nm in nms){
    message(nm)
    climate_change(nm, template)
}

# Visualize
lyr_nms <- c("bio1", "bio12", "lc")
np <- do.call(c, lapply(lyr_nms, function(lyr_nm){
    rast(file.path(sprintf("results/dir_change_%s.tif", lyr_nm))) %>% 
        subset(2)
})); names(np) <- c("BIO1", "BIO12", "LC")

g1 <- ggplot() +
    geom_spatraster(data = np) +
    geom_sf(data = st_union(bry), color = "black", fill = "transparent") +
    facet_wrap(~lyr, nrow = 3, strip.position = "left") +
    scale_fill_whitebox_c(
        name = "(% of species)",
        palette = "muted",
        n.breaks = 10,
        guide = guide_legend(reverse = TRUE)
    ) + 
    theme_void() +
    theme(strip.text = element_text(size = 10),
          legend.key.size = unit(0.4, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(rep(0.1, 4)), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

pn <- do.call(c, lapply(lyr_nms, function(lyr_nm){
    rast(file.path(sprintf("results/dir_change_%s.tif", lyr_nm))) %>% 
        subset(3)
})); names(pn) <- c("BIO1", "BIO12", "LC")

g2 <- ggplot() +
    geom_spatraster(data = pn) +
    geom_sf(data = st_union(bry), color = "black", fill = "transparent") +
    facet_wrap(~lyr, nrow = 3, strip.position = "left") +
    scale_fill_whitebox_c(
        name = "(% of species)",
        palette = "muted",
        n.breaks = 10,
        guide = guide_legend(reverse = TRUE)
    ) + 
    theme_void() +
    theme(strip.text = element_text(size = 10, color = "transparent"),
          legend.key.size = unit(0.4, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(rep(0.1, 4)), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggarrange(g1, g2, ncol = 2, labels = c("A", "B"), 
          font.label = list(size = 10, face = "plain"),
          hjust = -20,
          common.legend = TRUE, legend = "right")

ggsave("figures/change_direction.png",
       width = 7, height = 4.5, dpi = 500)

# One by one
lyr_nm <- "bio12"
dir_change <- rast(file.path(sprintf("results/dir_change_%s.tif", lyr_nm)))
names(dir_change) <- c("A", "B", "C", "D")
val_change <- rast(file.path(sprintf("results/val_change_%s.tif", lyr_nm)))
names(val_change) <- c("A", "B")
ggplot() +
    geom_spatraster(data = dir_change) +
    geom_sf(data = st_union(bry), color = "black", fill = "transparent") +
    facet_wrap(~lyr, nrow = 2, strip.position = "left") +
    scale_fill_whitebox_c(
        name = "(% of species)",
        palette = "muted",
        n.breaks = 10,
        guide = guide_legend(reverse = TRUE)
    ) + 
    theme_void() +
    theme(strip.text = element_text(size = 12),
          legend.key.size = unit(0.4, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(rep(0.1, 4)), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave(sprintf("figures/%s_change_direction.png", lyr_nm),
       width = 8, height = 3.5, dpi = 500)

ggplot() +
    geom_spatraster(data = val_change) +
    geom_sf(data = st_union(bry), color = "black", fill = "transparent") +
    facet_wrap(~lyr, nrow = 2, strip.position = "left") +
    scale_fill_whitebox_c(
        name = "Shapley value\nchange",
        palette = "muted",
        n.breaks = 10, limits = c(-0.2, 0.25),
        guide = guide_legend(reverse = TRUE)
    ) + 
    theme_void() +
    theme(strip.text = element_text(size = 12),
          legend.key.size = unit(0.4, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(rep(0.1, 4)), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave(sprintf("figures/%s_change.png", lyr_nm),
       width = 6, height = 4.5, dpi = 500)
