library(stringr)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(rnaturalearth)
library(tidyterra)
library(forcats)
sf_use_s2(FALSE)

# Get model evaluation
root_dir <- here()
dst_dir <- "results/sdm"

species_list <- read.csv(
    file.path(here(), "data/occurrences", "species_qualified_sdm.csv")) %>% 
    arrange(num_occ)
species_list <- species_list$species

evals <- lapply(species_list, function(sp){
    dr <- file.path(root_dir, "results/sdm", sp)
    eval <- read.csv(file.path(dr, sprintf("cv_eval_%s.csv", sp)))
    eval <- eval %>% 
        mutate(accuracy = (tp + tn) / (tp + tn + fp + fn))
    eval <- eval %>% select(-fold) %>% group_by(type) %>% 
        summarise_all("mean") %>% 
        mutate(species = sp)
    if (all(eval$auc >= 0.7)){
        eval
    } else NULL
}) %>% bind_rows()

# 1992 species left, these models can well explain the species distribution
# Visualize the results
evals <- evals %>% 
    select(type, tss, auc, accuracy, species) %>% 
    pivot_longer(2:4, names_to = "metric.eval", values_to = "value")
evals <- evals %>% 
    mutate(metric.eval = factor(
        metric.eval,
        levels = c("accuracy", "auc", "tss"),
        labels = c("Accuracy", "AUC", "TSS")))
ggplot(data = evals %>% filter(type == "testing"), 
       aes(x = value, after_stat(density), fill = metric.eval)) +
    geom_density(alpha = 0.30) +
    xlab("Metric value(0 - 1)") +
    ylab("Density") +
    scale_fill_aaas(name = "Metric") +
    theme_pubclean(base_size = 14)

ggsave("docs/figures/model_eval.png",
       width = 4.5, height = 4.5, dpi = 500)

# Subset the species
species_list <- unique(evals$species)

# Variables
vars <- data.frame(
    var = names(rast(file.path("data/variables/Env", "AllEnv.tif"))),
    num = 0)
for (sp in species_list) {
    dr <- file.path(root_dir, "data/variables/variable_list")
    vars_sl <- read.csv(file.path(dr, sprintf("%s.csv", sp)))
    vars_sl <- na.omit(vars_sl$var_uncorrelated)
    vars[vars$var %in% vars_sl, 'num'] <- 
        vars[vars$var %in% vars_sl, 'num'] + 1
}

vars <- vars %>% arrange(-num) %>% slice(1:10) %>% 
    mutate(var = ifelse(var == "forest", "Forest coverage",
                        ifelse(var == "human_impact", "Human impact",
                               gsub("bio", "BIO", var)))) %>% 
    mutate(var = fct_reorder(var, -num))

ggplot(data = vars, 
       aes(x = var, num)) +
    geom_col(fill = "black") +
    xlab("Variable") +
    ylab("No. of species") +
    scale_fill_aaas(name = "Metric") +
    theme_pubclean(base_size = 14) +
    theme(axis.text.x = element_text(color = "black",
                                     angle = 90, vjust = 0.5, hjust=1))

ggsave("docs/figures/vars_selected.png",
       width = 4.5, height = 5, dpi = 500)

# Climate change analysis
template <- rast(file.path("data/variables/Env", "AllEnv.tif"), lyrs = 1)
values(template) <- 0

climate_change <- function(lyr_nm, scenario, template){
    # Current
    species_use <- lapply(species_list, function(sp){
        var_list <- read.csv(
            file.path("data/variables/variable_list",
                      sprintf("%s.csv", sp))) %>% 
            pull(var_uncorrelated) %>% na.omit()
        if (lyr_nm %in% var_list){
            sp
        } else NULL
    }) %>% unlist() %>% na.omit()
    
    fnames <- lapply(species_use, function(sp){
        list.files(file.path(dst_dir, sp), pattern = "shap_base", 
                   full.names = TRUE)
    }) %>% unlist()
    
    ## Direction change
    dir_change <- do.call(c, lapply(fnames, function(fname){
        cur <- rast(fname) %>% subset(lyr_nm)
        fut <- rast(gsub("base", scenario, fname)) %>% subset(lyr_nm)
        
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
        fut <- rast(gsub("base", scenario, fname)) %>% subset(lyr_nm)
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

nms <- c("bio1", "bio12", "forest")
scenario <- "ssp370_2041-2070"
for (nm in nms){
    message(nm)
    climate_change(nm, template)
}

# Visualize
lyr_nms <- c("bio1", "bio12", "forest")
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
lyr_nm <- "bio1"
plot_crs <- "+proj=peirce_q +lon_0=25 +shape=diamond"
dir_change <- rast(file.path(sprintf("results/dir_change_%s.tif", lyr_nm)))
val_change <- rast(file.path(sprintf("results/val_change_%s.tif", lyr_nm)))

dir_change <- project(subset(dir_change, 3), plot_crs)
ggplot() +
    geom_spatraster(data = dir_change, maxcell = 5e+06) +
    scale_fill_grass_c(
        name = "(% of species)",
        palette = "inferno", 
        n.breaks = 10,
        guide = guide_legend(reverse = TRUE)
    ) + 
    theme_void() +
    theme(strip.text = element_text(size = 12),
          legend.key.size = unit(0.4, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(rep(0.1, 4)), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          legend.position = "none")

ggsave(sprintf("docs/figures/%s_change_direction_pn.png", lyr_nm),
       width = 32, height = 32, dpi = 500)

val_change <- project(subset(val_change, 1), plot_crs)
ggplot() +
    geom_spatraster(data = val_change, maxcell = 5e+06) +
    scale_fill_whitebox_c(
        name = "(% of species)",
        palette = "muted", 
        n.breaks = 10, limits = c(-0.2, 0.25),
        guide = guide_legend(reverse = TRUE)
    ) + 
    theme_void() +
    theme(strip.text = element_text(size = 12),
          legend.key.size = unit(0.4, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(rep(0.1, 4)), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          legend.position = "none")

ggsave(sprintf("docs/figures/%s_change_s.png", lyr_nm),
       width = 32, height = 32, dpi = 500)

ggplot() +
    geom_spatraster(data = subset(dir_change, 2:3)) +
    scale_fill_grass_c(
        name = "(% of species)",
        palette = "inferno", 
        n.breaks = 10,
        guide = guide_legend(reverse = TRUE)
    ) + 
    facet_wrap(~lyr, nrow = 1) +
    theme_void() +
    theme(strip.text = element_text(size = 12),
          legend.key.size = unit(0.4, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(rep(0.1, 4)), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave(sprintf("docs/figures/%s_change_direction.png", lyr_nm),
       width = 9, height = 3, dpi = 500)

ggplot() +
    geom_spatraster(data = subset(dir_change, 2)) +
    scale_fill_grass_c(
        name = "(% of species)",
        palette = "inferno", 
        n.breaks = 10,
        guide = guide_legend(reverse = TRUE)
    ) + 
    theme_void() +
    theme(strip.text = element_text(size = 12),
          legend.key.size = unit(0.4, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(rep(0.1, 4)), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave(sprintf("docs/figures/%s_change_direction.png", lyr_nm),
       width = 9, height = 3, dpi = 500)

ggplot() +
    geom_spatraster(data = val_change) +
    scale_fill_whitebox_c(
        name = "Shapley value\nchange",
        palette = "muted",
        n.breaks = 10, limits = c(-0.2, 0.25),
        guide = guide_legend(reverse = TRUE)
    ) + 
    facet_wrap(~lyr, nrow = 1) +
    theme_void() +
    theme(strip.text = element_text(size = 12),
          legend.key.size = unit(0.4, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(rep(0.1, 4)), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave(sprintf("docs/figures/%s_change_s.png", lyr_nm),
       width = 9, height = 3, dpi = 500)

