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
library(tidytext)
library(biscale)
library(Hmisc)
library(here)
sf_use_s2(FALSE)
library(ggtext)
library(showtext)
font_add_google('Merriweather')
showtext_auto()
showtext_opts(dpi = 500)
mask <- terra::mask

#### Model evaluation ####
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

write.csv(evals, "results/model_evaluation.csv", row.names = FALSE)

ggplot(data = evals %>% filter(type == "testing"), 
       aes(x = value, after_stat(density), fill = metric.eval)) +
    geom_density(alpha = 0.8) +
    xlab("Metric value(0 - 1)") +
    ylab("Density") +
    scale_fill_brewer(name = "Metric", palette = "Dark2") +
    theme_pubclean(base_size = 12, base_family = "Merriweather") +
    theme(axis.text = element_text(color = "black", family = 'Merriweather'))

ggsave("docs/figures/model_eval.png",
       width = 4, height = 3, dpi = 500)

#### Variable selection ####
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

vars <- vars %>% mutate(ratio = num / length(species_list)) %>% 
    filter(ratio > 0.1) %>% 
    mutate(var = ifelse(var == "forest", "Forest coverage",
                        ifelse(var == "human_impact", "Human impact",
                               gsub("bio", "BIO", var)))) %>% 
    mutate(var = ifelse(var == "grassland", "Grassland", var)) %>% 
    mutate(var = fct_reorder(var, num))

write.csv(vars, "results/variable_selection.csv", row.names = FALSE)

ggplot(data = vars, 
       aes(x = num, xend = 0,
           y = var, yend = var)) +
    geom_segment() +
    geom_point() +
    xlab("No. of species") + ylab("") +
    theme_pubclean(base_size = 12, base_family = "Merriweather") +
    theme(axis.text.x = element_text(
        color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(vjust = -1))

ggsave("docs/figures/vars_selected.png",
       width = 4, height = 3.5, dpi = 500)

# Before removing correlated variables
vars <- data.frame(
    var = names(rast(file.path("data/variables/Env", "AllEnv.tif"))),
    num = 0)
iter <- 30
for (sp in species_list) {
    dr <- file.path(root_dir, "data/variables/variable_list")
    vars_sl <- read.csv(file.path(dr, sprintf("%s_allruns.csv", sp)))
    
    ## Rearrange variables based on the voting results
    vars_selected <- vars_sl %>% group_by(var) %>% 
        summarise(n = n()) %>% 
        filter(n >= floor(iter * 0.5)) %>% # 50% agree
        arrange(-n) %>% pull(var)
    
    if (length(vars_selected) == 0){
        vars_selected <- vars_sl %>% group_by(var) %>% 
            summarise(n = n()) %>% 
            filter(n == max(n)) %>% 
            pull(var)
    }
    
    vars[vars$var %in% vars_selected, 'num'] <- 
        vars[vars$var %in% vars_selected, 'num'] + 1
}

vars <- vars %>% mutate(ratio = num / length(species_list)) %>% 
    filter(ratio > 0.1) %>% 
    mutate(var = ifelse(var == "forest", "Forest coverage",
                        ifelse(var == "human_impact", "Human impact",
                               gsub("bio", "BIO", var)))) %>% 
    mutate(var = ifelse(var == "grassland", "Grassland", var)) %>% 
    mutate(var = fct_reorder(var, num))

write.csv(vars, "results/variable_selection_raw.csv", row.names = FALSE)

ggplot(data = vars, 
       aes(x = num, xend = 0,
           y = var, yend = var)) +
    geom_segment() +
    geom_point() +
    xlab("No. of species") + ylab("") +
    theme_pubclean(base_size = 12, base_family = "Merriweather") +
    theme(axis.text.x = element_text(
        color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(vjust = -1))

ggsave("docs/figures/vars_selected_raw.png",
       width = 4, height = 4.5, dpi = 500)

#### Impact of drivers ####
##### Affected area ####
data_dir <- "results/climate_change"
time_periods <- c("2011-2040", "2041-2070", "2071-2100")

# All drivers
var_list <- lapply(species_list, function(sp){
    var_list <- read.csv(
        file.path(root_dir, "data/variables/variable_list",
                  sprintf("%s.csv", sp))) %>% 
        pull(var_uncorrelated) %>% na.omit()
}) %>% unlist()

var_list <- table(var_list) / length(species_list) * 100
var_list <- sort(var_list[var_list > 10], decreasing = TRUE)
features <- names(var_list)

features_periods <- lapply(time_periods, function(time_period){
    do.call(rbind, lapply(features, function(driver){
        fnames <- list.files(
            data_dir, pattern = sprintf("%s.tif", time_period), full.names = TRUE)
        fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
        
        dir_fnames <- fnames[str_detect(fnames, "dir_change")]
        num_fnames <- fnames[str_detect(fnames, "num_species")]
        
        # S to U
        stous <- do.call(c, lapply(1:length(dir_fnames), function(i){
            rast(dir_fnames[i])[[3]] / rast(num_fnames[i])[[1]]
        })); names(stous) <- c("SSP126", "SSP370", "SSP580")
        stous <- trim(stous)
        
        # U to S
        utoss <- do.call(c, lapply(1:length(dir_fnames), function(i){
            rast(dir_fnames[i])[[2]] / rast(num_fnames[i])[[1]]
        })); names(utoss) <- c("SSP126", "SSP370", "SSP580")
        utoss <- trim(utoss)
        
        all_area <- global(
            mask(cellSize(stous, unit = "km"), stous), 
            "sum", na.rm = TRUE)[[1]]
        
        stous[stous < 0.01] <- NA
        utoss[utoss < 0.01] <- NA
        stous_area <- global(
            mask(cellSize(stous, unit = "km"), stous),
            "sum", na.rm = TRUE)[[1]]
        utoss_area <- global(
            mask(cellSize(stous, unit = "km"), utoss),
            "sum", na.rm = TRUE)[[1]]
        
        data.frame(driver = driver, 
                   stou_percent = stous_area / all_area * 100,
                   utos_percent = utoss_area / all_area * 100) %>% 
            mutate(scenario = c("SSP126", "SSP370", "SSP585"))
    })) %>% mutate(time_period = time_period)
}) %>% bind_rows()

write.csv(features_periods, "results/affect_area.csv", row.names = FALSE)

# Manipulate data for plotting
pts <- features_periods %>% 
    mutate(driver = ifelse(
        driver == "forest", "Forest coverage",
        ifelse(driver == "human_impact", "Human impact",
               gsub("bio", "BIO", driver)))) %>% 
    mutate(driver = ifelse(driver == "grassland", "Grassland", driver))

drivers <- pts %>% filter(time_period == "2011-2040") %>% 
    group_by(driver) %>% summarise(stou_percent = mean(stou_percent)) %>% 
    arrange(stou_percent) %>% pull(driver)

pts <- pts %>% 
    pivot_longer(2:3, names_to = "type", values_to = "value") %>% 
    mutate(type = factor(
        type, levels = c("utos_percent", "stou_percent"),
        labels = c("From unsuitable to suitable", 
                   "From suitable to unsuitable"))) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))
segments <- pts %>% 
    group_by(driver, type, time_period) %>% 
    summarise(across(where(is.numeric), max)) %>% 
    pivot_wider(names_from = time_period, values_from = value) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

# Plot
ggplot() +
    geom_linerange(
        data = segments,
        aes(x = driver, ymin = `2011-2040`, ymax = `2041-2070`, group = type),
        col = 'grey60', position = position_dodge(.8)) +
    geom_linerange(
        data = segments,
        aes(x = driver, ymin = `2041-2070`, ymax = `2071-2100`, group = type),
        col = 'grey60', position = position_dodge(.8)) +
    coord_flip() +
    geom_point(data = pts,
               aes(y = value, x = driver, col = time_period, group = type),
               color = "white", fill = "white", size = 2, 
               position = position_dodge(.8)) +
    geom_point(data = pts,
               aes(y = value, x = driver, col = time_period, 
                   group = type, shape = type, alpha = scenario), size = 2, 
               position = position_dodge(.8)) +
    labs(y = "Percentage of affected area",
         x = element_blank()) +
    scale_color_brewer(name = "Time period", palette = "Dark2") +
    scale_shape_manual(name = "", values = c(1, 16)) +
    scale_alpha_manual(name = "Scenario", values = c(0.4, 0.7, 1.0)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) +
    theme_pubclean(base_size = 11, base_family = 'Merriweather') +
    theme(axis.text = element_text(
        color = "black", size = 11),
        strip.background = element_blank(),
        legend.text = element_text(size = 11),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.justification = c(1.1, 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        axis.title.x = element_text(vjust = -1),
        legend.box = 'vertical')

ggsave("docs/figures/affected_area.png",
       width = 5.5, height = 7, dpi = 500)

##### Affected area simple ####
# Manipulate data for plotting
features_periods <- features_periods %>% 
    group_by(driver, time_period) %>% 
    summarise(utos_percent = mean(utos_percent),
              stou_percent = mean(stou_percent))
pts <- features_periods %>% 
    mutate(driver = ifelse(
        driver == "forest", "Forest coverage",
        ifelse(driver == "human_impact", "Human impact",
               gsub("bio", "BIO", driver)))) %>% 
    mutate(driver = ifelse(driver == "grassland", "Grassland coverage", driver))

# Get the driver rank for visualization
drivers <- pts %>% filter(time_period == "2011-2040") %>% 
    arrange(stou_percent) %>% pull(driver)

# Manipulate the data for plotting
pts <- pts %>% 
    pivot_longer(3:4, names_to = "type", values_to = "value") %>% 
    mutate(type = factor(
        type, levels = c("utos_percent", "stou_percent"),
        labels = c("From unsuitable to suitable", 
                   "From suitable to unsuitable"))) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

# For transition between time stamps
segments <- pts %>% 
    pivot_wider(names_from = time_period, values_from = value) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

# Plot
ggplot() +
    geom_linerange(
        data = segments,
        aes(x = driver, ymin = `2011-2040`, ymax = `2041-2070`, group = type),
        col = 'grey60', position = position_dodge(.8)) +
    geom_linerange(
        data = segments,
        aes(x = driver, ymin = `2041-2070`, ymax = `2071-2100`, group = type),
        col = 'grey60', position = position_dodge(.8)) +
    coord_flip() +
    geom_point(data = pts,
               aes(y = value, x = driver, group = type),
               color = "white", fill = "white", size = 2.2, 
               position = position_dodge(.8)) +
    geom_point(data = pts,
               aes(y = value, x = driver, color = time_period,
                   group = type, shape = type), size = 2, 
               position = position_dodge(.8)) +
    geom_point(data = pts %>% filter(time_period != "2011-2040"),
               aes(y = value - 1.2, x = driver, group = type), color = "grey60",
               shape = ">", size = 3, 
               position = position_dodge(.8)) +
    labs(y = "Percentage of affected area",
         x = element_blank()) +
    scale_color_brewer(name = "Time period", palette = "Dark2") +
    scale_shape_manual(name = "", values = c(1, 16)) +
    scale_alpha_manual(name = "Scenario", values = c(0.4, 0.7, 1.0)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) +
    theme_pubclean(base_size = 11, base_family = 'Merriweather') +
    theme(axis.text = element_text(
        color = "black", size = 11),
        strip.background = element_blank(),
        legend.text = element_text(size = 11),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.justification = c(1.3, 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        axis.title.x = element_text(vjust = -1),
        legend.box = 'vertical')

ggsave("docs/figures/affected_area_simple.png",
       width = 5.5, height = 7, dpi = 500)

##### Affected area polar ####
pts <- features_periods %>% 
    mutate(driver = ifelse(
        driver == "forest", "Forest\ncoverage",
        ifelse(driver == "human_impact", "Human\nimpact",
               gsub("bio", "BIO", driver)))) %>% 
    mutate(driver = ifelse(driver == "grassland", "Grassland\ncoverage", driver))

# Get the driver rank for visualization
drivers <- pts %>% 
    group_by(driver) %>% 
    summarise(utos_percent = max(utos_percent),
              stou_percent = max(stou_percent)) %>% 
    mutate(stou_rank = rank(-stou_percent), 
           utos_rank = rank(-utos_percent)) %>% 
    mutate(SqRank = (stou_rank^2) + (utos_rank^2)/2) %>% 
    mutate(RankOrder = rank(SqRank)) %>% 
    arrange(RankOrder)

pts <- pts %>% 
    mutate(driver = factor(
        driver, levels = drivers$driver, labels = drivers$driver))

ref_values <- pts %>% filter(time_period == "2011-2040")
labels2 <- pts %>% filter(time_period == "2041-2070")
labels2$utos_percent <- (labels2$utos_percent - ref_values$utos_percent) / 2 + 
    ref_values$utos_percent
labels2$stou_percent <- (labels2$stou_percent - ref_values$stou_percent) / 2 + 
    ref_values$stou_percent

legend_data <- data.frame(y = c(rep(1, 3), rep(2, 3)), 
                          x = rep(1:3, 2)) %>% 
    mutate(y = factor(
        y, levels = 1:2, labels = c("From unsuitable to suitable",
                                    "From suitable to unsuitable")),
        x = factor(x, levels = 1:3,
                   labels = c("2011-2040", "2041-2070", "2071-2100"))) %>% 
    mutate(color = y, alpha = x)

lgd <- ggplot() +
    geom_tile(data = legend_data, 
              aes(x = x, y = y, fill = factor(color),
                  alpha = alpha), color = "white", linewidth = 2) +
    scale_fill_manual(values = c("#018571", '#a6611a')) +
    scale_alpha_manual(values = c(1, 0.6, 0.4)) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(
              family = 'Merriweather', size = 8, color = "black"))

g <- ggplot() + 
    geom_col(data = pts %>% filter(time_period == "2071-2100"), 
             color = "white", fill = "#a6611a", alpha = 0.4,
             aes(x = driver, y = stou_percent)) +
    geom_col(data = pts %>% filter(time_period == "2041-2070"), 
             color = "white", fill = "#a6611a", alpha = 0.6,
             aes(x = driver, y = stou_percent)) +
    geom_col(data = pts %>% filter(time_period == "2011-2040"), 
             color = "white", fill = "#a6611a",
             aes(x = driver, y = stou_percent)) +
    geom_col(data = pts %>% filter(time_period == "2071-2100"), 
             color = "white", fill = "#018571", alpha = 0.4,
             aes(x = driver, y = -utos_percent)) +
    geom_col(data = pts %>% filter(time_period == "2041-2070"), 
             color = "white", fill = "#018571", alpha = 0.6,
             aes(x = driver, y = -utos_percent)) +
    geom_col(data = pts %>% filter(time_period == "2011-2040"),
             aes(x = driver, y = -utos_percent), 
             color = "white", fill = "#018571") +
    # text for suitable to unsuitable
    geom_text(data = pts %>% filter(time_period == "2071-2100"), 
              aes(x = driver, y = stou_percent + 5, 
                  label = round(stou_percent, 0)), 
              color = '#a6611a', size = 2, family = "Merriweather") +
    geom_text(data = labels2, 
              aes(x = driver, y = stou_percent, 
                  label = round(stou_percent, 0)), 
              color = 'black', size = 2, family = "Merriweather") +
    geom_text(data = pts %>% filter(time_period == "2011-2040"), 
              aes(x = driver, y = stou_percent / 2, 
                  label = round(stou_percent, 0)), 
              color = 'white', size = 2, family = "Merriweather") +
    # text for unsuitable to suitable
    geom_text(data = pts %>% filter(time_period == "2071-2100"), 
              aes(x = driver, y = -utos_percent - 5, 
                  label = round(utos_percent, 0)), 
              color = '#018571', size = 2, family = "Merriweather") +
    geom_text(data = labels2, 
              aes(x = driver, y = -utos_percent, 
                  label = round(utos_percent, 0)), 
              color = 'black', size = 2, family = "Merriweather") +
    geom_text(data = pts %>% filter(time_period == "2011-2040"), 
              aes(x = driver, y = -utos_percent / 2, 
                  label = round(utos_percent, 0)), 
              color = 'white', size = 2, family = "Merriweather") +
    labs(x = "", y = "") +
    scale_y_continuous(limits = c(-100, 90)) +
    coord_polar() + theme_minimal(base_family = "Merriweather") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(color = "black", size = 8),
          plot.margin = unit(c(-0.5, -1, -1, -1), "cm"))

ggarrange(lgd, g, nrow = 2, heights = c(0.2, 1))

ggsave("docs/figures/affected_area_polar.png",
       width = 5.5, height = 5.5, dpi = 500, bg = "white")

##### Affected species ####
features_periods <- lapply(time_periods, function(time_period){
    do.call(rbind, lapply(features, function(driver){
        fnames <- list.files(
            data_dir, pattern = sprintf("%s.tif", time_period), full.names = TRUE)
        fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
        
        dir_fnames <- fnames[str_detect(fnames, "dir_change")]
        num_fnames <- fnames[str_detect(fnames, "num_species")]
        
        # S to U
        stous <- do.call(c, lapply(1:length(dir_fnames), function(i){
            rast(dir_fnames[i])[[3]] / rast(num_fnames[i])[[1]]
        })); names(stous) <- c("SSP126", "SSP370", "SSP580")
        stous <- trim(stous)
        
        # U to S
        utoss <- do.call(c, lapply(1:length(dir_fnames), function(i){
            rast(dir_fnames[i])[[2]] / rast(num_fnames[i])[[1]]
        })); names(utoss) <- c("SSP126", "SSP370", "SSP580")
        utoss <- trim(utoss)
        
        stous_sp <- c(cellSize(stous, unit = "km"), stous)
        stous_sp <- values(stous_sp) %>% na.omit() %>% as.data.frame()
        
        ssp126 <- stous_sp %>% select(area, SSP126) %>% filter(SSP126 > 0.01)
        ssp370 <- stous_sp %>% select(area, SSP370) %>% filter(SSP370 > 0.01)
        ssp580 <- stous_sp %>% select(area, SSP580) %>% filter(SSP580 > 0.01)
        
        stous_sds <- c(wtd.var(ssp126$SSP126, ssp126$area),
                       wtd.var(ssp370$SSP370, ssp370$area),
                       wtd.var(ssp580$SSP580, ssp580$area))
        
        stous_means <- c(wtd.mean(ssp126$SSP126, ssp126$area),
                         wtd.mean(ssp370$SSP370, ssp370$area),
                         wtd.mean(ssp580$SSP580, ssp580$area))
        
        utoss_sp <- c(cellSize(utoss, unit = "km"), utoss)
        utoss_sp <- values(utoss_sp) %>% na.omit() %>% as.data.frame()
        
        ssp126 <- utoss_sp %>% select(area, SSP126) %>% filter(SSP126 > 0.01)
        ssp370 <- utoss_sp %>% select(area, SSP370) %>% filter(SSP370 > 0.01)
        ssp580 <- utoss_sp %>% select(area, SSP580) %>% filter(SSP580 > 0.01)
        
        utoss_sds <- c(wtd.var(ssp126$SSP126, ssp126$area),
                       wtd.var(ssp370$SSP370, ssp370$area),
                       wtd.var(ssp580$SSP580, ssp580$area))
        
        utoss_means <- c(wtd.mean(ssp126$SSP126, ssp126$area),
                         wtd.mean(ssp370$SSP370, ssp370$area),
                         wtd.mean(ssp580$SSP580, ssp580$area))
        
        data.frame(driver = driver, 
                   stou_sp_mean = stous_means * 100,
                   stou_sp_sd = stous_sds * 100,
                   utos_sp_mean = utoss_means * 100,
                   utos_sp_sd = utoss_sds * 100) %>% 
            mutate(scenario = c("SSP126", "SSP370", "SSP585"))
    })) %>% mutate(time_period = time_period)
}) %>% bind_rows()

write.csv(features_periods, "results/affect_species.csv", row.names = FALSE)

# Manipulate data for plotting
pts <- features_periods %>% 
    select(driver, stou_sp_mean, utos_sp_mean, scenario, time_period) %>% 
    mutate(driver = ifelse(
        driver == "forest", "Forest coverage",
        ifelse(driver == "human_impact", "Human impact",
               gsub("bio", "BIO", driver)))) %>% 
    mutate(driver = ifelse(driver == "grassland", "Grassland", driver))

drivers <- pts %>% filter(time_period == "2011-2040") %>% 
    group_by(driver) %>% summarise(stou_sp_mean = min(stou_sp_mean)) %>% 
    arrange(stou_sp_mean) %>% pull(driver)

pts <- pts %>% 
    pivot_longer(2:3, names_to = "type", values_to = "value") %>% 
    mutate(type = factor(
        type, levels = c("utos_sp_mean", "stou_sp_mean"),
        labels = c("From unsuitable to suitable", 
                   "From suitable to unsuitable"))) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

segments <- pts %>% 
    group_by(driver, type) %>% 
    summarise(vmin = min(value), vmax = max(value)) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

# Plot
ggplot() +
    geom_linerange(
        data = segments,
        aes(x = driver, ymin = vmin, ymax = vmax, group = type),
        col = 'grey60', position = position_dodge(.8)) +
    coord_flip() +
    geom_point(data = pts,
               aes(y = value, x = driver, col = time_period, group = type),
               color = "white", fill = "white", size = 2, 
               position = position_dodge(.8)) +
    geom_point(data = pts,
               aes(y = value, x = driver, col = time_period, 
                   group = type, shape = type, alpha = scenario), size = 2, 
               position = position_dodge(.8)) +
    labs(y = "Percentage of affected species per pixel",
         x = element_blank()) +
    scale_color_brewer(name = "Time period", palette = "Dark2") +
    scale_shape_manual(name = "", values = c(1, 16)) +
    scale_alpha_manual(name = "Scenario", values = c(0.4, 0.7, 1.0)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) +
    theme_pubclean(base_size = 11, base_family = 'Merriweather') +
    theme(axis.text = element_text(
        color = "black", size = 11),
        strip.background = element_blank(),
        legend.text = element_text(size = 11),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.justification = c(1.3, 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        axis.title.x = element_text(vjust = -1),
        legend.box = 'vertical')

ggsave("docs/figures/affected_species.png",
       width = 5.5, height = 7, dpi = 500)

# The table
tb <- features_periods
tb$stou_sp <- sprintf("%.1f±%.1f", tb$stou_sp_mean, tb$stou_sp_sd)
tb$utos_sp <- sprintf("%.1f±%.1f", tb$utos_sp_mean, tb$utos_sp_sd)
tb <- tb %>% select(driver, scenario, time_period, stou_sp, utos_sp)

##### Affected species simple ####
# Manipulate data for plotting
features_periods <- features_periods %>% 
    group_by(driver, time_period) %>% 
    summarise(stou_sp_mean = mean(stou_sp_mean),
              utos_sp_mean = mean(utos_sp_mean))

pts <- features_periods %>% 
    mutate(driver = ifelse(
        driver == "forest", "Forest coverage",
        ifelse(driver == "human_impact", "Human impact",
               gsub("bio", "BIO", driver)))) %>% 
    mutate(driver = ifelse(driver == "grassland", "Grassland coverage", driver))

drivers <- pts %>% filter(time_period == "2011-2040") %>% 
    arrange(stou_sp_mean) %>% pull(driver)

pts <- pts %>% 
    pivot_longer(3:4, names_to = "type", values_to = "value") %>% 
    mutate(type = factor(
        type, levels = c("utos_sp_mean", "stou_sp_mean"),
        labels = c("From unsuitable to suitable", 
                   "From suitable to unsuitable"))) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

segments <- pts %>% 
    pivot_wider(names_from = time_period, values_from = value) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

# arrows
pts_arrows_pos <- pts %>% 
    pivot_wider(names_from = time_period, values_from = value) %>% 
    mutate(range1 = `2041-2070` - `2011-2040`,
           range2 = `2071-2100` - `2041-2070`) %>% 
    mutate(`2041-2070` = ifelse(range1 < 1, NA, `2041-2070`),
           `2071-2100` = ifelse(range2 < 1, NA, `2071-2100`))

pts_arrows_neg <- pts %>% 
    pivot_wider(names_from = time_period, values_from = value) %>% 
    mutate(range1 = `2041-2070` - `2011-2040`,
           range2 = `2071-2100` - `2041-2070`) %>% 
    mutate(`2041-2070` = ifelse(range1 <= -1, `2041-2070`, NA),
           `2071-2100` = ifelse(range2 <= -1, `2071-2100`, NA))

# Plot
ggplot() +
    geom_linerange(
        data = segments,
        aes(x = driver, ymin = `2011-2040`, ymax = `2041-2070`, group = type),
        col = 'grey60', position = position_dodge(.8)) +
    geom_linerange(
        data = segments,
        aes(x = driver, ymin = `2041-2070`, ymax = `2071-2100`, group = type),
        col = 'grey60', position = position_dodge(.8)) +
    coord_flip() +
    geom_point(data = pts,
               aes(y = value, x = driver, group = type),
               color = "white", fill = "white", size = 2.2,
               position = position_dodge(.8)) +
    geom_point(data = pts,
               aes(y = value, x = driver, color = time_period,
                   group = type, shape = type), size = 2, 
               position = position_dodge(.8)) +
    geom_point(data = pts_arrows_pos,
               aes(y = `2041-2070` - 0.5, x = driver, group = type), 
               color = "grey60", shape = ">", size = 3,
               position = position_dodge(.8)) +
    geom_point(data = pts_arrows_pos,
               aes(y = `2071-2100` - 0.5, x = driver, group = type), 
               color = "grey60", shape = ">", size = 3,
               position = position_dodge(.8)) +
    geom_point(data = pts_arrows_neg,
               aes(y = `2041-2070` + 0.5, x = driver, group = type), 
               color = "grey60", shape = "<", size = 3,
               position = position_dodge(.8)) +
    geom_point(data = pts_arrows_neg,
               aes(y = `2071-2100` + 0.5, x = driver, group = type), 
               color = "grey60", shape = "<", size = 3,
               position = position_dodge(.8)) +
    labs(y = "Percentage of affected species per pixel",
         x = element_blank()) +
    scale_color_brewer(name = "Time period", palette = "Dark2") +
    scale_shape_manual(name = "", values = c(1, 16)) +
    scale_alpha_manual(name = "Scenario", values = c(0.4, 0.7, 1.0)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) +
    theme_pubclean(base_size = 11, base_family = 'Merriweather') +
    theme(axis.text = element_text(
        color = "black", size = 11),
        strip.background = element_blank(),
        legend.text = element_text(size = 11),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.justification = c(1.3, 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        axis.title.x = element_text(vjust = -1),
        legend.box = 'vertical')

ggsave("docs/figures/affected_species_simple.png",
       width = 5.5, height = 7, dpi = 500)

##### Affected species polar ####
pts <- features_periods %>% 
    mutate(driver = ifelse(
        driver == "forest", "Forest\ncoverage",
        ifelse(driver == "human_impact", "Human\nimpact",
               gsub("bio", "BIO", driver)))) %>% 
    mutate(driver = ifelse(driver == "grassland", "Grassland\ncoverage", driver))

# Get the driver rank for visualization
drivers <- pts %>% 
    group_by(driver) %>% 
    summarise(stou_sp_mean = max(stou_sp_mean),
              utos_sp_mean = max(utos_sp_mean)) %>% 
    mutate(stou_rank = rank(-stou_sp_mean), 
           utos_rank = rank(-utos_sp_mean)) %>% 
    mutate(SqRank = (stou_rank^2) + (utos_rank^2)/2) %>% 
    mutate(RankOrder = rank(SqRank)) %>% 
    arrange(RankOrder)

pts <- pts %>% 
    mutate(driver = factor(
        driver, levels = drivers$driver, labels = drivers$driver))

bio14 <- pts %>% filter(driver == "BIO14")
pts[pts$driver == "BIO14", "utos_sp_mean"] <- NA

ref_values <- pts %>% filter(time_period == "2011-2040")
labels2 <- pts %>% filter(time_period == "2041-2070")
labels2$utos_sp_mean <- (labels2$utos_sp_mean - ref_values$utos_sp_mean) / 2 + 
    ref_values$utos_sp_mean
labels2$stou_sp_mean <- (labels2$stou_sp_mean - ref_values$stou_sp_mean) / 2 + 
    ref_values$stou_sp_mean

legend_data <- data.frame(y = c(rep(1, 3), rep(2, 3)), 
                          x = rep(1:3, 2)) %>% 
    mutate(y = factor(
        y, levels = 1:2, labels = c("From unsuitable to suitable",
                                    "From suitable to unsuitable")),
        x = factor(x, levels = 1:3,
                   labels = c("2011-2040", "2041-2070", "2071-2100"))) %>% 
    mutate(color = y, alpha = x)

lgd <- ggplot() +
    geom_tile(data = legend_data, 
              aes(x = x, y = y, fill = factor(color),
                  alpha = alpha), color = "white", linewidth = 2) +
    scale_fill_manual(values = c("#018571", '#a6611a')) +
    scale_alpha_manual(values = c(1, 0.6, 0.4)) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(
              family = 'Merriweather', size = 8, color = "black"))

g <- ggplot() + 
    geom_col(data = pts %>% filter(time_period == "2071-2100"), 
             color = "white", fill = "#a6611a", alpha = 0.4,
             aes(x = driver, y = stou_sp_mean)) +
    geom_col(data = pts %>% filter(time_period == "2041-2070"), 
             color = "white", fill = "#a6611a", alpha = 0.6,
             aes(x = driver, y = stou_sp_mean)) +
    geom_col(data = pts %>% filter(time_period == "2011-2040"), 
             color = "white", fill = "#a6611a",
             aes(x = driver, y = stou_sp_mean)) +
    geom_col(data = pts %>% filter(time_period == "2071-2100"), 
             color = "white", fill = "#018571", alpha = 0.4,
             aes(x = driver, y = -utos_sp_mean)) +
    geom_col(data = pts %>% filter(time_period == "2041-2070"), 
             color = "white", fill = "#018571", alpha = 0.6,
             aes(x = driver, y = -utos_sp_mean)) +
    geom_col(data = pts %>% filter(time_period == "2011-2040"),
             aes(x = driver, y = -utos_sp_mean), 
             color = "white", fill = "#018571") +
    # text for suitable to unsuitable
    geom_text(data = pts %>% filter(time_period == "2071-2100"), 
              aes(x = driver, y = stou_sp_mean + 3, 
                  label = round(stou_sp_mean, 0)), 
              color = '#a6611a', size = 2, family = "Merriweather") +
    geom_text(data = labels2, 
              aes(x = driver, y = stou_sp_mean, 
                  label = round(stou_sp_mean, 0)), 
              color = 'black', size = 2, family = "Merriweather") +
    geom_text(data = pts %>% filter(time_period == "2011-2040"), 
              aes(x = driver, y = stou_sp_mean / 2, 
                  label = round(stou_sp_mean, 0)), 
              color = 'white', size = 2, family = "Merriweather") +
    # text for unsuitable to suitable
    geom_text(data = pts %>% filter(time_period == "2071-2100"), 
              aes(x = driver, y = -utos_sp_mean - 3, 
                  label = round(utos_sp_mean, 0)), 
              color = '#018571', size = 2, family = "Merriweather") +
    geom_text(data = labels2, 
              aes(x = driver, y = -utos_sp_mean, 
                  label = round(utos_sp_mean, 0)), 
              color = 'black', size = 2, family = "Merriweather") +
    geom_text(data = pts %>% filter(time_period == "2011-2040"), 
              aes(x = driver, y = -utos_sp_mean / 2, 
                  label = round(utos_sp_mean, 0)), 
              color = 'white', size = 2, family = "Merriweather") +
    geom_col(data = bio14 %>% filter(time_period == "2011-2040"),
             aes(x = driver, y = -utos_sp_mean), 
             color = "white", fill = "#018571") +
    geom_col(data = bio14 %>% filter(time_period == "2041-2070"),
             aes(x = driver, y = -utos_sp_mean), 
             color = "white", fill = "white") +
    geom_col(data = bio14 %>% filter(time_period == "2041-2070"),
             aes(x = driver, y = -utos_sp_mean), 
             color = "white", fill = "#018571", alpha = 0.6) +
    geom_col(data = bio14 %>% filter(time_period == "2071-2100"),
             aes(x = driver, y = -utos_sp_mean), 
             color = "white", fill = "#018571", alpha = 0.4) +
    geom_text(data = bio14 %>% filter(time_period == "2041-2070"), 
              aes(x = driver, y = -utos_sp_mean / 2, 
                  label = round(utos_sp_mean, 0)), 
              color = 'white', size = 2, family = "Merriweather") +
    geom_text(data = bio14 %>% filter(time_period == "2011-2040"), 
              aes(x = driver, y = -utos_sp_mean - 3, 
                  label = round(utos_sp_mean, 0)), 
              color = '#018571', size = 2, family = "Merriweather") +
    geom_text(data = bio14 %>% filter(time_period == "2071-2100"), 
              aes(x = driver, y = -utos_sp_mean, 
                  label = round(utos_sp_mean, 0)), 
              color = 'black', size = 2, family = "Merriweather") +
    labs(x = "", y = "") +
    scale_y_continuous(limits = c(-40, 40)) +
    coord_polar() + theme_minimal(base_family = "Merriweather") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(color = "black", size = 8),
          plot.margin = unit(c(-0.5, -1, -1, -1), "cm"))

ggarrange(lgd, g, nrow = 2, heights = c(0.2, 1))

ggsave("docs/figures/affected_species_polar.png",
       width = 5.5, height = 5.5, dpi = 500, bg = "white")

#### SHAP value change ####
features_periods <- lapply(time_periods, function(time_period){
    do.call(rbind, lapply(features, function(driver){
        fnames <- list.files(
            data_dir, pattern = sprintf("%s.tif", time_period), full.names = TRUE)
        fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
        
        val_fnames <- fnames[str_detect(fnames, sprintf("val_change_%s", driver))]
        num_fnames <- fnames[str_detect(fnames, "num_species")]
        
        # Suitable
        ss <- do.call(c, lapply(1:length(val_fnames), function(i){
            rast(val_fnames[i])[[1]] / rast(num_fnames[i])[[2]]
        })); names(ss) <- c("SSP126", "SSP370", "SSP580")
        ss <- trim(ss)
        
        # Unsuitable
        uss <- do.call(c, lapply(1:length(val_fnames), function(i){
            rast(val_fnames[i])[[2]] / rast(num_fnames[i])[[3]]
        })); names(uss) <- c("SSP126", "SSP370", "SSP580")
        uss <- trim(uss)
        
        ss_val <- c(cellSize(ss, unit = "km"), ss)
        ss_val <- values(ss_val) %>% na.omit() %>% as.data.frame()
        
        ssp126 <- ss_val %>% select(area, SSP126) %>%
            filter(SSP126 >= quantile(SSP126, 0.01) & SSP126 <= quantile(SSP126, 0.99))
        ssp370 <- ss_val %>% select(area, SSP370) %>%
            filter(SSP370 >= quantile(SSP370, 0.01) & SSP370 <= quantile(SSP370, 0.99))
        ssp580 <- ss_val %>% select(area, SSP580) %>%
            filter(SSP580 >= quantile(SSP580, 0.01) & SSP580 <= quantile(SSP580, 0.99))
        
        ss_sds <- c(wtd.var(ssp126$SSP126, ssp126$area),
                    wtd.var(ssp370$SSP370, ssp370$area),
                    wtd.var(ssp580$SSP580, ssp580$area))
        
        ss_means <- c(wtd.mean(ssp126$SSP126, ssp126$area),
                      wtd.mean(ssp370$SSP370, ssp370$area),
                      wtd.mean(ssp580$SSP580, ssp580$area))
        
        uss_val <- c(cellSize(uss, unit = "km"), uss)
        uss_val <- values(uss_val) %>% na.omit() %>% as.data.frame()
        
        ssp126 <- uss_val %>% select(area, SSP126) %>%
            filter(SSP126 >= quantile(SSP126, 0.01) & SSP126 <= quantile(SSP126, 0.99))
        ssp370 <- uss_val %>% select(area, SSP370) %>%
            filter(SSP370 >= quantile(SSP370, 0.01) & SSP370 <= quantile(SSP370, 0.99))
        ssp580 <- uss_val %>% select(area, SSP580) %>%
            filter(SSP580 >= quantile(SSP580, 0.01) & SSP580 <= quantile(SSP580, 0.99))
        
        uss_sds <- c(wtd.var(ssp126$SSP126, ssp126$area),
                     wtd.var(ssp370$SSP370, ssp370$area),
                     wtd.var(ssp580$SSP580, ssp580$area))
        
        uss_means <- c(wtd.mean(ssp126$SSP126, ssp126$area),
                       wtd.mean(ssp370$SSP370, ssp370$area),
                       wtd.mean(ssp580$SSP580, ssp580$area))
        
        data.frame(driver = driver, 
                   suitable_change = ss_means,
                   suitable_sd = ss_sds,
                   unsuitable_change = uss_means,
                   unsuitable_sd = uss_sds) %>% 
            mutate(scenario = c("SSP126", "SSP370", "SSP585"))
    })) %>% mutate(time_period = time_period)
}) %>% bind_rows()

write.csv(features_periods, "results/shap_value_change.csv", row.names = FALSE)

# Manipulate data for plotting
pts <- features_periods %>% 
    select(driver, suitable_change, unsuitable_change, scenario, time_period) %>% 
    mutate(driver = ifelse(driver == "forest", "Forest coverage",
                           ifelse(driver == "human_impact", "Human impact",
                                  gsub("bio", "BIO", driver)))) %>% 
    mutate(driver = ifelse(driver == "grassland", "Grassland", driver))

drivers <- pts %>% filter(time_period == "2011-2040") %>% 
    group_by(driver) %>% summarise(suitable_change = mean(suitable_change)) %>% 
    arrange(-suitable_change) %>% pull(driver)

pts <- pts %>% 
    pivot_longer(2:3, names_to = "type", values_to = "value") %>% 
    mutate(type = factor(
        type, levels = c("suitable_change", "unsuitable_change"),
        labels = c("Current suitable area", "Current unsuitable area"))) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

segments <- pts %>% 
    group_by(driver, type) %>% 
    summarise(vmin = min(value), vmax = max(value)) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

# Plot
ggplot() +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_linerange(
        data = segments,
        aes(x = driver, ymin = vmin, ymax = vmax, group = type),
        col = 'grey60', position = position_dodge(.8)) +
    coord_flip() +
    geom_point(data = pts,
               aes(y = value, x = driver, col = time_period, group = type),
               color = "white", fill = "white", size = 2, 
               position = position_dodge(.8)) +
    geom_point(data = pts,
               aes(y = value, x = driver, col = time_period, 
                   group = type, shape = type, alpha = scenario), size = 2, 
               position = position_dodge(.8)) +
    labs(y = "SHAP value change",
         x = element_blank()) +
    scale_color_brewer(name = "Time period", palette = "Dark2") +
    scale_shape_manual(name = "", values = c(16, 1)) +
    scale_alpha_manual(name = "Scenario", values = c(0.4, 0.7, 1.0)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) +
    theme_pubclean(base_size = 11, base_family = 'Merriweather') +
    theme(axis.text = element_text(
        color = "black", size = 11),
        strip.background = element_blank(),
        legend.text = element_text(size = 11),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.justification = c(-1, 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        axis.title.x = element_text(vjust = -1),
        legend.box = 'vertical')

ggsave("docs/figures/shap_change.png",
       width = 6.5, height = 6, dpi = 500)

# The table
tb <- features_periods
tb$suitable_change <- sprintf("%.3f±%.3f", tb$suitable_change, tb$suitable_sd)
tb$unsuitable_change <- sprintf("%.3f±%.3f", tb$unsuitable_change, tb$unsuitable_sd)
tb <- tb %>% select(driver, scenario, time_period, suitable_change, unsuitable_change)

##### Simple version ####
# Manipulate data for plotting
features_periods <- features_periods %>% 
    group_by(driver, time_period) %>% 
    summarise(suitable_change = mean(suitable_change),
              unsuitable_change = mean(unsuitable_change))

pts <- features_periods %>% 
    mutate(driver = ifelse(driver == "forest", "Forest coverage",
                           ifelse(driver == "human_impact", "Human impact",
                                  gsub("bio", "BIO", driver)))) %>% 
    mutate(driver = ifelse(driver == "grassland", "Grassland coverage", driver))

drivers <- pts %>% filter(time_period == "2011-2040") %>% 
    arrange(-suitable_change) %>% pull(driver)

pts <- pts %>% 
    pivot_longer(3:4, names_to = "type", values_to = "value") %>% 
    mutate(type = factor(
        type, levels = c("suitable_change", "unsuitable_change"),
        labels = c("Current suitable area", "Current unsuitable area"))) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

segments <- pts %>% 
    pivot_wider(names_from = time_period, values_from = value) %>% 
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

# arrows
pts_arrows_pos <- pts %>% 
    pivot_wider(names_from = time_period, values_from = value) %>% 
    mutate(range1 = `2041-2070` - `2011-2040`,
           range2 = `2071-2100` - `2041-2070`) %>% 
    mutate(`2041-2070` = ifelse(range1 < 0.003, NA, `2041-2070`),
           `2071-2100` = ifelse(range2 < 0.003, NA, `2071-2100`))

pts_arrows_neg <- pts %>% 
    pivot_wider(names_from = time_period, values_from = value) %>% 
    mutate(range1 = `2041-2070` - `2011-2040`,
           range2 = `2071-2100` - `2041-2070`) %>% 
    mutate(`2041-2070` = ifelse(range1 <= -0.003, `2041-2070`, NA),
           `2071-2100` = ifelse(range2 <= -0.003, `2071-2100`, NA))

# Plot
ggplot() +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_linerange(
        data = segments,
        aes(x = driver, ymin = `2011-2040`, ymax = `2041-2070`, group = type),
        col = 'grey60') +
    geom_linerange(
        data = segments,
        aes(x = driver, ymin = `2041-2070`, ymax = `2071-2100`, group = type),
        col = 'grey60') +
    coord_flip() +
    geom_point(data = pts,
               aes(y = value, x = driver, col = time_period, group = type),
               color = "white", fill = "white", size = 2.2) +
    geom_point(data = pts,
               aes(y = value, x = driver, col = time_period, 
                   group = type, shape = type), size = 2) +
    geom_point(data = pts_arrows_pos,
               aes(y = `2041-2070` - 0.0015, x = driver, group = type), 
               color = "grey60", shape = ">", size = 3) +
    geom_point(data = pts_arrows_pos,
               aes(y = `2071-2100` - 0.0015, x = driver, group = type), 
               color = "grey60", shape = ">", size = 3) +
    geom_point(data = pts_arrows_neg,
               aes(y = `2041-2070` + 0.0015, x = driver, group = type), 
               color = "grey60", shape = "<", size = 3) +
    geom_point(data = pts_arrows_neg,
               aes(y = `2071-2100` + 0.0015, x = driver, group = type), 
               color = "grey60", shape = "<", size = 3) +
    labs(y = "SHAP value change",
         x = element_blank()) +
    scale_color_brewer(name = "Time period", palette = "Dark2") +
    scale_shape_manual(name = "", values = c(16, 1)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) +
    theme_pubclean(base_size = 11, base_family = 'Merriweather') +
    theme(axis.text = element_text(
        color = "black", size = 11),
        legend.text = element_text(size = 11),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.justification = c(-3, 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        axis.title.x = element_text(vjust = -1, hjust = 0.32),
        plot.subtitle = element_text(size = 10, vjust = -25, hjust = 2.7),
        legend.box = 'vertical')

ggsave("docs/figures/shap_change_simple.png",
       width = 6.5, height = 6, dpi = 500)

#### Individual impact ####
# Set interested driver and projection for visualization
plot_crs <- "ESRI:54030"

# World boundary without antarctic
bry <- ne_countries(scale = "large") %>%
    filter(continent != "Antarctica") %>%
    st_union() %>% st_transform(plot_crs)

for (driver in features){
    # Figure of temperature from suitable to unsuitable
    values_periods <- lapply(time_periods, function(time_period){
        fnames <- list.files(
            data_dir, pattern = sprintf("%s.tif", time_period), 
            full.names = TRUE)
        fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
        
        dir_fnames <- fnames[str_detect(fnames, "dir_change")]
        num_fnames <- fnames[str_detect(fnames, "num_species")]
        
        lyrs <- do.call(c, lapply(1:length(dir_fnames), function(i){
            project(rast(dir_fnames[i])[[3]] / rast(num_fnames[i])[[1]] * 100, 
                    plot_crs)
        })); names(lyrs) <- c("SSP126", "SSP370", "SSP580")
        lyrs <- trim(lyrs)
        
        means <- mean(lyrs, na.rm = TRUE)
        sds <- stdev(lyrs, na.rm = TRUE, pop = TRUE)
        
        # Put values into data.frame
        data <- as.data.frame(c(means, sds), xy = TRUE)
        names(data) <- c("x", "y", "Mean", "SD")
        
        data %>% mutate(time_period = time_period)
    })
    values_periods <- values_periods %>% bind_rows()
    
    ## Convert values to biviariates
    values_bivi <- bi_class(
        values_periods, x = Mean, y = SD, 
        style = "fisher", dim = 4, dig_lab = 4)
    breaks <- bi_class_breaks(
        values_periods, x = Mean, y = SD, style = "fisher", 
        dim = 4, dig_lab = 4, split = TRUE)
    breaks$bi_x <- round(breaks$bi_x, 0)
    breaks$bi_y <- round(breaks$bi_y, 0)
    
    lgd <- ggplotGrob(
        bi_legend(pal = "PurpleOr",
                  pad_color = "white",
                  flip_axes = FALSE,
                  rotate_pal = FALSE,
                  dim = 4,
                  xlab = "Mean (%)",
                  ylab = "SD (%)",
                  breaks = breaks,
                  pad_width = 0.5,
                  size = 10) +
            theme(axis.text = element_text(
                size = 8, color = "black", family = 'Merriweather'),
                axis.title = element_text(
                    size = 10, color = "black", family = 'Merriweather'),
                panel.background = element_rect(fill='transparent'),
                plot.background = element_rect(fill='transparent', color = NA),
                axis.title.x = element_text(vjust = 3.5),
                axis.title.y = element_text(vjust = -3)))
    
    gs <- lapply(time_periods, function(tp){
        values <- values_bivi %>% filter(time_period == tp)
        
        # Make the figure
        ggplot() +
            geom_sf(data = bry, fill = "#d3d3d3", color = "#d3d3d3") + 
            geom_tile(
                data = values, 
                mapping = aes(x = x, y = y, fill = bi_class), 
                show.legend = FALSE) +
            coord_sf(crs = plot_crs) +
            bi_scale_fill(
                pal = "PurpleOr", dim = 4, 
                flip_axes = FALSE, rotate_pal = FALSE) +
            annotation_custom(grob = lgd, 
                              xmin = -20208181, xmax = -7808181,
                              ymin = -10342702, ymax = 2057298) +
            theme_void() + 
            theme(plot.margin = unit(c(-0.2, -15, -0.2, -15), "cm"))
    })
    
    ggarrange(plotlist = gs, nrow = 3, 
              labels = c("A", "B", "C"), hjust = -2, vjust = 3,
              font.label = list(
                  size = 12, face = "bold", family = 'Merriweather'))
    
    ggsave(sprintf("docs/figures/%s_sus_mean_sd.png", driver),
           width = 6, height = 8, dpi = 500)
    
    # Figure of temperature from unsuitable to suitable
    values_periods <- lapply(time_periods, function(time_period){
        fnames <- list.files(
            data_dir, pattern = sprintf("%s.tif", time_period), 
            full.names = TRUE)
        fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
        
        dir_fnames <- fnames[str_detect(fnames, "dir_change")]
        num_fnames <- fnames[str_detect(fnames, "num_species")]
        
        lyrs <- do.call(c, lapply(1:length(dir_fnames), function(i){
            project(rast(dir_fnames[i])[[2]] / rast(num_fnames[i])[[1]] * 100, 
                    plot_crs)
        })); names(lyrs) <- c("SSP126", "SSP370", "SSP580")
        lyrs <- trim(lyrs)
        
        means <- mean(lyrs, na.rm = TRUE)
        sds <- stdev(lyrs, na.rm = TRUE, pop = TRUE)
        
        # Put values into data.frame
        data <- as.data.frame(c(means, sds), xy = TRUE)
        names(data) <- c("x", "y", "Mean", "SD")
        
        data %>% mutate(time_period = time_period)
    })
    values_periods <- values_periods %>% bind_rows()
    
    ## Convert values to biviariates
    values_bivi <- bi_class(
        values_periods, x = Mean, y = SD, 
        style = "fisher", dim = 4, dig_lab = 4)
    breaks <- bi_class_breaks(
        values_periods, x = Mean, y = SD, style = "fisher", 
        dim = 4, dig_lab = 4, split = TRUE)
    breaks$bi_x <- round(breaks$bi_x, 0)
    breaks$bi_y <- round(breaks$bi_y, 0)
    
    lgd <- ggplotGrob(
        bi_legend(pal = "PurpleOr",
                  pad_color = "white",
                  flip_axes = FALSE,
                  rotate_pal = FALSE,
                  dim = 4,
                  xlab = "Mean (%)",
                  ylab = "SD (%)",
                  breaks = breaks,
                  pad_width = 0.5,
                  size = 10) +
            theme(axis.text = element_text(
                size = 8, color = "black", family = 'Merriweather'),
                axis.title = element_text(
                    size = 10, color = "black", family = 'Merriweather'),
                panel.background = element_rect(fill='transparent'),
                plot.background = element_rect(fill='transparent', color = NA),
                axis.title.x = element_text(vjust = 3.5),
                axis.title.y = element_text(vjust = -3)))
    
    gs <- lapply(time_periods, function(tp){
        values <- values_bivi %>% filter(time_period == tp)
        
        # Make the figure
        ggplot() +
            geom_sf(data = bry, fill = "#d3d3d3", color = "#d3d3d3") + 
            geom_tile(
                data = values, 
                mapping = aes(x = x, y = y, fill = bi_class), 
                show.legend = FALSE) +
            coord_sf(crs = plot_crs) +
            bi_scale_fill(
                pal = "PurpleOr", dim = 4, 
                flip_axes = FALSE, rotate_pal = FALSE) +
            annotation_custom(grob = lgd, 
                              xmin = -20208181, xmax = -7808181,
                              ymin = -10342702, ymax = 2057298) +
            theme_void() + 
            theme(plot.margin = unit(c(-0.2, -15, -0.2, -15), "cm"))
    })
    
    ggarrange(plotlist = gs, nrow = 3, 
              labels = c("A", "B", "C"), hjust = -2, vjust = 3,
              font.label = list(
                  size = 12, face = "bold", family = 'Merriweather'))
    
    ggsave(sprintf("docs/figures/%s_uss_mean_sd.png", driver),
           width = 6, height = 8, dpi = 500)
}

#### Variable change - SHAP change ####
## Set parameters
driver <- "bio5"
s <- 1
time_period <- "2011-2040"
nm <- ifelse(driver == "forest", "Forest coverage",
             ifelse(driver == "human_impact", "Human impact",
                    gsub("bio", "BIO", driver)))
nm <- ifelse(nm == "grassland", "Grassland", nm)

sus <- ifelse(s == 1, "s", "us")

# Suitable area
shap_changes <- do.call(c, lapply(time_periods, function(time_period){
    fnames <- list.files(
        data_dir, pattern = sprintf("%s.tif", time_period), full.names = TRUE)
    fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
    
    val_fnames <- fnames[str_detect(fnames, sprintf("val_change_%s", driver))]
    num_fnames <- fnames[str_detect(fnames, "num_species")]
    
    ss <- do.call(c, lapply(1:length(val_fnames), function(i){
        rast(val_fnames[i])[[s]] / rast(num_fnames[i])[[s + 1]]
    })); names(ss) <- c("SSP126", "SSP370", "SSP580")
    
    mean(ss, na.rm = TRUE)
})); names(shap_changes) <- time_periods

val_current <- rast("data/variables/Env/AllEnv.tif", lyrs = driver)
val_changes <- do.call(c, lapply(time_periods, function(time_period){
    temp_future <- c(
        rast(sprintf("data/variables/OtherEnvMean/ssp126_%s.tif", 
                     time_period), lyrs = driver),
        rast(sprintf("data/variables/OtherEnvMean/ssp370_%s.tif", 
                     time_period), lyrs = driver),
        rast(sprintf("data/variables/OtherEnvMean/ssp585_%s.tif", 
                     time_period), lyrs = driver))
    (mean(temp_future, na.rm = TRUE) - val_current)
})); names(val_changes) <- time_periods

## Convert values to biviariates
## Calculate statistics
data <- c(cellSize(shap_changes[[time_period]], unit = "km"), 
          val_changes[[time_period]], shap_changes[[time_period]])
data <- as.data.frame(data, xy = TRUE) %>% na.omit()
names(data) <- c("x", "y", "area", "variable", "SHAP")

data$SHAP_br <- ifelse(data$SHAP <= 0, 1, 2)
data$SHAP_br <- factor(data$SHAP_br, levels = c(1, 2))
data$variable_br <- ifelse(data$variable <= 0, 1, 2)
data$variable_br <- factor(data$variable_br, levels = c(1, 2))
data <- bi_class(data, x = variable_br, y = SHAP_br, dim = 2)
total_area <- sum(data$area)
statistics <- data %>% 
    # filter(SHAP >= quantile(SHAP, 0.01) & 
    #            SHAP <= quantile(SHAP, 0.99)) %>% 
    dplyr::group_by(bi_class) %>% 
    summarise(variable_sd = wtd.var(variable, area),
              SHAP_sd = wtd.var(SHAP, area),
              SHAP = wtd.mean(SHAP, area),
              variable = wtd.mean(variable, area),
              area = sum(area) / total_area * 100)
statistics$area <- sprintf("%.2f%%", statistics$area)

digit_shap <- ceiling(log10(1 / mean(abs(statistics$SHAP))))
digit_var <- ceiling(log10(1 / mean(abs(statistics$variable))))

statistics$SHAP <- sprintf(
    "%.2f\u00B1%.2f", statistics$SHAP * 10^digit_shap, 
    statistics$SHAP_sd * 10^digit_shap)
statistics$variable <- sprintf(
    "%.2f\u00B1%.2f", statistics$variable * 10^digit_var, 
    statistics$variable_sd * 10^digit_var)
statistics <- statistics %>% select(bi_class, variable, SHAP, area)

color_table <- data.frame(
    category = c("1-1", "2-1", "1-2", "2-2"),
    color = c("#d3d3d3", "#52b6b6", "#ad5b9c", "#434e87"))

## Build table into a nice ggtextable() to visualize it
statistics <- as.data.frame(statistics)

colors <- sapply(statistics$bi_class, function(x){
    color_table %>% filter(category == x) %>% pull(color)
})

names(statistics) <- c(
    "Category~(Var-SHAP)", 
    sprintf("%s~change~(10^-%s)", nm, digit_var), 
    sprintf("SHAP~change~(10^-%s)", digit_shap), "Area (%)")

tbl <- ggtexttable(
    statistics[, -1], rows = NULL, 
    theme = ttheme(
        padding = unit(c(1, 2), "mm"),
        colnames.style = colnames_style(
            size = 8, fill = "transparent", parse = TRUE, 
            face = "plain", font = "Merriweather"),
        tbody.style = tbody_style(
            size = 8, color = "black", font = "Merriweather", 
            fill = colors)))

## Visualization
data <- c(cellSize(shap_changes[[time_period]], unit = "km"), 
          val_changes[[time_period]], shap_changes[[time_period]])
data <- project(data, plot_crs)
data <- as.data.frame(data, xy = TRUE) %>% na.omit()
names(data) <- c("x", "y", "area", "variable", "SHAP")

# data <- data %>% 
#     filter(SHAP >= quantile(SHAP, 0.01) & 
#                SHAP <= quantile(SHAP, 0.99))

if (all(data$SHAP <= 0)) {
    data[nrow(data) + 1, ] <- c(rep(NA, 3), 1, 1)
    data[nrow(data) + 1, ] <- c(rep(NA, 3), -1, 1)
}
if (all(data$SHAP > 0)) {
    data[nrow(data) + 1, ] <- c(rep(NA, 3), 1, -1)
    data[nrow(data) + 1, ] <- c(rep(NA, 2), -1, -1)
}
if (all(data$variable <= 0)) {
    data[nrow(data) + 1, ] <- c(rep(NA, 3), 1, 1)
    data[nrow(data) + 1, ] <- c(rep(NA, 3), 1, -1)
}
if (all(data$variable > 0)) {
    data[nrow(data) + 1, ] <- c(rep(NA, 3), -1, 1)
    data[nrow(data) + 1, ] <- c(rep(NA, 3), -1, -1)
}

data$SHAP <- ifelse(data$SHAP <= 0, 1, 2)
data$SHAP <- as.factor(data$SHAP)

data$variable <- ifelse(data$variable <= 0, 1, 2)
data$variable <- as.factor(data$variable)

data_classified <- bi_class(data, x = variable, y = SHAP, dim = 2)
breaks <- bi_class_breaks(data, x = variable, y = SHAP, 
                          dim = 2, dig_lab = 2, split = TRUE)
breaks$bi_x <- c("Drop", "Rise")
breaks$bi_y <- c("Drop", "Rise")

lgd <- ggplotGrob(
    bi_legend(
        pal = "DkBlue2",
        pad_color = "white",
        flip_axes = FALSE,
        rotate_pal = FALSE,
        dim = 2,
        breaks = breaks,
        xlab = nm,
        ylab = "SHAP change",
        pad_width = 0.5,
        size = 10) +
        theme(axis.text = element_text(
            size = 8, color = "black", family = 'Merriweather'),
            axis.title = element_text(
                size = 8, color = "black", family = 'Merriweather'),
            panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', 
                                           color = NA),
            axis.title.x = element_text(vjust = 4),
            axis.title.y = element_text(vjust = -7.5)))

# Create the bivariate map using ggplot2
ggplot() +
    geom_sf(data = bry, fill = "transparent", 
            color = "#d3d3d3", linewidth = 0.1) + 
    geom_tile(data = data_classified, 
              mapping = aes(x = x, y = y, fill = bi_class), 
              show.legend = FALSE) +
    bi_scale_fill(
        pal = "DkBlue2", dim = 2,
        flip_axes = FALSE, rotate_pal = FALSE) +
    annotation_custom(grob = lgd, 
                      xmin = -28108181, xmax = -6908181,
                      ymin = -15342702, ymax = 5557298) +
    annotation_custom(grob = ggplotGrob(tbl), 
                      xmin = -10108181, xmax = 18208181,
                      ymin = -13542702, ymax = 107298) +
    theme_void() + 
    theme(plot.margin = unit(c(-2.5, -0.4, -1, -0.6), "cm"))

ggsave(sprintf("docs/figures/%s_%s_shap_%s.png", 
               sus, driver, time_period),
       width = 6, height = 3.1, dpi = 500)
