# Load settings
source("R/figures/setting.R")

##### Affected area ####
areas_periods <- lapply(time_periods, function(time_period){
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

write.csv(areas_periods, file.path(tbl_dir, "affect_area.csv"), 
          row.names = FALSE)

# Across variable groups: temperature, precipitation and land cover
tbl_vals <- areas_periods %>% 
    mutate(group = as.integer(str_extract(driver, "[0-9]+"))) %>% 
    mutate(group = ifelse(is.na(group), "Landcover", 
                          ifelse (group < 12, "Temperature", 
                                  "Precipitation"))) %>% 
    group_by(scenario, time_period, group) %>% 
    summarise(stou_sd = sd(stou_percent), 
              stou_percent = mean(stou_percent),
              utos_sd = sd(utos_percent),
              utos_percent = mean(utos_percent),
              .groups = "drop") %>% 
    mutate(P2N = sprintf("%.1f\u00B1%.1f", stou_percent, stou_sd),
           N2P = sprintf("%.1f\u00B1%.1f", utos_percent, utos_sd)) %>% 
    select(scenario, time_period, group, P2N, N2P) %>% 
    arrange(time_period) %>% 
    pivot_wider(names_from = c(scenario, time_period), 
                values_from = c(P2N, N2P))

write.csv(tbl_vals, file.path(tbl_dir, "affect_area_groups.csv"), 
          row.names = FALSE, fileEncoding = "UTF-16LE")

##### Affected species ####
species_periods <- lapply(time_periods, function(time_period){
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

write.csv(species_periods, file.path(tbl_dir, "affect_species.csv"), 
          row.names = FALSE)

# Across variable groups: temperature, precipitation and land cover
tbl_vals <- species_periods %>% 
    mutate(group = as.integer(str_extract(driver, "[0-9]+"))) %>% 
    mutate(group = ifelse(is.na(group), "Landcover", 
                          ifelse (group < 12, "Temperature", 
                                  "Precipitation"))) %>% 
    group_by(scenario, time_period, group) %>% 
    summarise(stou_sd = sd(stou_sp_mean), 
              stou_percent = mean(stou_sp_mean),
              utos_sd = sd(utos_sp_mean),
              utos_percent = mean(utos_sp_mean),
              .groups = "drop") %>% 
    mutate(P2N = sprintf("%.1f\u00B1%.1f", stou_percent, stou_sd),
           N2P = sprintf("%.1f\u00B1%.1f", utos_percent, utos_sd)) %>% 
    select(scenario, time_period, group, P2N, N2P) %>% 
    arrange(time_period) %>% 
    pivot_wider(names_from = c(scenario, time_period), 
                values_from = c(P2N, N2P))

write.csv(tbl_vals, file.path(tbl_dir, "affect_species_groups.csv"),
          row.names = FALSE, fileEncoding = "UTF-16LE")

# Clean a bit
rm(tbl_vals); gc()

for (ssp in c("SSP126", "SSP370", "SSP585")){
    for(tp in time_periods){
        area_pts <- areas_periods %>% 
            filter(scenario == ssp & time_period == tp) %>% 
            mutate(driver = ifelse(
                driver == "forest", "FOR",
                ifelse(driver == "human_impact", "HLU",
                       gsub("bio", "BIO", driver)))) %>% 
            mutate(driver = ifelse(
                driver == "grassland", "GRA", driver))
        
        # Get the driver rank for visualization
        drivers <- area_pts %>% 
            mutate(stou_rank = rank(-stou_percent), 
                   utos_rank = rank(-utos_percent)) %>% 
            mutate(SqRank = (stou_rank^2) + (utos_rank^2)/2) %>% 
            mutate(RankOrder = rank(SqRank)) %>% 
            arrange(-RankOrder)
        
        area_pts <- area_pts %>% 
            mutate(driver = factor(
                driver, levels = drivers$driver, labels = drivers$driver))
        
        g1 <- ggplot() + 
            geom_col(data = area_pts, 
                     color = "white", fill = "#a6611a",
                     aes(x = driver, y = stou_percent)) +
            geom_col(data = area_pts,
                     aes(x = driver, y = -utos_percent), 
                     color = "white", fill = "#018571") +
            # text for suitable to unsuitable
            geom_text(data = area_pts, 
                      aes(x = driver, y = stou_percent + 1, 
                          label = sprintf("%.0f%s", stou_percent, "%")), 
                      color = 'black', size = 2.7, family = "Merriweather", 
                      hjust = 0) +
            # text for unsuitable to suitable
            geom_text(data = area_pts, 
                      aes(x = driver, y = -utos_percent - 1, 
                          label = sprintf("%.0f%s", utos_percent, "%")), 
                      color = 'black', size = 2.7, family = "Merriweather", 
                      hjust = 1) +
            labs(x = "", y = "") + coord_flip() + ylim(-95, 95) +
            geom_hline(yintercept = 0, color = 'black', linewidth = 1) +
            theme_pubclean(base_family = "Merriweather", base_size = 11) +
            theme(axis.text.x = element_blank(),
                  panel.grid.major.y = element_line(color = "white"),
                  axis.ticks.x = element_blank(),
                  axis.text.y = element_text(color = "black"),
                  plot.margin = unit(c(0, 0, -0.5, 0), "cm"))
        
        # Species
        species_pts <- species_periods %>% 
            filter(scenario == ssp & time_period == tp) %>% 
            mutate(driver = ifelse(
                driver == "forest", "FOR",
                ifelse(driver == "human_impact", "HLU",
                       gsub("bio", "BIO", driver)))) %>% 
            mutate(driver = ifelse(
                driver == "grassland", "GRA", driver))
        
        # Get the driver rank for visualization
        drivers <- species_pts %>% 
            mutate(stou_rank = rank(-stou_sp_mean), 
                   utos_rank = rank(-utos_sp_mean)) %>% 
            mutate(SqRank = (stou_rank^2) + (utos_rank^2) / 2) %>% 
            mutate(RankOrder = rank(SqRank)) %>% 
            arrange(-RankOrder)
        
        species_pts <- species_pts %>% 
            mutate(driver = factor(
                driver, levels = drivers$driver, labels = drivers$driver))
        
        g2 <- ggplot() + 
            geom_col(data = species_pts, 
                     color = "white", fill = "#a6611a",
                     aes(x = driver, y = stou_sp_mean)) +
            geom_errorbar(data = species_pts, 
                          aes(x = driver, ymin = stou_sp_mean - stou_sp_sd, 
                              ymax = stou_sp_mean + stou_sp_sd), 
                          width = 0.3) +
            geom_col(data = species_pts,
                     aes(x = driver, y = -utos_sp_mean), 
                     color = "white", fill = "#018571") +
            geom_errorbar(data = species_pts, 
                          aes(x = driver, ymin = -utos_sp_mean + utos_sp_sd, 
                              ymax = -utos_sp_mean - utos_sp_sd), 
                          width = 0.3)+
            # text for suitable to unsuitable
            geom_text(data = species_pts, 
                      aes(x = driver, y = stou_sp_mean + stou_sp_sd + 1, 
                          label = sprintf("%.0f\u00B1%.0f%s", 
                                          stou_sp_mean, stou_sp_sd, "%")), 
                      color = 'black', size = 2.7, family = "Merriweather", 
                      hjust = 0) +
            # text for unsuitable to suitable
            geom_text(data = species_pts, 
                      aes(x = driver, y = -utos_sp_mean - utos_sp_sd - 1, 
                          label = sprintf("%.0f\u00B1%.0f%s", 
                                          utos_sp_mean, utos_sp_sd, "%")), 
                      color = 'black', size = 2.7, family = "Merriweather", 
                      hjust = 1) +
            labs(x = "", y = "") + coord_flip() + ylim(-50, 50) +
            geom_hline(yintercept = 0, color = 'black', linewidth = 1) +
            theme_pubclean(base_family = "Merriweather", base_size = 11) +
            theme(axis.text.x = element_blank(),
                  panel.grid.major.y = element_line(color = "white"),
                  axis.ticks.x = element_blank(),
                  axis.text.y = element_text(color = "black"),
                  plot.margin = unit(c(0, 0, -0.5, 0), "cm"))
        
        img <- image_read(file.path(fig_dir, "fig2_flow.png"))
        g <- image_ggplot(img, interpolate = TRUE)
        
        ann1 <- ggplot() + 
            geom_text(aes(x = 0, y = 0, label = "bold(b)"), 
                      parse = TRUE, size = 4.2, hjust = -2.4) +
            theme_void(base_family = "Merriweather")
        
        ann2 <- ggplot() + 
            geom_text(aes(x = 0, y = 0, label = "bold(c)"), 
                      parse = TRUE, size = 4.2, hjust = -2.4) +
            theme_void(base_family = "Merriweather")
        
        ggarrange(
            g, ggarrange(ann1, ann2, ncol = 2), ggarrange(g1, g2, ncol = 2), 
            nrow = 3, heights = c(1, 0.2, 2.8), labels = c("a", "", ""),
            vjust = 4, 
            font.label = list(
                size = 11, color = "black", 
                face = "bold", family = "Merriweather"))
        
        fname <- ifelse(ssp == "SSP370", "docs/figures/Figure1_turnover_area.png",
                        sprintf("docs/figures/Figure_s_%s_turnover_area.png", ssp))
        
        fname <- "/Users/leisong/downloads/test.png"
        ggsave(fname, width = 6.5, height = 4, dpi = 500, bg = "white")
    }
}
