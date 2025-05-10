# Load settings
source("R/figures/setting.R")

#### Define layer for low, middle and high latitude areas ###
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

#### Affected area ####
areas_periods <- lapply(c("High", "Middle", "Low", "Global"), function(lat){
    lapply(time_periods, function(time_period){
        do.call(rbind, lapply(features, function(driver){
            fnames <- list.files(
                cc_dir, pattern = sprintf("%s.tif", time_period), 
                full.names = TRUE)
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
            
            if (lat == "Global"){
                stous_area <- global(
                    mask(cellSize(stous, unit = "km"), stous),
                    "sum", na.rm = TRUE)[[1]]
                utoss_area <- global(
                    mask(cellSize(stous, unit = "km"), utoss),
                    "sum", na.rm = TRUE)[[1]]
            } else {
                # Load area
                msk <- areas %>% filter(area == lat)
                
                stous_area <- global(
                    mask(mask(cellSize(stous, unit = "km"), stous), msk),
                    "sum", na.rm = TRUE)[[1]]
                utoss_area <- global(
                    mask(mask(cellSize(stous, unit = "km"), utoss), msk),
                    "sum", na.rm = TRUE)[[1]]
            }
            
            data.frame(driver = driver, 
                       stou_percent = stous_area / all_area * 100,
                       utos_percent = utoss_area / all_area * 100) %>% 
                mutate(scenario = c("SSP126", "SSP370", "SSP585"))
        })) %>% mutate(time_period = time_period)
    }) %>% bind_rows() %>% mutate(area = lat)
}) %>% bind_rows()

write.csv(
    areas_periods, file.path(tbl_dir, "supporting_table_turnover_global.csv"), 
    row.names = FALSE)

#### Affected species ####
species_periods <- lapply(c("High", "Middle", "Low", "Global"), function(lat){
    mclapply(time_periods, function(time_period){
        do.call(rbind, lapply(features, function(driver){
            fnames <- list.files(
                cc_dir, pattern = sprintf("%s.tif", time_period), full.names = TRUE)
            fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
            
            dir_fnames <- fnames[str_detect(fnames, "dir_change")]
            num_fnames <- fnames[str_detect(fnames, "num_species")]
            
            # S to U
            stous <- do.call(c, lapply(1:length(dir_fnames), function(i){
                rast(dir_fnames[i])[[3]] / rast(num_fnames[i])[[1]]
            })); names(stous) <- c("SSP126", "SSP370", "SSP580")
            
            # U to S
            utoss <- do.call(c, lapply(1:length(dir_fnames), function(i){
                rast(dir_fnames[i])[[2]] / rast(num_fnames[i])[[1]]
            })); names(utoss) <- c("SSP126", "SSP370", "SSP580")
            
            if (lat != "Global"){
                # Load area
                msk <- areas %>% filter(area == lat)
                
                stous <- mask(stous, msk)
                utoss <- mask(utoss, msk)}
            
            # Trim the layers to remove unnecessary computation
            stous <- trim(stous)
            utoss <- trim(utoss)
            
            # Calculate the case of S to U
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
            
            # Calculate the case of U to S
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
    }, mc.cores = 3) %>% bind_rows() %>% mutate(area = lat)
}) %>% bind_rows()

write.csv(
    species_periods, file.path(tbl_dir, "supporting_table_turnover_local.csv"), 
    row.names = FALSE)

#### Make figures and tables ####
##### Supplementary Table 2 ####

areas_periods_tosave <- areas_periods %>% 
    select(driver, area, scenario, time_period, stou_percent, utos_percent) %>%
    mutate(driver = ifelse(
        driver == "forest", "FOR",
        ifelse(driver == "human_impact", "HLU",
               gsub("bio", "BIO", driver)))) %>%
    mutate(driver = ifelse(
        driver == "grassland", "GRA", driver)) %>% 
    mutate(area = ifelse(
        area == "Global", area, paste(area, "altitude"))) %>% as_tibble() %>% 
    rename("Variable" = driver, "Region" = area, "Scenario" = scenario,
           "Time period" = time_period, "Global area of P2N (%)" = stou_percent,
           "Global area of N2P (%)" = utos_percent)
write.csv(areas_periods_tosave, file.path(tbl_dir, "supplementary_table2.csv"), 
          row.names = FALSE)
rm(areas_periods_tosave)

##### Supplementary Table 3 ####
species_periods_tosave <- species_periods %>% 
    select(driver, area, scenario, time_period, stou_sp_mean, 
           stou_sp_sd, utos_sp_mean, utos_sp_sd) %>%
    mutate(driver = ifelse(
        driver == "forest", "FOR",
        ifelse(driver == "human_impact", "HLU",
               gsub("bio", "BIO", driver)))) %>%
    mutate(driver = ifelse(
        driver == "grassland", "GRA", driver)) %>% 
    mutate(area = ifelse(
        area == "Global", area, paste(area, "altitude"))) %>% as_tibble() %>% 
    rename("Variable" = driver, "Region" = area, "Scenario" = scenario,
           "Time period" = time_period, 
           "Local intensity of P2N (% of species, mean)" = stou_sp_mean,
           "Local intensity of of P2N (% of species, SD)" = stou_sp_sd,
           "Local intensity of N2P (% of species, mean)" = utos_sp_mean, 
           "Local intensity of N2P (% of species, SD)" = utos_sp_sd)
write.csv(species_periods_tosave, file.path(tbl_dir, "supplementary_table3.csv"), 
          row.names = FALSE)
rm(species_periods_tosave)

##### Extended Data Table 2 ####
# Across variable groups: temperature, precipitation and land cover
area_vals <- areas_periods %>% filter(area == "Global") %>% 
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
    mutate(type = "Global area") %>% 
    pivot_longer(c(8, 9), names_to = "turnover") %>% 
    select(type, turnover, group, scenario, time_period, value) %>% 
    pivot_wider(names_from = c(time_period, scenario), 
                values_from = value)

# Across variable groups: temperature, precipitation and land cover
species_vals <- species_periods %>% filter(area == "Global") %>% 
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
    mutate(type = "Local intensity") %>% 
    pivot_longer(c(8, 9), names_to = "turnover") %>% 
    select(type, turnover, group, scenario, time_period, value) %>% 
    pivot_wider(names_from = c(time_period, scenario), 
                values_from = value)

nms_in_order <- lapply(time_periods, function(x) {
    sprintf("%s_%s", x, c("SSP126", "SSP370", "SSP585"))}) %>% unlist()

rbind(area_vals, species_vals) %>% 
    arrange(type, turnover, group) %>% 
    select(all_of(c("type", "turnover", "group", nms_in_order))) %>% 
    rename("Type" = type, "Turnover" = turnover, 
           "Variable\ngroup" = group) %>% 
    flextable() %>% separate_header(split = "_") %>% 
    autofit() %>% align(align = "center", part = "all") %>% 
    font(fontname = "Merriweather", part = "all") %>% 
    bold(part = "header") %>% 
    bold(i = c(3, 6, 9, 10), j = 3:12) %>% 
    merge_at(i = 1:6, j = 1) %>% 
    merge_at(i = 1:3, j = 2) %>% merge_at(i = 4:6, j = 2) %>%
    merge_at(i = 7:12, j = 1) %>% 
    merge_at(i = 7:10, j = 2) %>% merge_at(i = 11:12, j = 2) %>% 
    # Add some lines
    border(i = c(3, 9), j = 3:12, 
           border.bottom = fp_border(color = "gray")) %>% 
    border(i = 6, j = 2:12, 
           border.bottom = fp_border(color = "gray")) %>%
    save_as_image(file.path(fig_dir, "extended_data_table2.png"), res = 500)

# Clean a bit
rm(area_vals, species_vals, nms_in_order); gc()

##### Figure 2 and Extended Data Fig.2-3 ####
figs <- lapply(c("SSP126", "SSP370", "SSP585"), function(ssp){
    figs <- lapply(time_periods, function(tp){
        if (ssp == "SSP370" & tp == "2041-2070"){
            titles <- c("(b) Global area", "(c) Area by region", 
                        "(d) Local intensity", "(e) Intensity by region")
        } else{
            titles <- c(
                sprintf("(%s %s)\nGlobal area", tolower(ssp), tp), 
                sprintf("(%s %s)\nArea by region", tolower(ssp), tp),
                sprintf("(%s %s)\nLocal intensity", tolower(ssp), tp), 
                sprintf("(%s %s)\nIntensity by region", tolower(ssp), tp))}
        
        if (tp == "2071-2100"){
            xlims1 <- c(-100, 100)
            xlims2 <- c(-65, 65)
        } else{
            xlims1 <- c(-95, 95)
            xlims2 <- c(-55, 55)
        } 
        
        area_pts_all <- areas_periods %>% 
            filter(scenario == ssp & time_period == tp) %>% 
            mutate(driver = ifelse(
                driver == "forest", "FOR",
                ifelse(driver == "human_impact", "HLU",
                       gsub("bio", "BIO", driver)))) %>% 
            mutate(driver = ifelse(
                driver == "grassland", "GRA", driver))
        
        # Get the driver rank for visualization
        drivers <- area_pts_all %>% filter(area == "Global") %>% 
            mutate(stou_rank = rank(-stou_percent), 
                   utos_rank = rank(-utos_percent)) %>% 
            mutate(SqRank = (stou_rank^2) + (utos_rank^2)/2) %>% 
            mutate(RankOrder = rank(SqRank)) %>% 
            arrange(-RankOrder)
        
        area_pts <- area_pts_all %>% filter(area == "Global") %>% 
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
            annotate(geom = "text", y = c(-80, 80), x = c(2, 2), 
                     label = c("\U2190N2P", "P2N\U2192"), fontface = "bold",
                     family = "Merriweather", color = "black", size = 2.5) +
            labs(x = "", y = "") + ggtitle(titles[1]) + coord_flip() + 
            scale_y_continuous(labels = function(x){paste0(abs(x), "%")},
                               limits = xlims1) +
            geom_hline(yintercept = 0, color = 'black', linewidth = 1) +
            theme_pubclean(base_family = "Merriweather", base_size = 11) +
            theme(panel.grid.major.y = element_line(color = "white"),
                  panel.grid.major.x = element_line(
                      linetype = "dotted", color = "lightgrey"),
                  axis.text.y = element_text(color = "black"),
                  plot.title = element_text(
                      family = "Merriweather", size = 11,
                      face = "bold", hjust = 0.5),
                  plot.margin = unit(c(0, 0, -0.2, 0), "cm"))
        
        # By regions
        area_regions_pts <- area_pts_all %>% filter(area != "Global") %>% 
            mutate(driver = factor(
                driver, levels = drivers$driver, labels = drivers$driver)) %>%
            mutate(area_p2n = paste0(area, "_p2n"), area_n2p = area) %>% 
            mutate(area_n2p = factor(
                area_n2p, levels = c("Low", "Middle", "High"))) %>% 
            mutate(area_p2n = factor(
                area_p2n, levels = c("Low_p2n", "Middle_p2n", "High_p2n")))
        
        g2 <- ggplot(area_regions_pts) +
            geom_tile(aes(x = area_n2p, y = driver, fill = utos_percent),
                      color = "white") + 
            scale_fill_bs5(name = "N2P turnover\n(%)", "green",
                           guide = guide_colorbar(order = 1)) +
            new_scale_fill() +
            geom_tile(aes(x = area_p2n, y = driver, fill = stou_percent),
                      color = "white") + 
            scale_fill_bs5(name = "P2N turnover\n(%)", "orange",
                           guide = guide_colorbar(order = 2)) +
            scale_x_discrete(labels = rep(c("L", "M", "H"), 2)) +
            geom_vline(xintercept = 3.5, color = 'black', linewidth = 1) +
            labs(x = "", y = "") + ggtitle(titles[2]) + coord_equal() +
            theme_pubclean(base_family = "Merriweather", base_size = 11) +
            theme(panel.grid.major.y = element_line(color = "white"),
                  axis.text.y = element_text(color = "black"),
                  legend.position = "right", 
                  legend.direction = "vertical",
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 8),
                  legend.key.height = unit(3, "mm"),
                  legend.key.width = unit(2, "mm"),
                  plot.title = element_text(
                      family = "Merriweather", size = 11,
                      face = "bold", hjust = 0.5),
                  axis.text.x = element_text(color = "black"),
                  plot.margin = unit(c(0, 0, -0.2, 0), "cm"))
        
        g1 <- ggarrange(g1, g2, ncol = 2)
        
        # Species
        species_pts_all <- species_periods %>% 
            filter(scenario == ssp & time_period == tp) %>% 
            mutate(driver = ifelse(
                driver == "forest", "FOR",
                ifelse(driver == "human_impact", "HLU",
                       gsub("bio", "BIO", driver)))) %>% 
            mutate(driver = ifelse(
                driver == "grassland", "GRA", driver))
        
        # Get the driver rank for visualization
        drivers <- species_pts_all %>% filter(area == "Global") %>% 
            mutate(stou_rank = rank(-stou_sp_mean), 
                   utos_rank = rank(-utos_sp_mean)) %>% 
            mutate(SqRank = (stou_rank^2) + (utos_rank^2) / 2) %>% 
            mutate(RankOrder = rank(SqRank)) %>% 
            arrange(-RankOrder)
        
        species_pts <- species_pts_all %>% filter(area == "Global") %>% 
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
            annotate(geom = "text", y = c(-45, 45), x = c(2, 2), 
                     label = c("\U2190N2P", "P2N\U2192"), fontface = "bold",
                     family = "Merriweather", color = "black", size = 2.5) +
            labs(x = "", y = "") + ggtitle(titles[3]) + coord_flip() + 
            scale_y_continuous(labels = function(x){paste0(abs(x), "%")},
                               limits = xlims2) +
            geom_hline(yintercept = 0, color = 'black', linewidth = 1) +
            theme_pubclean(base_family = "Merriweather", base_size = 11) +
            theme(panel.grid.major.y = element_line(color = "white"),
                  panel.grid.major.x = element_line(
                      linetype = "dotted", color = "lightgrey"),
                  axis.text.y = element_text(color = "black"),
                  plot.title = element_text(
                      family = "Merriweather", size = 11,
                      face = "bold", hjust = 0.5),
                  plot.margin = unit(c(0, 0, -0.2, 0), "cm"))
        
        # By regions
        species_regions_pts <- species_pts_all %>% filter(area != "Global") %>% 
            mutate(driver = factor(
                driver, levels = drivers$driver, labels = drivers$driver)) %>%
            mutate(area_p2n = paste0(area, "_p2n"), area_n2p = area) %>% 
            mutate(area_n2p = factor(
                area_n2p, levels = c("Low", "Middle", "High"))) %>% 
            mutate(area_p2n = factor(
                area_p2n, levels = c("Low_p2n", "Middle_p2n", "High_p2n")))
        
        g3 <- ggplot(species_regions_pts) +
            geom_tile(aes(x = area_n2p, y = driver, fill = stou_sp_mean),
                      color = "white") + 
            scale_fill_bs5(name = "N2P turnover\n(%)", "green",
                           guide = guide_colorbar(order = 1)) +
            new_scale_fill() +
            geom_tile(aes(x = area_p2n, y = driver, fill = utos_sp_mean),
                      color = "white") + 
            scale_fill_bs5(name = "P2N turnover\n(%)", "orange",
                           guide = guide_colorbar(order = 2)) +
            scale_x_discrete(labels = rep(c("L", "M", "H"), 2)) +
            geom_vline(xintercept = 3.5, color = 'black', linewidth = 1) +
            labs(x = "", y = "") + ggtitle(titles[4]) + coord_equal() +
            theme_pubclean(base_family = "Merriweather", base_size = 11) +
            theme(panel.grid.major.y = element_line(color = "white"),
                  axis.text.y = element_text(color = "black"),
                  legend.position = "right", 
                  legend.direction = "vertical",
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 8),
                  legend.key.height = unit(3, "mm"),
                  legend.key.width = unit(2, "mm"),
                  plot.title = element_text(
                      family = "Merriweather", size = 11,
                      face = "bold", hjust = 0.5),
                  axis.text.x = element_text(color = "black"),
                  plot.margin = unit(c(0, 0, 0, 0), "cm"))
        
        g2 <- ggarrange(g2, g3, ncol = 2)
        
        list("area" = g1, "intensity" = g2)
    })
    
    # Add name
    names(figs) <- time_periods
    figs
}); names(figs) <- ssps

###### Figure 2 ####

img <- image_read(file.path(fig_dir, "fig2_flow.png"))
g <- image_ggplot(img, interpolate = TRUE)

ggarrange(g, figs[[2]][[2]][[1]], figs[[2]][[2]][[2]], 
          nrow = 3, heights = c(1.4, 3, 3))

ggsave(file.path(fig_dir, "Figure2_global_patterns.png"), 
       width = 6.5, height = 7, dpi = 500, bg = "white")

###### Extended Data Fig.2 ####

figs <- do.call(c, figs)

plotlist <- lapply(names(figs), function(x) {
    if (x == "ssp370.2041-2070"){
        NULL
    } else figs[[x]]$area})
ggarrange(plotlist = plotlist, nrow = 3, ncol = 3)

ggsave(file.path(fig_dir, "extended_data_fig2.png"), 
       width = 18, height = 11, dpi = 500, bg = "white")

###### Extended Data Fig.3 ####

plotlist <- lapply(names(figs), function(x) {
    if (x == "ssp370.2041-2070"){
        NULL
    } else figs[[x]]$intensity})
ggarrange(plotlist = plotlist, nrow = 3, ncol = 3)

ggsave(file.path(fig_dir, "extended_data_fig3.png"), 
       width = 18, height = 11, dpi = 500, bg = "white")
