# Load settings
source("R/figures/setting.R")

#### Affected area ####
areas_periods <- lapply(c("High", "Middle", "Low", "Global"), function(lat){
    mclapply(time_periods, function(time_period){
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
            })); names(stous) <- c("SSP126", "SSP370", "SSP585")
            stous <- trim(stous)
            
            # U to S
            utoss <- do.call(c, lapply(1:length(dir_fnames), function(i){
                rast(dir_fnames[i])[[2]] / rast(num_fnames[i])[[1]]
            })); names(utoss) <- c("SSP126", "SSP370", "SSP585")
            utoss <- trim(utoss)
            
            all_area <- global(
                mask(cellSize(stous, unit = "km"), stous), 
                "sum", na.rm = TRUE)[[1]]
            
            stous[stous <= 0] <- NA
            utoss[utoss <= 0] <- NA
            
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
    }, mc.cores = 3) %>% bind_rows() %>% mutate(area = lat)
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
            })); names(stous) <- c("SSP126", "SSP370", "SSP585")
            
            # U to S
            utoss <- do.call(c, lapply(1:length(dir_fnames), function(i){
                rast(dir_fnames[i])[[2]] / rast(num_fnames[i])[[1]]
            })); names(utoss) <- c("SSP126", "SSP370", "SSP585")
            
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
            
            ssp126 <- stous_sp %>% select(area, SSP126) %>% filter(SSP126 > 0)
            ssp370 <- stous_sp %>% select(area, SSP370) %>% filter(SSP370 > 0)
            ssp585 <- stous_sp %>% select(area, SSP585) %>% filter(SSP585 > 0)
            
            stous_1q <- c(
                weighted.quantile(ssp126$SSP126, ssp126$area, 0.25),
                weighted.quantile(ssp370$SSP370, ssp370$area, 0.25),
                weighted.quantile(ssp585$SSP585, ssp585$area, 0.25))
            
            stous_3q <- c(
                weighted.quantile(ssp126$SSP126, ssp126$area, 0.75),
                weighted.quantile(ssp370$SSP370, ssp370$area, 0.75),
                weighted.quantile(ssp585$SSP585, ssp585$area, 0.75))
            
            stous_medians <- c(
                weighted.median(ssp126$SSP126, ssp126$area, type = 4),
                weighted.median(ssp370$SSP370, ssp370$area, type = 4),
                weighted.median(ssp585$SSP585, ssp585$area, type = 4))
            
            # Calculate the case of U to S
            utoss_sp <- c(cellSize(utoss, unit = "km"), utoss)
            utoss_sp <- values(utoss_sp) %>% na.omit() %>% as.data.frame()
            
            ssp126 <- utoss_sp %>% select(area, SSP126) %>% filter(SSP126 > 0)
            ssp370 <- utoss_sp %>% select(area, SSP370) %>% filter(SSP370 > 0)
            ssp585 <- utoss_sp %>% select(area, SSP585) %>% filter(SSP585 > 0)
            
            utoss_1q <- c(
                weighted.quantile(ssp126$SSP126, ssp126$area, 0.25),
                weighted.quantile(ssp370$SSP370, ssp370$area, 0.25),
                weighted.quantile(ssp585$SSP585, ssp585$area, 0.25))
            
            utoss_3q <- c(
                weighted.quantile(ssp126$SSP126, ssp126$area, 0.75),
                weighted.quantile(ssp370$SSP370, ssp370$area, 0.75),
                weighted.quantile(ssp585$SSP585, ssp585$area, 0.75))
            
            utoss_medians <- c(
                weighted.median(ssp126$SSP126, ssp126$area, type = 4),
                weighted.median(ssp370$SSP370, ssp370$area, type = 4),
                weighted.median(ssp585$SSP585, ssp585$area, type = 4))
            
            data.frame(driver = driver, 
                       stou_sp_median = stous_medians * 100,
                       stou_sp_1q = stous_1q * 100,
                       stou_sp_3q = stous_3q * 100,
                       utos_sp_median = utoss_medians * 100,
                       utos_sp_1q = utoss_1q * 100,
                       utos_sp_3q = utoss_3q * 100) %>% 
                mutate(scenario = c("SSP126", "SSP370", "SSP585"))
        })) %>% mutate(time_period = time_period)
    }, mc.cores = 3) %>% bind_rows() %>% mutate(area = lat)
}) %>% bind_rows()

write.csv(
    species_periods, file.path(tbl_dir, "supporting_table_turnover_local.csv"), 
    row.names = FALSE)

#### Make figures and tables ####
##### Supplementary Table 3 ####

areas_periods_tosave <- areas_periods %>% 
    select(driver, area, scenario, time_period, stou_percent, utos_percent) %>%
    mutate(driver = case_when(
        driver == "forest" ~ "FOR",
        driver == "human_impact" ~ "HLU",
        driver == "grassland" ~ "GRA",
        str_detect(driver, "bio") ~ toupper(driver))) %>% 
    mutate(area = ifelse(
        area == "Global", area, paste(area, "altitude"))) %>% as_tibble() %>% 
    rename("Variable" = driver, "Region" = area, "Scenario" = scenario,
           "Time period" = time_period, "Global area of P2N (%)" = stou_percent,
           "Global area of N2P (%)" = utos_percent)
write.csv(areas_periods_tosave, file.path(tbl_dir, "supplementary_table3.csv"), 
          row.names = FALSE)
rm(areas_periods_tosave)

##### Supplementary Table 4 ####
species_periods_tosave <- species_periods %>% 
    select(driver, area, scenario, time_period, stou_sp_median, 
           stou_sp_1q, stou_sp_3q, utos_sp_median, utos_sp_1q, utos_sp_3q) %>%
    mutate(driver = case_when(
        driver == "forest" ~ "FOR",
        driver == "human_impact" ~ "HLU",
        driver == "grassland" ~ "GRA",
        str_detect(driver, "bio") ~ toupper(driver))) %>% 
    mutate(area = ifelse(
        area == "Global", area, paste(area, "altitude"))) %>% as_tibble() %>% 
    rename("Variable" = driver, "Region" = area, "Scenario" = scenario,
           "Time period" = time_period, 
           "Local intensity of P2N\n(% of species, median)" = stou_sp_median,
           "Local intensity of of P2N\n(% of species, Q1)" = stou_sp_1q,
           "Local intensity of of P2N\n(% of species, Q3)" = stou_sp_3q,
           "Local intensity of N2P\n(% of species, median)" = utos_sp_median, 
           "Local intensity of N2P\n(% of species, Q1)" = utos_sp_1q,
           "Local intensity of N2P\n(% of species, Q3)" = utos_sp_3q)
write.csv(species_periods_tosave, file.path(tbl_dir, "supplementary_table4.csv"), 
          row.names = FALSE)
rm(species_periods_tosave)

##### Extended Data Table 3 ####
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
    summarise(stou_sd = sd(stou_sp_median), 
              stou_percent = mean(stou_sp_median),
              utos_sd = sd(utos_sp_median),
              utos_percent = mean(utos_sp_median),
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
    bg(bg = "white", part = "all") %>% 
    autofit() %>% align(align = "center", part = "all") %>% 
    font(fontname = "Merriweather", part = "all") %>% 
    bold(part = "header") %>% 
    bold(i = c(3, 6), j = 3:12) %>% 
    bold(i = 7, j = c(4, 10)) %>% bold(i = 9, j = c(5:9, 11:12)) %>% 
    bold(i = 10, j = c(4:8, 10)) %>% bold(i = 12, j = c(9, 11:12)) %>% 
    merge_at(i = 1:6, j = 1) %>% 
    merge_at(i = 1:3, j = 2) %>% merge_at(i = 4:6, j = 2) %>%
    merge_at(i = 7:12, j = 1) %>% 
    merge_at(i = 7:10, j = 2) %>% merge_at(i = 11:12, j = 2) %>% 
    # Add some lines
    border(i = c(3, 9), j = 3:12, 
           border.bottom = fp_border(color = "gray")) %>% 
    border(i = 6, j = 2:12, 
           border.bottom = fp_border(color = "gray")) %>%
    save_as_image(file.path(fig_dir, "extended_data_table3.png"), res = 500)

# Clean a bit
rm(area_vals, species_vals, nms_in_order); gc()

##### Figure 2 and Extended Data Fig.3 ####
figs <- lapply(c("SSP126", "SSP370", "SSP585"), function(ssp){
    figs <- lapply(time_periods, function(tp){
        if (ssp == "SSP370" & tp == "2041-2070"){
            titles <- c("(b) Global area", "(c) Local intensity")
        } else{
            titles <- c(
                sprintf("(%s %s)\nGlobal area", tolower(ssp), tp), 
                sprintf("(%s %s)\nLocal intensity", tolower(ssp), tp))}
        
        if (tp == "2071-2100"){
            xlims1 <- c(-100, 101)
            xlims2 <- c(-65, 68)
        } else{
            xlims1 <- c(-95, 95)
            xlims2 <- c(-55, 55)
        }
        
        if (ssp == "SSP126" & tp == "2041-2070"){
            xlims2 <- c(-55, 68)
        }
        
        if (ssp == "SSP585" & tp == "2071-2100"){
            xlims1 <- c(-100, 105)
        }
        
        if (ssp == "SSP585" & tp == "2041-2070"){
            xlims1 <- c(-100, 100)
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
        a_drivers <- area_pts_all %>% filter(area == "Global") %>% 
            mutate(stou_rank = rank(-stou_percent), 
                   utos_rank = rank(-utos_percent)) %>% 
            mutate(SqRank = (stou_rank^2) + (utos_rank^2)/2) %>% 
            mutate(RankOrder = rank(SqRank)) %>% 
            arrange(-RankOrder)
        
        area_pts <- area_pts_all %>% filter(area == "Global") %>% 
            mutate(driver = factor(
                driver, levels = a_drivers$driver, labels = a_drivers$driver))
        
        # Species
        species_pts_all <- species_periods %>% 
            filter(scenario == ssp & time_period == tp) %>% 
            mutate(driver = case_when(
                driver == "forest" ~ "FOR",
                driver == "human_impact" ~ "HLU",
                driver == "grassland" ~ "GRA",
                str_detect(driver, "bio") ~ toupper(driver)))
        
        # Get the driver rank for visualization
        s_drivers <- species_pts_all %>% filter(area == "Global") %>% 
            mutate(stou_rank = rank(-stou_sp_median), 
                   utos_rank = rank(-utos_sp_median)) %>% 
            mutate(SqRank = (stou_rank^2) + (utos_rank^2) / 2) %>% 
            mutate(RankOrder = rank(SqRank)) %>% 
            arrange(-RankOrder)
        
        species_pts <- species_pts_all %>% filter(area == "Global") %>% 
            mutate(driver = factor(
                driver, levels = s_drivers$driver, labels = s_drivers$driver))
        
        # Get the changes between these two
        ord_change <- data.frame(
            area = a_drivers$driver, intensity = s_drivers$driver) %>% 
            mutate(area_order = 1:nrow(.))
        
        a_cols <- lapply(1:nrow(ord_change), function(i){
            ord_change %>% slice(i) %>% 
                mutate(intensity_order = which(ord_change$intensity == .$area))
        }) %>% bind_rows() %>% 
            mutate(change = intensity_order - area_order)
        
        a_cols <- ifelse(a_cols$change >= 0, "#1d3557", "#e63946")
        
        # Plots
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
                      color = 'black', size = 2.5, family = "Merriweather", 
                      hjust = 0) +
            # text for unsuitable to suitable
            geom_text(data = area_pts, 
                      aes(x = driver, y = -utos_percent - 1, 
                          label = sprintf("%.0f%s", utos_percent, "%")), 
                      color = 'black', size = 2.5, family = "Merriweather", 
                      hjust = 1) +
            annotate(geom = "text", y = c(xlims1[1] * 0.8, xlims1[2] * 0.8), 
                     x = c(2, 2), 
                     label = c("\U2190N2P", "P2N\U2192"), fontface = "bold",
                     family = "Merriweather", color = "black", size = 2.5) +
            labs(x = "", y = "% of terrestrial area") + 
            ggtitle(titles[1]) + coord_flip() + 
            scale_y_continuous(labels = function(x){paste0(abs(x), "%")},
                               limits = xlims1) +
            scale_x_discrete(position = "top") +
            geom_hline(yintercept = 0, color = 'black', linewidth = 1) +
            theme_pubclean(base_family = "Merriweather", base_size = 11) +
            theme(panel.grid.major.y = element_line(color = "white"),
                  panel.grid.major.x = element_line(
                      linetype = "dotted", color = "lightgrey"),
                  axis.text.x = element_text(color = "black"),
                  axis.text.y = element_text(color = a_cols),
                  plot.title = element_text(
                      family = "Merriweather", size = 11,
                      face = "bold", hjust = 0.5),
                  plot.margin = unit(c(0.1, -0.3, 0.1, 0.5), "cm"))
        
        s_cols <- lapply(1:nrow(ord_change), function(i){
            ord_change %>% slice(i) %>% 
                mutate(intensity_order = which(ord_change$area == .$intensity))
        }) %>% bind_rows() %>% 
            mutate(change = intensity_order - area_order)
        
        s_cols <- ifelse(s_cols$change > 0, "#e63946", "#1d3557")
        
        g2 <- ggplot() + 
            geom_col(data = species_pts, 
                     color = "white", fill = "#a6611a",
                     aes(x = driver, y = stou_sp_3q)) +
            geom_col(data = species_pts, 
                     color = "white", fill = "white",
                     aes(x = driver, y = stou_sp_1q)) +
            geom_errorbar(data = species_pts, 
                          aes(x = driver, ymin = stou_sp_median, 
                              ymax = stou_sp_median), width = 0.8) +
            geom_col(data = species_pts,
                     aes(x = driver, y = -utos_sp_3q), 
                     color = "white", fill = "#018571") +
            geom_col(data = species_pts,
                     aes(x = driver, y = -utos_sp_1q), 
                     color = "white", fill = "white") +
            geom_errorbar(data = species_pts, 
                          aes(x = driver, ymin = -utos_sp_median, 
                              ymax = -utos_sp_median), width = 0.8) +
            # text for suitable to unsuitable
            geom_text(data = species_pts, 
                      aes(x = driver, y = stou_sp_3q + 1, 
                          label = sprintf("%.0f-%.0f%s", 
                                          stou_sp_1q, stou_sp_3q, "%")), 
                      color = 'black', size = 2.5, family = "Merriweather", 
                      hjust = 0) +
            # text for unsuitable to suitable
            geom_text(data = species_pts, 
                      aes(x = driver, y = -utos_sp_3q - 1, 
                          label = sprintf("%.0f-%.0f%s", 
                                          utos_sp_1q, utos_sp_3q, "%")), 
                      color = 'black', size = 2.5, family = "Merriweather", 
                      hjust = 1) +
            annotate(geom = "text", y = c(xlims2[1] * 0.8, xlims2[2] * 0.8), 
                     x = c(2, 2), 
                     label = c("\U2190N2P", "P2N\U2192"), fontface = "bold",
                     family = "Merriweather", color = "black", size = 2.5) +
            labs(x = "", y = "% of species") + 
            ggtitle(titles[2]) + coord_flip() + 
            scale_y_continuous(labels = function(x){paste0(abs(x), "%")},
                               limits = xlims2) +
            geom_hline(yintercept = 0, color = 'black', linewidth = 1) +
            theme_pubclean(base_family = "Merriweather", base_size = 11) +
            theme(panel.grid.major.y = element_line(color = "white"),
                  panel.grid.major.x = element_line(
                      linetype = "dotted", color = "lightgrey"),
                  axis.text.x = element_text(color = "black"),
                  axis.text.y = element_text(color = s_cols),
                  plot.title = element_text(
                      family = "Merriweather", size = 11,
                      face = "bold", hjust = 0.5),
                  plot.margin = unit(c(0.1, 0.5, 0.1, -0.2), "cm"))
        
        dirs <- lapply(1:nrow(ord_change), function(i){
            ord_change %>% slice(i) %>% 
                mutate(intensity_order = which(ord_change$intensity == .$area))
        }) %>% bind_rows() %>% 
            mutate(area = factor(
                area, levels = area, labels = area))
        
        d_cols <- ifelse(a_cols == "black", "#1d3557", a_cols)
        
        g3 <- ggplot(dirs) + 
            geom_point(aes(x = 0, y = area_order), 
                       color = a_cols, size = 0.6) +
            geom_point(aes(x = 1, y = intensity_order), 
                       color = a_cols, size = 0.6) +
            ggtitle(titles[1]) + labs(y = "", x = " ") + 
            geom_segment(
                aes(x = 0, y = area_order, 
                    xend = 1, yend = intensity_order), 
                color = d_cols) + 
            theme_pubclean(base_family = "Merriweather", base_size = 11) +
            theme(panel.grid.major.y = element_line(color = "white"),
                  panel.grid.major.x = element_line(color = "white"),
                  axis.text.x = element_text(color = "white"),
                  axis.ticks = element_blank(),
                  axis.text.y = element_blank(),
                  plot.title = element_text(
                      family = "Merriweather", size = 11,
                      face = "bold", hjust = 0.5, color = "white"),
                  plot.margin = unit(c(0, 0, 0, -0.4), "cm"))
        
        ggarrange(g1, g3, g2, nrow = 1, widths = c(4, 0.8, 4))
    })
    
    # Add name
    names(figs) <- time_periods
    figs
}); names(figs) <- ssps

###### Figure 2 ####

img <- image_read(file.path(fig_dir, "fig2_flow.png"))
g <- image_ggplot(img, interpolate = TRUE)

ggarrange(g, NULL, figs[[2]][[2]], 
          nrow = 3, heights = c(1.5, 0.1, 5),
          labels = c("\n(a)", "", ""),
          font.label = list(
              size = 11, color = "black", 
              face = "bold", family = "Merriweather"))

ggsave(file.path(fig_dir, "Figure2_global_patterns.png"), 
       width = 6.5, height = 4.6, dpi = 500, bg = "white")

###### Extended Data Fig.3 ####

figs <- do.call(c, figs)

plotlist <- lapply(names(figs), function(x) {
    if (x == "ssp370.2041-2070"){
        NULL
    } else figs[[x]]})
ggarrange(plotlist = plotlist, nrow = 3, ncol = 3)

ggsave(file.path(fig_dir, "extended_data_fig3.png"), 
       width = 18, height = 11, dpi = 500, bg = "white")

# Clean up
rm(list = ls()); gc()
