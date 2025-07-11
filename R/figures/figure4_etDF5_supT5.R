# Load settings
source("R/figures/setting.R")

#### SHAP magnitude shifts ####
magnitude_periods <- lapply(
    c("High", "Middle", "Low", "Global"), function(lat){
    lapply(time_periods, function(time_period){
        do.call(rbind, lapply(features, function(driver){
            fnames <- list.files(
                cc_dir, pattern = sprintf("%s.tif", time_period), 
                full.names = TRUE)
            fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
            
            val_fnames <- fnames[
                str_detect(fnames, sprintf("val_change_%s", driver))]
            num_fnames <- fnames[str_detect(fnames, "num_species")]
            
            # Suitable
            ss <- do.call(c, lapply(1:length(val_fnames), function(i){
                msk <- rast(num_fnames[i])[[2]]
                msk[msk == 0] <- NA
                rast(val_fnames[i])[[1]] / msk
            })); names(ss) <- c("SSP126", "SSP370", "SSP585")
            ss <- trim(ss)
            
            # Unsuitable
            uss <- do.call(c, lapply(1:length(val_fnames), function(i){
                msk <- rast(num_fnames[i])[[3]]
                msk[msk == 0] <- NA
                rast(val_fnames[i])[[2]] / msk
            })); names(uss) <- c("SSP126", "SSP370", "SSP585")
            uss <- trim(uss)
            
            ss_val <- c(cellSize(ss, unit = "km"), ss)
            
            # Mask it if not global
            if (lat != "Global"){
                # Load area
                msk <- areas %>% filter(area == lat)
                
                ss_val <- trim(mask(ss_val, msk))
            }
            
            ss_val <- values(ss_val) %>% na.omit() %>% as.data.frame()
            ssp126 <- ss_val %>% select(area, SSP126)
            ssp370 <- ss_val %>% select(area, SSP370)
            ssp585 <- ss_val %>% select(area, SSP585)
            ss_1q <- c(
                weighted.quantile(ssp126$SSP126, ssp126$area, 0.25),
                weighted.quantile(ssp370$SSP370, ssp370$area, 0.25),
                weighted.quantile(ssp585$SSP585, ssp585$area, 0.25))
            ss_3q <- c(
                weighted.quantile(ssp126$SSP126, ssp126$area, 0.75),
                weighted.quantile(ssp370$SSP370, ssp370$area, 0.75),
                weighted.quantile(ssp585$SSP585, ssp585$area, 0.75))
            ss_medians <- c(
                weighted.quantile(ssp126$SSP126, ssp126$area, 0.5),
                weighted.quantile(ssp370$SSP370, ssp370$area, 0.5),
                weighted.quantile(ssp585$SSP585, ssp585$area, 0.5))
            
            uss_val <- c(cellSize(uss, unit = "km"), uss)
            
            # Mask it if not global
            if (lat != "Global"){
                # Load area
                msk <- areas %>% filter(area == lat)
                
                uss_val <- trim(mask(uss_val, msk))
            }
            
            uss_val <- values(uss_val) %>% na.omit() %>% as.data.frame()
            ssp126 <- uss_val %>% select(area, SSP126)
            ssp370 <- uss_val %>% select(area, SSP370)
            ssp585 <- uss_val %>% select(area, SSP585)
            uss_1q <- c(
                weighted.quantile(ssp126$SSP126, ssp126$area, 0.25),
                weighted.quantile(ssp370$SSP370, ssp370$area, 0.25),
                weighted.quantile(ssp585$SSP585, ssp585$area, 0.25))
            uss_3q <- c(
                weighted.quantile(ssp126$SSP126, ssp126$area, 0.75),
                weighted.quantile(ssp370$SSP370, ssp370$area, 0.75),
                weighted.quantile(ssp585$SSP585, ssp585$area, 0.75))
            uss_medians <- c(
                weighted.quantile(ssp126$SSP126, ssp126$area, 0.5),
                weighted.quantile(ssp370$SSP370, ssp370$area, 0.5),
                weighted.quantile(ssp585$SSP585, ssp585$area, 0.5))
            
            data.frame(driver = driver,
                       favorable_median = ss_medians,
                       favorable_1q = ss_1q,
                       favorable_3q = ss_3q,
                       unfavorable_median = uss_medians,
                       unfavorable_1q = uss_1q,
                       unfavorable_3q = uss_3q) %>%
                mutate(scenario = c("SSP126", "SSP370", "SSP585"))
        })) %>% mutate(time_period = time_period)
    }) %>% bind_rows() %>% mutate(area = lat)
}) %>% bind_rows()

write.csv(
    magnitude_periods, 
    file.path(tbl_dir, "supporting_table_magnitude_shift.csv"), 
    row.names = FALSE)

##### Supplementary Table 5 ####

magnitude_periods_tosave <- magnitude_periods %>% 
    pivot_longer(2:7, names_to = c("type", ".value"),
                 names_pattern = "(favorable|unfavorable)_(median|1q|3q)") %>% 
    mutate(driver = case_when(
        driver == "forest" ~ "FOR",
        driver == "human_impact" ~ "HLU",
        driver == "grassland" ~ "GRA",
        str_detect(driver, "bio") ~ toupper(driver)),
        type = str_to_title(type),
        area = ifelse(area == "Global", area, paste(area, "altitude"))) %>% 
    select(driver, area, scenario, time_period, type, median, `1q`, `3q`) %>%
    rename("Baseline" = type,
        "Variable" = driver, "Region" = area, "Scenario" = scenario,
        "Time period" = time_period, "SHAP value change (median)" = median,
        "SHAP value change (Q1)" = `1q`, "SHAP value change (Q3)" = `3q`)

write.csv(magnitude_periods_tosave, file.path(tbl_dir, "supplementary_table5.csv"), 
          row.names = FALSE)
rm(magnitude_periods_tosave)

#### Figure 4 and Extended Data Fig.5 ####

# Parameter setting
ssp <- "ssp370"
tp <- "2041-2070"

shift_pts <- magnitude_periods %>% filter(area == "Global") %>% 
    filter(scenario == toupper(ssp) & time_period == tp) %>% 
    mutate(driver = case_when(
        driver == "forest" ~ "FOR",
        driver == "human_impact" ~ "HLU",
        driver == "grassland" ~ "GRA",
        str_detect(driver, "bio") ~ toupper(driver)))

drivers <- shift_pts %>% 
    mutate(s_rank = rank(favorable_median),
           us_rank = rank(-unfavorable_median)) %>%
    mutate(SqRank = (s_rank^2) + (us_rank^2)/2) %>%
    mutate(RankOrder = rank(SqRank)) %>%
    arrange(-RankOrder) %>% pull(driver)

shift_pts <- shift_pts %>%
    select(driver, favorable_median, favorable_1q, favorable_3q,
           unfavorable_median, unfavorable_1q, unfavorable_3q) %>%
    pivot_longer(2:7, names_to = c("type", ".value"),
                 names_pattern = "(favorable|unfavorable)_(median|1q|3q)") %>%
    mutate(type = factor(
        type, levels = c("favorable", "unfavorable"),
        labels = c("Baseline-favorable area", "Baseline-unfavorable area"))) %>%
    mutate(driver = factor(driver, levels = drivers, labels = drivers))

# variable groups
grps <- data.frame(var = paste0("BIO", c(1, 5, 9, 3, 4, 2, 8, 7)),
                   group = "Temperature", color = "#e63946") %>% 
    rbind(data.frame(var = paste0("BIO", c(14, 12, 15, 18, 19)),
                     group = "Precipitation", color = "#1d3557")) %>% 
    rbind(data.frame(var = c("FOR", "HLU", "GRA"),
                     group = "Land cover", color = "#fb8500")) %>% 
    mutate(group = factor(group, levels = c("Temperature", "Precipitation", "Land cover")))

# Make all figures
figs <- lapply(c("SSP126", "SSP370", "SSP585"), function(ssp){
    figs <- lapply(time_periods, function(tp){
        if (ssp == "SSP370" & tp == "2041-2070"){
            titles <- c("(b) Magnitude shift", "(c) Shift by region")
        } else{
            titles <- c(
                sprintf("(%s %s)\nMagnitude shift", tolower(ssp), tp), 
                sprintf("(%s %s)\nShift by region", tolower(ssp), tp))}
        
        shift_pts <- magnitude_periods %>% filter(area == "Global") %>% 
            filter(scenario == toupper(ssp) & time_period == tp) %>% 
            mutate(driver = case_when(
                driver == "forest" ~ "FOR",
                driver == "human_impact" ~ "HLU",
                driver == "grassland" ~ "GRA",
                str_detect(driver, "bio") ~ toupper(driver)))
        
        drivers <- shift_pts %>% 
            mutate(s_rank = rank(favorable_median),
                   us_rank = rank(-unfavorable_median)) %>%
            mutate(SqRank = (s_rank^2) + (us_rank^2)/2) %>%
            mutate(RankOrder = rank(SqRank)) %>%
            arrange(-RankOrder) %>% pull(driver)
        
        shift_pts <- shift_pts %>%
            select(driver, favorable_median, favorable_1q, favorable_3q,
                   unfavorable_median, unfavorable_1q, unfavorable_3q) %>%
            pivot_longer(2:7, names_to = c("type", ".value"),
                         names_pattern = "(favorable|unfavorable)_(median|1q|3q)") %>%
            mutate(type = factor(
                type, levels = c("favorable", "unfavorable"),
                labels = c("Baseline-favorable area", "Baseline-unfavorable area"))) %>%
            mutate(driver = factor(driver, levels = drivers, labels = drivers))
        
        a_cols <- grps %>% arrange(match(var, drivers)) %>% pull(color)
        
        # Plot
        g1 <- ggplot(data = shift_pts) +
            geom_hline(yintercept = 0, color = "white") +
            geom_hline(yintercept = 0, color = "#e74c3c", linetype = "dotted") +
            geom_errorbar(
                data = shift_pts,
                aes(x = driver, ymin = `1q`, ymax = `3q`,
                    color = type, group = type),
                width = 0, linewidth = 1.5, alpha = 0.5,
                position = position_dodge(width = 1)) +
            geom_point(
                data = shift_pts,
                aes(y = median, x = driver, shape = type, group = type),
                color ="white", size = 2, position = position_dodge(width = 1)) +
            geom_point(
                data = shift_pts,
                aes(y = median, x = driver, shape = type,
                    color = type, group = type),
                size = 1.6, position = position_dodge(width = 1)) +
            scale_color_manual("", values = c('#a6611a', "#018571")) +
            scale_shape_manual(name = "", values = c(20, 18)) +
            coord_flip() + ggtitle(titles[1]) +
            labs(y = "\u2206 SHAP",
                 x = element_blank()) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            theme_pubclean(base_size = 11, base_family = 'Merriweather') +
            theme(axis.text = element_text(
                color = "black"),
                axis.text.y = element_text(color = a_cols),
                legend.text = element_text(size = 8),
                panel.grid.major.y = element_blank(),
                panel.grid.major.x = element_line(
                    linetype = "dotted", color = "lightgrey"),
                legend.position = "top",
                legend.key.height = unit(1, "mm"),
                legend.key.width = unit(1, "mm"),
                legend.spacing = unit(rep(0), "cm"),
                plot.margin = unit(c(0, 0.5, 0.2, 0.5), "cm"),
                plot.title = element_text(
                    family = "Merriweather", size = 11,
                    face = "bold", hjust = 0.5))
        
        # By regions
        shift_region_pts <- magnitude_periods %>% filter(area != "Global") %>% 
            filter(scenario == toupper(ssp) & time_period == tp) %>% 
            mutate(driver = case_when(
                driver == "forest" ~ "FOR",
                driver == "human_impact" ~ "HLU",
                driver == "grassland" ~ "GRA",
                str_detect(driver, "bio") ~ toupper(driver))) %>% 
            mutate(driver = factor(
                driver, levels = drivers, labels = drivers)) %>%
            mutate(area_p = paste0(area, "_p"), area_n = area) %>% 
            mutate(area_n = factor(
                area_n, levels = c("Low", "Middle", "High"))) %>% 
            mutate(area_p = factor(
                area_p, levels = c("Low_p", "Middle_p", "High_p")))
        
        fill_max <- max(shift_region_pts$favorable_median, 
                        shift_region_pts$unfavorable_median)
        fill_min <- min(shift_region_pts$favorable_median, 
                        shift_region_pts$unfavorable_median)
        
        g2 <- ggplot(shift_region_pts) +
            geom_tile(aes(x = area_p, y = driver, fill = favorable_median),
                      color = "white") + 
            scale_fill_continuous_diverging(
                palette = "Blue-Red 3",
                name = "Baseline\nfavorable area", 
                guide = guide_colorbar(order = 1), 
                limits = c(fill_min, fill_max)) +
            new_scale_fill() +
            geom_tile(aes(x = area_n, y = driver, fill = unfavorable_median),
                      color = "white") + 
            scale_fill_continuous_diverging(
                palette = "Purple-Green",
                name = "Baseline\nunfavorable area", 
                guide = guide_colorbar(order = 2),
                limits = c(fill_min, fill_max)) +
            scale_x_discrete(labels = rep(c("L", "M", "H"), 2)) +
            geom_vline(xintercept = 3.5, color = 'black', linewidth = 1) +
            labs(x = "Region", y = "") + coord_equal() + ggtitle(titles[2]) +
            theme_pubclean(base_family = "Merriweather", base_size = 11) +
            theme(panel.grid.major.y = element_line(color = "white"),
                  axis.text = element_text(color = "black"),
                  axis.text.y = element_text(color = a_cols),
                  plot.title = element_text(
                      family = "Merriweather", size = 11,
                      face = "bold", hjust = 0.5),
                  legend.position = "right", 
                  legend.direction = "vertical",
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 8),
                  legend.key.height = unit(3, "mm"),
                  legend.key.width = unit(2, "mm"),
                  plot.margin = unit(c(0, 0, 0.2, 0), "cm"))
        
        list("shift" = g1, "region" = g2)
    })
    
    # Add name
    names(figs) <- time_periods
    figs
}); names(figs) <- ssps

##### Figure 4 ####
img <- image_read(file.path(fig_dir, "fig4_flow.jpg"))
g <- image_ggplot(img, interpolate = TRUE)

ggarrange(g,
          ggarrange(NULL, as_ggplot(ggpubr::get_legend(figs[[2]][[2]][[1]])), 
                    NULL, nrow = 1, widths = c(1, 10, 10)),
          ggarrange(figs[[2]][[2]][[1]] + 
                        annotate(
                            geom = "text", y = 0.04, 
                            x = 7, label = "Temperature", color = "#e63946", 
                            family = "Merriweather", size = 3) +
                        annotate(
                            geom = "text", y = 0.04, x = 5.5, 
                            label = "Precipitation", color = "#1d3557", 
                            family = "Merriweather", size = 3) +
                        annotate(
                            geom = "text", y = 0.04, x = 4, 
                            label = "Landcover", color = "#fb8500", 
                            family = "Merriweather", size = 3) +
                        theme(legend.position = "none", 
                              plot.title = element_blank(),
                              plot.margin = unit(c(0, 0.2, 0.2, 0.8), "cm")), 
                    figs[[2]][[2]][[2]] +
                        theme(plot.title = element_blank()), ncol = 2, 
                    labels = c("(b)", "(c)"),
                    font.label = list(
                        size = 11, color = "black", 
                        face = "bold", family = "Merriweather")), 
          nrow = 3, heights = c(4.2, 1, 10),
          labels = c("(a)", ""),
          font.label = list(
              size = 11, color = "black", 
              face = "bold", family = "Merriweather"))

ggsave(file.path(fig_dir, "Figure4_shifts.png"), 
       width = 6.5, height = 4.4, dpi = 500, bg = "white")

##### Extended Data Fig.5 ####

figs <- lapply(ssps, function(ssp){
    figs <- lapply(time_periods, function(time_period){
        
        if (ssp == "ssp370" & time_period == "2041-2070"){
            texts <- ggplot() + 
                geom_point(aes(x = rep(1, 3), y = 1:3), color = "white",
                           show.legend = FALSE) + 
                annotate(
                geom = "text", y = 2.5, 
                x = 1, label = "Temperature", color = "#e63946", 
                family = "Merriweather", size = 5) +
                annotate(
                    geom = "text", y = 2, x = 1, 
                    label = "Precipitation", color = "#1d3557", 
                    family = "Merriweather", size = 5) +
                annotate(
                    geom = "text", y = 1.5, x = 1, 
                    label = "Landcover", color = "#fb8500", 
                    family = "Merriweather", size = 5) +
                theme_void()
            
            ggarrange(as_ggplot(ggpubr::get_legend(
                figs[[ssp]][[time_period]][[1]] +
                    theme(legend.position = "left",
                          legend.text = element_text(size = 13),
                          legend.title = element_text(size = 13)) +
                    guides(
                        color = guide_legend(
                            override.aes = list(size = 4, linewidth = 2))))),
                texts, nrow = 2)
        } else {
            g1 <- figs[[ssp]][[time_period]][[1]] +
                theme(legend.position = "none")
            g2 <- figs[[ssp]][[time_period]][[2]]
            
            ggarrange(g1, g2, ncol = 2)
        }
    })
    
    names(figs) <- time_periods
    figs
}); names(figs) <- ssps

figs <- do.call(c, figs)

ggarrange(plotlist = figs, nrow = 3, ncol = 3)

ggsave(file.path(fig_dir, "extended_data_fig5.png"), 
       width = 18, height = 11, dpi = 500, bg = "white")

# Clean up
rm(list = ls()); gc()
