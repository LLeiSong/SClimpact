# Load settings
source("R/figures/setting.R")

#### Figure 3 and Extended Data Fig.4 ####
# Load intensity patterns
species_periods <- read.csv(
    file.path(tbl_dir, "supporting_table_turnover_local.csv"))

figs_regs <- lapply(c("SSP126", "SSP370", "SSP585"), function(ssp){
    figs <- lapply(time_periods, function(tp){
        if (ssp == "SSP370"){
            if (tp == "2041-2070"){
                titles <- "(a) Intensity by region"
            } else if(tp == "2011-2040"){
                titles <- sprintf(
                    "a. Intensity by region (%s %s)", tolower(ssp), tp)
            } else{
                titles <- sprintf(
                    "f. Intensity by region (%s %s)", tolower(ssp), tp)
            }
            
        } else{
            titles <- sprintf("Intensity by region (%s %s)", tolower(ssp), tp)}
        
        # Species
        species_pts_all <- species_periods %>% 
            filter(scenario == ssp & time_period == tp) %>% 
            mutate(driver = case_when(
                driver == "forest" ~ "FOR",
                driver == "human_impact" ~ "HLU",
                driver == "grassland" ~ "GRA",
                str_detect(driver, "bio") ~ toupper(driver)))
        
        # Get the driver rank for visualization
        drivers <- species_pts_all %>% filter(area == "Global") %>% 
            mutate(stou_rank = rank(-stou_sp_median), 
                   utos_rank = rank(-utos_sp_median)) %>% 
            mutate(SqRank = (stou_rank^2) + (utos_rank^2) / 2) %>% 
            mutate(RankOrder = rank(SqRank)) %>% 
            arrange(RankOrder)
        
        # By regions
        species_regions_pts <- species_pts_all %>% filter(area != "Global") %>% 
            mutate(driver = factor(
                driver, levels = drivers$driver, labels = drivers$driver)) %>%
            mutate(area_p2n = paste0(area, "_p2n"), area_n2p = area) %>% 
            mutate(area_n2p = factor(
                area_n2p, levels = c("Low", "Middle", "High"))) %>% 
            mutate(area_p2n = factor(
                area_p2n, levels = c("Low_p2n", "Middle_p2n", "High_p2n")))
        
        fill_max <- max(species_regions_pts$stou_sp_median, 
                        species_regions_pts$utos_sp_median)
        
        ggplot(species_regions_pts) +
            geom_tile(aes(x = driver, y = area_n2p, fill = stou_sp_median),
                      color = "white") + 
            scale_fill_bs5(name = "Favoring (%)", "green",
                           guide = guide_colorbar(order = 2),
                           limits = c(0, fill_max)) +
            new_scale_fill() +
            geom_tile(aes(x = driver, y = area_p2n, fill = utos_sp_median),
                      color = "white") + 
            scale_fill_bs5(name = "Disfavoring (%)", "orange",
                           guide = guide_colorbar(order = 1),
                           limits = c(0, fill_max)) +
            scale_y_discrete(labels = rep(c("L", "M", "H"), 2)) +
            geom_hline(yintercept = 3.5, color = 'black', linewidth = 1) +
            labs(x = "", y = "") + ggtitle(titles) + 
            theme_pubclean(base_size = 12) +
            theme(panel.grid.major.y = element_line(color = "white"),
                  axis.text = element_text(color = "black", size = 8),
                  legend.position = "right", 
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 8),
                  legend.key.height = unit(2.3, "mm"),
                  legend.key.width = unit(1.8, "mm"),
                  legend.margin = margin(0, 0, -0.1, 0, unit = "cm"),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  plot.title = element_text(
                      size = 12, face = "bold", hjust = 0.5),
                  plot.margin = unit(c(0.24, 0.2, -0.2, 0), "cm"))
    })
    
    # Add name
    names(figs) <- time_periods
    figs
}); names(figs_regs) <- ssps

# Here we care about bio1 and bio5. bio12, and bio14.
# Use two as a group, otherwise RAM will explode.
driver_to_plots <- list(c("bio1", "bio5"), c("bio12", "bio14"))
ssp <- "ssp370"
for (time_period in time_periods){
    if (time_period == "2041-2070"){
        # This chunk is not very flexible, which only fit for this case.
        titles <- letters[2:5]
        titles_list <- split(titles, ceiling(seq_along(titles) / 2))
        titles_list <- lapply(1:2, function(x){
            titles <- titles_list[[x]]
            titles <- split(titles, ceiling(seq_along(titles) / 1))
            names(titles) <- driver_to_plots[[x]]
            titles})
    } else if (time_period == "2011-2040"){
        # This chunk is not very flexible, which only fit for this case.
        titles <- letters[2:5]
        titles_list <- split(titles, ceiling(seq_along(titles) / 2))
        titles_list <- lapply(1:2, function(x){
            titles <- titles_list[[x]]
            titles <- split(titles, ceiling(seq_along(titles) / 1))
            names(titles) <- driver_to_plots[[x]]
            titles})
    } else{
        # This chunk is not very flexible, which only fit for this case.
        titles <- letters[7:10]
        titles_list <- split(titles, ceiling(seq_along(titles) / 2))
        titles_list <- lapply(1:2, function(x){
            titles <- titles_list[[x]]
            titles <- split(titles, ceiling(seq_along(titles) / 1))
            names(titles) <- driver_to_plots[[x]]
            titles})
    }
    
    figs <- lapply(1:length(driver_to_plots), function(i){
        driver_to_plot <- driver_to_plots[[i]]
        
        if (!is.null(titles_list)){
            titles <- titles_list[[i]]
        } else titles <- NULL
        
        figs <- lapply(driver_to_plot, function(driver){
            # files
            fnames <- list.files(
                cc_dir, pattern = sprintf("%s.tif", time_period), 
                full.names = TRUE)
            fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
            fnames <- fnames[str_detect(fnames, sprintf("_%s_", ssp))]
            
            dir_fname <- fnames[str_detect(fnames, "dir_change")]
            num_fname <- fnames[str_detect(fnames, "num_species")]
            
            # S to U
            lyr_s2u <- project(
                rast(dir_fname)[[3]] / rast(num_fname)[[1]] * 100, plot_crs)
            # U to S
            lyr_u2s <- project(
                rast(dir_fname)[[2]] / rast(num_fname)[[1]] * 100, plot_crs)
            lyrs <- trim(c(lyr_s2u, lyr_u2s))
            
            # Put values into data.frame
            values_periods <- as.data.frame(lyrs, xy = TRUE) %>% 
                mutate(time_period = time_period)
            names(values_periods) <- c("x", "y", "P2N", "N2P", "time_period")
            
            ## Convert values to biviariates
            set.seed(123)
            values_bivi <- bi_class(
                values_periods, x = P2N, y = N2P, 
                style = "fisher", dim = 3, dig_lab = 3)
            set.seed(123)
            breaks <- bi_class_breaks(
                values_periods, x = P2N, y = N2P, style = "fisher", 
                dim = 3, dig_lab = 3, split = TRUE)
            breaks$bi_x <- round(breaks$bi_x, 0)
            breaks$bi_y <- round(breaks$bi_y, 0)
            
            lgd <- ggplotGrob(
                bi_legend(pal = "PurpleOr",
                          pad_color = "white",
                          flip_axes = FALSE,
                          rotate_pal = FALSE,
                          dim = 3,
                          xlab = "Disfavoring (%)",
                          ylab = "Favoring (%)",
                          breaks = breaks,
                          pad_width = 0.3,
                          size = 6) +
                    theme(axis.text = element_text(
                        size = 6, color = "black"),
                        axis.title = element_text(
                            size = 6, color = "black"),
                        panel.background = element_rect(fill='transparent'),
                        plot.background = element_rect(fill='transparent', color = NA),
                        axis.title.x = element_text(
                            margin = margin(t = 0.5, unit = "mm")),
                        axis.title.y = element_text(
                            margin = margin(r = -1, unit = "mm"))))
            
            values <- values_bivi
            rm(values_bivi, values_periods); gc()
            
            # Make the figure
            g <- ggplot() +
                geom_sf(data = bry, fill = "#d3d3d3", color = "#d3d3d3") + 
                geom_tile(
                    data = values, 
                    mapping = aes(x = x, y = y, fill = bi_class), 
                    show.legend = FALSE) +
                coord_sf(crs = plot_crs) +
                bi_scale_fill(
                    pal = "PurpleOr", dim = 3, 
                    flip_axes = FALSE, rotate_pal = FALSE) +
                annotation_custom(grob = lgd, 
                                  xmin = -18708181, xmax = -6808181,
                                  ymin = -8642702, ymax = 3357298) +
                theme_void() + 
                ggtitle(ifelse(
                    time_period != "2041-2070", 
                    sprintf("%s. %s (%s, %s)", 
                            titles[[driver]][1], toupper(driver), ssp, time_period), 
                    sprintf("(%s) %s", titles[[driver]][1], toupper(driver)))) +
                theme(plot.margin = unit(c(-0.2, -0.3, 0, -0.3), "cm"),
                      plot.title = element_text(size = 12,
                          face = "bold", hjust = 0.5, 
                          margin = margin(b = -5)))
            
            rm(values); gc()
            g
        })
        
        ggarrange(plotlist = figs, ncol = 2)
    })
    
    ggarrange(figs_regs[["ssp370"]][[time_period]],
              ggarrange(plotlist = figs, nrow = 2),
              nrow = 2, heights = c(1.6, 3.4))
    
    if (time_period == "2041-2070"){
        fname <- file.path(fig_dir, "figure4_spatial.png")
    } else {
        fname <- file.path(
            fig_dir, sprintf("extended_data_fig4_%s.png", time_period))
    }
    
    ggsave(fname, width = 6.5, height = 4.9, dpi = 500, bg = "white")
    rm(figs); gc()
}

##### Extended Data Fig.4 ####
fnames <- list.files(fig_dir, pattern = "extended_data_fig4_",
                     full.names = TRUE)
figs <- do.call(c, lapply(fnames, image_read))
fig <- image_append(figs, stack = TRUE)
image_write(fig, file.path(fig_dir, "extended_data_fig4.png"))

# Clean up
fnames <- list.files(
    fig_dir, pattern = "extended_data_fig4", full.names = TRUE)
fnames <- fnames[!str_detect(fnames, "extended_data_fig4.png")]
file.remove(fnames)

#### Supplementary Figures 5-13 ####
driver_to_plots <- c(unlist(driver_to_plots),
    paste0("bio", c(2:4, 7:9, 15, 18, 19)), 
           "forest", "grassland", "human_impact")
driver_to_plots <- split(driver_to_plots, 
                         ceiling(seq_along(driver_to_plots) / 4))

for (time_period in time_periods){
    for (ssp in ssps){
        for (i in 1:length(driver_to_plots)){
            driver_to_plot <- driver_to_plots[[i]]
            
            figs <- lapply(driver_to_plot, function(driver){
                # Get name
                title_nm <- case_when(
                    driver == "forest" ~ "FOR",
                    driver == "human_impact" ~ "HLU",
                    driver == "grassland" ~ "GRA",
                    str_detect(driver, "bio") ~ toupper(driver))
                
                # files
                fnames <- list.files(
                    cc_dir, pattern = sprintf("%s.tif", time_period), 
                    full.names = TRUE)
                fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
                fnames <- fnames[str_detect(fnames, sprintf("_%s_", ssp))]
                
                dir_fname <- fnames[str_detect(fnames, "dir_change")]
                num_fname <- fnames[str_detect(fnames, "num_species")]
                
                # S to U
                lyr_s2u <- project(
                    rast(dir_fname)[[3]] / rast(num_fname)[[1]] * 100, plot_crs)
                # U to S
                lyr_u2s <- project(
                    rast(dir_fname)[[2]] / rast(num_fname)[[1]] * 100, plot_crs)
                lyrs <- trim(c(lyr_s2u, lyr_u2s))
                
                # Put values into data.frame
                values_periods <- as.data.frame(lyrs, xy = TRUE) %>% 
                    mutate(time_period = time_period)
                names(values_periods) <- c("x", "y", "P2N", "N2P", "time_period")
                
                ## Convert values to biviariates
                set.seed(123)
                values_bivi <- bi_class(
                    values_periods, x = P2N, y = N2P, 
                    style = "fisher", dim = 3, dig_lab = 3)
                set.seed(123)
                breaks <- bi_class_breaks(
                    values_periods, x = P2N, y = N2P, style = "fisher", 
                    dim = 3, dig_lab = 3, split = TRUE)
                breaks$bi_x <- round(breaks$bi_x, 0)
                breaks$bi_y <- round(breaks$bi_y, 0)
                
                lgd <- ggplotGrob(
                    bi_legend(pal = "PurpleOr",
                              pad_color = "white",
                              flip_axes = FALSE,
                              rotate_pal = FALSE,
                              dim = 3,
                              xlab = "Disfavoring (%)",
                              ylab = "Favoring (%)",
                              breaks = breaks,
                              pad_width = 0.3,
                              size = 6) +
                        theme(axis.text = element_text(
                            size = 6, color = "black"),
                            axis.title = element_text(
                                size = 6, color = "black"),
                            panel.background = element_rect(fill='transparent'),
                            plot.background = element_rect(fill='transparent', color = NA),
                            axis.title.x = element_text(
                                margin = margin(t = 0.5, unit = "mm")),
                            axis.title.y = element_text(
                                margin = margin(r = -1, unit = "mm"))))
                
                values <- values_bivi
                rm(values_bivi, values_periods); gc()
                
                # Make the figure
                g <- ggplot() +
                    geom_sf(data = bry, fill = "#d3d3d3", color = "#d3d3d3") + 
                    geom_tile(
                        data = values, 
                        mapping = aes(x = x, y = y, fill = bi_class), 
                        show.legend = FALSE) +
                    coord_sf(crs = plot_crs) +
                    bi_scale_fill(
                        pal = "PurpleOr", dim = 3, 
                        flip_axes = FALSE, rotate_pal = FALSE) +
                    annotation_custom(grob = lgd, 
                                      xmin = -18708181, xmax = -6808181,
                                      ymin = -8642702, ymax = 3357298) +
                    theme_void() + 
                    ggtitle(sprintf("%s (%s %s)", title_nm, ssp, time_period)) +
                    theme(plot.margin = unit(c(-0.2, -0.3, 0, -0.3), "cm"),
                          plot.title = element_text(size = 12,
                                                    face = "bold", hjust = 0.5, 
                                                    margin = margin(b = -5)))
                
                rm(values); gc()
                g
            })
            
            if (i == 1){
                ggarrange(figs_regs[[ssp]][[time_period]],
                          ggarrange(plotlist = figs, ncol = 2, nrow = 2),
                          nrow = 2, heights = c(1.6, 3.4))
                
                fname <- file.path(
                    fig_dir, sprintf("Figure_s_sp_%s_%s_%s.png", i - 1, ssp, time_period))
                
                ggsave(fname, width = 6.5, height = 5, dpi = 500, bg = "white")
            } else {
                ggarrange(plotlist = figs, ncol = 2, nrow = 2)
                
                fname <- file.path(
                    fig_dir, sprintf("Figure_s_sp_%s_%s_%s.png", i - 1, ssp, time_period))
                
                ggsave(fname, width = 6.5, height = 3.4, dpi = 500, bg = "white")
            }
            
            rm(figs); gc()
        }
    }
}

# Clean up
rm(list = ls()); gc()
