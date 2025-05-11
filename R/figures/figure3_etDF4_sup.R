# Load settings
source("R/figures/setting.R")

#### Figure 3 and Extended Data Fig.4 ####
# Here we care about bio1 and bio5. bio12, and bio14.
# Use two as a group, otherwise RAM will explode.
driver_to_plots <- list(c("bio1", "bio5"), c("bio12", "bio14"))

for (tp in time_periods){
    if (tp == "2041-2070"){
        # This chunk is not very flexible, which only fit for this case.
        titles <- letters[1:8]
        titles_list <- split(titles, ceiling(seq_along(titles) / 4))
        titles_list <- lapply(1:2, function(x){
            titles <- titles_list[[x]]
            titles <- split(titles, ceiling(seq_along(titles) / 2))
            names(titles) <- driver_to_plots[[x]]
            titles})
    } else titles_list <- NULL
    
    for (i in 1:length(driver_to_plots)){
        driver_to_plot <- driver_to_plots[[i]]
        
        if (!is.null(titles_list)){
            titles <- titles_list[[i]]
        } else titles <- NULL
        
        figs <- lapply(driver_to_plot, function(driver){
            # Figure of var from suitable to unsuitable
            # S to U
            values_s2u <- lapply(time_periods, function(time_period){
                fnames <- list.files(
                    cc_dir, pattern = sprintf("%s.tif", time_period), 
                    full.names = TRUE)
                fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
                
                dir_fnames <- fnames[str_detect(fnames, "dir_change")]
                num_fnames <- fnames[str_detect(fnames, "num_species")]
                
                # S to U
                lyrs <- do.call(c, lapply(1:length(dir_fnames), function(i){
                    project(rast(dir_fnames[i])[[3]] / rast(num_fnames[i])[[1]] * 100, 
                            plot_crs)
                })); names(lyrs) <- c("SSP126", "SSP370", "SSP580")
                lyrs <- trim(lyrs)
                
                means <- mean(lyrs, na.rm = TRUE)
                ranges <- range(lyrs, na.rm = TRUE)
                ranges <- ranges$range_max - ranges$range_min
                
                # Put values into data.frame
                data <- as.data.frame(c(means, ranges), xy = TRUE)
                names(data) <- c("x", "y", "Mean", "Range")
                
                data %>% mutate(time_period = time_period)
            })
            values_s2u <- values_s2u %>% bind_rows() %>% mutate(type = "P2N")
            
            # U to S
            values_u2s <- lapply(time_periods, function(time_period){
                fnames <- list.files(
                    cc_dir, pattern = sprintf("%s.tif", time_period), 
                    full.names = TRUE)
                fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
                
                dir_fnames <- fnames[str_detect(fnames, "dir_change")]
                num_fnames <- fnames[str_detect(fnames, "num_species")]
                
                # U to S
                lyrs <- do.call(c, lapply(1:length(dir_fnames), function(i){
                    project(rast(dir_fnames[i])[[2]] / rast(num_fnames[i])[[1]] * 100, 
                            plot_crs)
                })); names(lyrs) <- c("SSP126", "SSP370", "SSP580")
                lyrs <- trim(lyrs)
                
                means <- mean(lyrs, na.rm = TRUE)
                ranges <- range(lyrs, na.rm = TRUE)
                ranges <- ranges$range_max - ranges$range_min
                
                # Put values into data.frame
                data <- as.data.frame(c(means, ranges), xy = TRUE)
                names(data) <- c("x", "y", "Mean", "Range")
                
                data %>% mutate(time_period = time_period)
            })
            values_u2s <- values_u2s %>% bind_rows() %>% mutate(type = "N2P")
            
            values_periods <- rbind(values_s2u, values_u2s)
            rm(values_s2u, values_u2s)
            
            ## Convert values to biviariates
            values_bivi <- bi_class(
                values_periods, x = Mean, y = Range, 
                style = "fisher", dim = 4, dig_lab = 4)
            breaks <- bi_class_breaks(
                values_periods, x = Mean, y = Range, style = "fisher", 
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
                          ylab = "Range (%)",
                          breaks = breaks,
                          pad_width = 0.3,
                          size = 6) +
                    theme(axis.text = element_text(
                        size = 6, color = "black", family = 'Merriweather'),
                        axis.title = element_text(
                            size = 6, color = "black", family = 'Merriweather'),
                        panel.background = element_rect(fill='transparent'),
                        plot.background = element_rect(fill='transparent', color = NA),
                        axis.text.x = element_text(
                            angle = -45, vjust = 0.5, hjust = 0),
                        axis.title.x = element_text(
                            margin = margin(t = -26, unit = "mm")),
                        axis.title.y = element_text(
                            margin = margin(r = -26, unit = "mm"))))
            
            values <- values_bivi %>% filter(time_period == tp)
            rm(values_bivi, values_periods); gc()
            
            # Make the figure
            g_p2n <- ggplot() +
                geom_sf(data = bry, fill = "#d3d3d3", color = "#d3d3d3") + 
                geom_tile(
                    data = values %>% filter(type == "P2N"), 
                    mapping = aes(x = x, y = y, fill = bi_class), 
                    show.legend = FALSE) +
                coord_sf(crs = plot_crs) +
                bi_scale_fill(
                    pal = "PurpleOr", dim = 4, 
                    flip_axes = FALSE, rotate_pal = FALSE) +
                annotation_custom(grob = lgd, 
                                  xmin = -18708181, xmax = -6808181,
                                  ymin = -8642702, ymax = 3357298) +
                theme_void() + 
                ggtitle(ifelse(is.null(titles), 
                               sprintf("%s (P2N %s)", toupper(driver), tp), 
                               sprintf("(%s) %s-P2N", titles[[driver]][1], 
                                       toupper(driver)))) +
                theme(plot.margin = unit(c(-0.2, -0.3, 0, -0.3), "cm"),
                      plot.title = element_text(
                          family = "Merriweather", size = 10,
                          face = "bold", hjust = 0.5))
            
            g_n2p <- ggplot() +
                geom_sf(data = bry, fill = "#d3d3d3", color = "#d3d3d3") + 
                geom_tile(
                    data = values %>% filter(type == "N2P"), 
                    mapping = aes(x = x, y = y, fill = bi_class), 
                    show.legend = FALSE) +
                coord_sf(crs = plot_crs) +
                bi_scale_fill(
                    pal = "PurpleOr", dim = 4, 
                    flip_axes = FALSE, rotate_pal = FALSE) +
                annotation_custom(grob = lgd, 
                                  xmin = -18708181, xmax = -6808181,
                                  ymin = -8642702, ymax = 3357298) +
                theme_void() + 
                ggtitle(ifelse(is.null(titles), 
                               sprintf("%s (N2P %s)", toupper(driver), tp), 
                               sprintf("(%s) %s-N2P", titles[[driver]][2], 
                                       toupper(driver)))) +
                theme(plot.margin = unit(c(-0.2, -0.3, 0, -0.3), "cm"),
                      plot.title = element_text(
                          family = "Merriweather", size = 10,
                          face = "bold", hjust = 0.5))
            
            rm(values); gc()
            ggarrange(g_p2n, g_n2p, ncol = 2)
        })
        
        ggarrange(plotlist = figs, nrow = 2)
        
        if (tp == "2041-2070"){
            fname <- file.path(fig_dir, sprintf("figure3_spatial_%s.png", i))
        } else {
            fname <- file.path(
                fig_dir, sprintf("extended_data_fig4_%s_%s.png", i, tp))
        }
        
        ggsave(fname, width = 6.5, height = 3.4, dpi = 500)
        rm(figs); gc()
    }
}

##### Figure 3 ####
fnames <- file.path(fig_dir, sprintf("figure3_spatial_%s.png", 1:2))
figs <- do.call(c, lapply(fnames, image_read))
fig <- image_append(figs, stack = TRUE)
image_write(fig, file.path(fig_dir, "figure3_spatial.png"))

##### Extended Data Fig.4 ####
for (tp in c("2011-2040", "2071-2100")){
    fnames <- file.path(
        fig_dir, c(sprintf("extended_data_fig4_%s_%s.png", 1:2, tp)))
    
    figs <- do.call(c, lapply(fnames, image_read))
    fig <- image_append(figs, stack = TRUE)
    image_write(
        fig, file.path(fig_dir, sprintf("extended_data_fig4_%s.png", tp)))
}

# Append again
fnames <- file.path(
    fig_dir, sprintf("extended_data_fig4_%s.png", 
                     c("2011-2040", "2071-2100")))
figs <- do.call(c, lapply(fnames, image_read))
fig <- image_append(figs)
image_write(fig, file.path(fig_dir, "extended_data_fig4.png"))

#### Supplementary Figure ####
driver_to_plots <- c(
    paste0("bio", c(2:4, 7:9, 15, 18, 19)), 
           "forest", "grassland", "human_impact")
driver_to_plots <- split(driver_to_plots, 
                         ceiling(seq_along(driver_to_plots) / 2))

for (tp in time_periods){
    for (i in 1:length(driver_to_plots)){
        driver_to_plot <- driver_to_plots[[i]]
        
        figs <- lapply(driver_to_plot, function(driver){
            # Get name
            title_nm <- case_when(
                driver == "forest" ~ "FOR",
                driver == "human_impact" ~ "HLU",
                driver == "grassland" ~ "GRA",
                str_detect(driver, "bio") ~ toupper(driver))
            
            # Figure of var from suitable to unsuitable
            # S to U
            values_s2u <- lapply(time_periods, function(time_period){
                fnames <- list.files(
                    cc_dir, pattern = sprintf("%s.tif", time_period), 
                    full.names = TRUE)
                fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
                
                dir_fnames <- fnames[str_detect(fnames, "dir_change")]
                num_fnames <- fnames[str_detect(fnames, "num_species")]
                
                # S to U
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
            values_s2u <- values_s2u %>% bind_rows() %>% mutate(type = "P2N")
            
            # U to S
            values_u2s <- lapply(time_periods, function(time_period){
                fnames <- list.files(
                    cc_dir, pattern = sprintf("%s.tif", time_period), 
                    full.names = TRUE)
                fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
                
                dir_fnames <- fnames[str_detect(fnames, "dir_change")]
                num_fnames <- fnames[str_detect(fnames, "num_species")]
                
                # U to S
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
            values_u2s <- values_u2s %>% bind_rows() %>% mutate(type = "N2P")
            
            values_periods <- rbind(values_s2u, values_u2s); rm(values_s2u, values_u2s)
            
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
                          pad_width = 0.3,
                          size = 6) +
                    theme(axis.text = element_text(
                        size = 6, color = "black", family = 'Merriweather'),
                        axis.title = element_text(
                            size = 6, color = "black", family = 'Merriweather'),
                        panel.background = element_rect(fill='transparent'),
                        plot.background = element_rect(fill='transparent', color = NA),
                        axis.text.x = element_text(
                            angle = -45, vjust = 0.5, hjust = 0),
                        axis.title.x = element_text(
                            margin = margin(t = -23, unit = "mm")),
                        axis.title.y = element_text(
                            margin = margin(r = -16, unit = "mm"))))
            
            values <- values_bivi %>% filter(time_period == tp)
            rm(values_bivi, values_periods); gc()
            
            # Make the figure
            g_p2n <- ggplot() +
                geom_sf(data = bry, fill = "#d3d3d3", color = "#d3d3d3") + 
                geom_tile(
                    data = values %>% filter(type == "P2N"), 
                    mapping = aes(x = x, y = y, fill = bi_class), 
                    show.legend = FALSE) +
                coord_sf(crs = plot_crs) +
                bi_scale_fill(
                    pal = "PurpleOr", dim = 4, 
                    flip_axes = FALSE, rotate_pal = FALSE) +
                annotation_custom(grob = lgd, 
                                  xmin = -18708181, xmax = -6808181,
                                  ymin = -8642702, ymax = 3357298) +
                theme_void() + 
                ggtitle(sprintf("%s (P2N %s)", title_nm, tp)) +
                theme(plot.margin = unit(c(-0.2, -0.3, 0, -0.3), "cm"),
                      plot.title = element_text(
                          family = "Merriweather", size = 10,
                          face = "bold", hjust = 0.5))
            
            g_n2p <- ggplot() +
                geom_sf(data = bry, fill = "#d3d3d3", color = "#d3d3d3") + 
                geom_tile(
                    data = values %>% filter(type == "N2P"), 
                    mapping = aes(x = x, y = y, fill = bi_class), 
                    show.legend = FALSE) +
                coord_sf(crs = plot_crs) +
                bi_scale_fill(
                    pal = "PurpleOr", dim = 4, 
                    flip_axes = FALSE, rotate_pal = FALSE) +
                annotation_custom(grob = lgd, 
                                  xmin = -18708181, xmax = -6808181,
                                  ymin = -8642702, ymax = 3357298) +
                theme_void() + 
                ggtitle(sprintf("%s (N2P %s)", title_nm, tp)) +
                theme(plot.margin = unit(c(-0.2, -0.3, 0, -0.3), "cm"),
                      plot.title = element_text(
                          family = "Merriweather", size = 10,
                          face = "bold", hjust = 0.5))
            
            rm(values); gc()
            ggarrange(g_p2n, g_n2p, ncol = 2)
        })
        
        ggarrange(plotlist = figs, nrow = 2)
        
        fname <- file.path(
            fig_dir, sprintf("Figure_s_sp_%s_%s.png", i, tp))
        
        ggsave(fname, width = 6.5, height = 3.4, dpi = 500)
        rm(figs); gc()
    }
}
