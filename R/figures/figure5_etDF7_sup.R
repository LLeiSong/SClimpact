# Load settings
source("R/figures/setting.R")

#### Figure 5 ####
## Set parameters
drivers <- c("bio1", "forest")
ssp <- "ssp370"
s <- 1 # 1 or 2

for (time_period in time_periods){
    figs <- lapply(drivers, function(driver){
        nm <- case_when(
            driver == "forest" ~ "FOR",
            driver == "human_impact" ~ "HLU",
            driver == "grassland" ~ "GRA",
            str_detect(driver, "bio") ~ toupper(driver))
        
        sus <- ifelse(s == 1, "s", "us")
        lb <- ifelse(s == 1, "baseline-favorable", "baseline-unfavorable")
        
        if (time_period == "2041-2070"){
            plot_title <- sprintf(
                "(%s) %s over %s area", letters[1:2][drivers == driver], nm, lb)
        } else {
            plot_title <- sprintf("%s (over %s area %s)", nm, lb, time_period)
        }
        
        # Load layers
        # Suitable area
        shap_changes <- do.call(c, lapply(time_periods, function(time_period){
            fnames <- list.files(
                cc_dir, pattern = sprintf("%s_%s.tif", ssp, time_period), 
                full.names = TRUE)
            fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
            
            val_fnames <- fnames[str_detect(fnames, sprintf("val_change_%s", driver))]
            num_fnames <- fnames[str_detect(fnames, "num_species")]
            
            ss <- rast(val_fnames)[[s]] / rast(num_fnames)[[s + 1]]
            names(ss) <- ssp
            ss
        })); names(shap_changes) <- time_periods
        
        val_current <- rast("data/variables/Env/AllEnv.tif", lyrs = driver)
        val_changes <- do.call(c, lapply(time_periods, function(time_period){
            temp_future <- rast(
                sprintf("data/variables/OtherEnvMean/%s_%s.tif", 
                        ssp, time_period), lyrs = driver)
            temp_future - val_current
        })); names(val_changes) <- time_periods
        
        num_species <- do.call(c, lapply(time_periods, function(time_period){
            rast(sprintf("results/climate_change/num_species_%s_%s_%s.tif", 
                         driver, ssp, time_period))[[1]]
        })); names(num_species) <- time_periods
        
        ## Convert values to bivariates
        ## Calculate statistics
        data <- c(cellSize(shap_changes[[time_period]], unit = "km"), 
                  val_changes[[time_period]], shap_changes[[time_period]],
                  num_species[[time_period]])
        data <- as.data.frame(data, xy = TRUE) %>% na.omit()
        names(data) <- c("x", "y", "area", "variable", "SHAP", "num_species")
        
        # Remove non-zero values
        data <- data %>% filter(variable != 0) %>% filter(SHAP != 0)
        
        data$SHAP_br <- ifelse(data$SHAP < 0, 1, 2)
        data$SHAP_br <- factor(data$SHAP_br, levels = c(1, 2))
        data$variable_br <- ifelse(data$variable <= 0, 1, 2)
        data$variable_br <- factor(data$variable_br, levels = c(1, 2))
        data <- bi_class(data, x = variable_br, y = SHAP_br, dim = 2)
        total_area <- sum(data$area)
        statistics <- data %>% 
            dplyr::group_by(bi_class) %>% 
            summarise(variable_1q = weighted.quantile(variable, area, 0.25),
                      variable_3q = weighted.quantile(variable, area, 0.75),
                      variable_median = weighted.quantile(variable, area, 0.5),
                      SHAP_1q = weighted.quantile(SHAP, area, 0.25),
                      SHAP_3q = weighted.quantile(SHAP, area, 0.75),
                      SHAP_median = weighted.quantile(SHAP, area, 0.5),
                      area = sum(area) / total_area * 100) %>% 
            filter(area >= 0.01) # Otherwise it will show 0 in the table
        statistics$area <- sprintf("%.2f%%", statistics$area)
        
        digit_var <- ifelse(driver == "bio1", 0, 1)
        digit_shap <- 2
        
        statistics$SHAP <- sprintf(
            "%.2f (%.2f-%.2f)", statistics$SHAP_median * 10^digit_shap, 
            statistics$SHAP_1q * 10^digit_shap, statistics$SHAP_3q * 10^digit_shap)
        statistics$variable <- sprintf(
            "%.2f (%.2f-%.2f)", statistics$variable_median * 10^digit_var, 
            statistics$variable_1q * 10^digit_var, 
            statistics$variable_3q * 10^digit_var)
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
            ifelse(digit_var == 0, sprintf("%s change", nm), 
                   ifelse(digit_var > 0, sprintf("%s~change~(10^-%s)", nm, digit_var),
                          sprintf("%s~change~(10^%s)", nm, abs(digit_var)))), 
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
                  val_changes[[time_period]], shap_changes[[time_period]],
                  num_species[[time_period]])
        data <- project(data, plot_crs)
        data <- as.data.frame(data, xy = TRUE) %>% na.omit()
        names(data) <- c("x", "y", "area", "variable", "SHAP", "num_species")
        
        # Remove non-zero values
        data <- data %>% filter(variable != 0) %>% filter(SHAP != 0)
        data <- data %>% select(x, y, area, num_species, variable, SHAP)
        
        if (all(data$SHAP < 0)) {
            data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, 1)
            data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, 1)
        }
        if (all(data$SHAP > 0)) {
            data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, -1)
            data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, -1)
        }
        if (all(data$variable < 0)) {
            data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, 1)
            data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, -1)
        }
        if (all(data$variable > 0)) {
            data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, 1)
            data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, -1)
        }
        
        data$SHAP <- ifelse(data$SHAP < 0, 1, 2)
        data$SHAP <- as.factor(data$SHAP)
        
        data$variable <- ifelse(data$variable < 0, 1, 2)
        data$variable <- as.factor(data$variable)
        
        # species
        nums <- quantile(data$num_species, c(0.33, 0.66, 1.0), na.rm = TRUE)
        data$num_species <- ifelse(data$num_species <= nums[1], 1,
                                   ifelse(data$num_species <= nums[2], 2, 3))
        data$num_species <- as.factor(data$num_species)
        
        data_classified <- bi_class(data, x = variable, y = SHAP, dim = 2)
        breaks <- bi_class_breaks(data, x = variable, y = SHAP, 
                                  dim = 2, dig_lab = 2, split = TRUE)
        breaks$bi_x <- c("Drop", "Rise")
        breaks$bi_y <- c("Drop", "Rise")
        
        leg <- biscale:::bi_pal_pull(
            pal = "DkBlue2", dim = 2, 
            flip_axes = FALSE, 
            rotate_pal = FALSE)
        leg <- data.frame(
            bi_class = names(leg),
            bi_fill = leg)
        
        leg <- lapply(c(0.7, 0.8, 1.0), function(val){
            leg %>% mutate(value = val) %>% mutate(group = val)
        }) %>% bind_rows() %>% 
            mutate(group = as.factor(group)) %>% 
            mutate(bi_class = factor(
                bi_class, levels = c("2-2", "2-1", "1-1","1-2"),
                labels = c("2-2", "2-1", "1-1","1-2")))
        
        pol_lgd <- ggplotGrob(
            ggplot(leg,
                   aes(x = bi_class,
                       y = value, alpha = group,
                       fill = bi_class)) +
                geom_col(width = 1, color = "white", linewidth = 0.2) +
                geom_text(x = 0, y = 0, label = 0, 
                          size = 2, family = "Merriweather") +
                geom_text(x = 0.5, y = 1, label = round(nums[1], 0), 
                          size = 2, family = "Merriweather") +
                geom_text(x = 0.5, y = 1.8, label = round(nums[2], 0), 
                          size = 2, family = "Merriweather") +
                geom_text(x = 0.5, y = 2.5, label = round(nums[3], 0), 
                          size = 2, family = "Merriweather") +
                scale_alpha_manual("", values = c(1.0, 0.7, 0.4)) +
                ggtitle("Species richness") +
                scale_fill_manual(
                    "", values = c("#434e87", "#52b6b6", "#d3d3d3", "#ad5b9c")) +
                coord_polar() + theme_void() +
                theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_blank(),
                      plot.title = element_text(
                          hjust = 0.5, family = "Merriweather",
                          size = 8, margin = margin(0, 0, -16, 0))))
        
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
                pad_width = 0,
                size = 10) + 
                scale_fill_manual(values = rep("transparent", 4)) +
                annotation_custom(grob = pol_lgd,
                                  xmin = 0.4, xmax = 2.8,
                                  ymin = 0.4, ymax = 2.8) +
                theme(axis.text = element_text(
                    size = 8, color = "black", family = 'Merriweather'),
                    axis.title = element_text(
                        size = 8, color = "black", family = 'Merriweather'),
                    panel.background = element_rect(fill='transparent'),
                    plot.background = element_rect(fill='transparent', 
                                                   color = NA),
                    axis.title.x = element_text(
                        margin = margin(t = -15, unit = "mm")),
                    axis.title.y = element_text(
                        margin = margin(r = -45, unit = "mm")),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_blank(),
                    legend.position = "none"))
        
        # Create the bivariate map using ggplot2
        ggplot() +
            geom_sf(data = bry, fill = "transparent", 
                    color = "#d3d3d3", linewidth = 0.1) + 
            geom_tile(data = data_classified, 
                      mapping = aes(x = x, y = y, fill = bi_class, 
                                    alpha = num_species), 
                      show.legend = FALSE) +
            bi_scale_fill(
                pal = "DkBlue2", dim = 2,
                flip_axes = FALSE, rotate_pal = FALSE) +
            scale_alpha_manual("", values = c(0.4, 0.7, 1.0)) +
            ggtitle(plot_title) +
            annotation_custom(grob = lgd, 
                              xmin = -18108181, xmax = -6908181,
                              ymin = -11142702, ymax = 1457298) +
            annotation_custom(grob = ggplotGrob(tbl), 
                              xmin = -10108181, xmax = 18208181,
                              ymin = -13542702, ymax = 107298) +
            theme_void() + 
            theme(plot.margin = unit(c(-1, -0.4, 0.1, -0.4), "cm"),
                  plot.title = element_text(
                      family = "Merriweather", size = 10,
                      face = "bold", hjust = 0.5,
                      margin = margin(0, 0, -5, 0)))
    })
    
    ggarrange(plotlist = figs, nrow = 2)
    
    name <- ifelse(
        time_period == "2041-2070", 
        file.path(fig_dir, "Figure5_spatial_shift.png"),
        file.path(fig_dir, sprintf("extended_data_fig7_%s.png", time_period)))
    
    ggsave(name, width = 6, height = 6.5, dpi = 500)
}

#### Extended Data Fig.7 ####
fnames <- file.path(
    fig_dir, c(sprintf("extended_data_fig7_%s.png", 
                       c("2011-2040", "2071-2100"))))

figs <- do.call(c, lapply(fnames, image_read))
fig <- image_append(figs)
image_write(
    fig, file.path(fig_dir, "extended_data_fig7.png"))

# Clean up
file.remove(fnames); rm(ssp, drivers, s, figs, fig, fnames)

#### Supplementary ####
## Set parameters
ss <- 1:2 # 1 or 2

for (driver in features){
    nm <- case_when(
        driver == "forest" ~ "FOR",
        driver == "human_impact" ~ "HLU",
        driver == "grassland" ~ "GRA",
        str_detect(driver, "bio") ~ toupper(driver))
    
    for (ssp in c("ssp126", "ssp370", "ssp585")){
        for (time_period in time_periods){
            figs <- lapply(ss, function(s){
                sus <- ifelse(s == 1, "s", "us")
                lb <- ifelse(s == 1, "baseline-favorable", 
                             "baseline-unfavorable")
                
                plot_title <- sprintf("%s (over %s area, %s %s)", 
                                      nm, lb, ssp, time_period)
                
                # Load layers
                # Suitable area
                shap_changes <- do.call(c, lapply(time_periods, function(time_period){
                    fnames <- list.files(
                        cc_dir, pattern = sprintf("%s_%s.tif", ssp, time_period), 
                        full.names = TRUE)
                    fnames <- fnames[str_detect(fnames, sprintf("_%s_", driver))]
                    
                    val_fnames <- fnames[str_detect(fnames, sprintf("val_change_%s", driver))]
                    num_fnames <- fnames[str_detect(fnames, "num_species")]
                    
                    ss <- rast(val_fnames)[[s]] / rast(num_fnames)[[s + 1]]
                    names(ss) <- ssp
                    ss
                })); names(shap_changes) <- time_periods
                
                val_current <- rast("data/variables/Env/AllEnv.tif", lyrs = driver)
                val_changes <- do.call(c, lapply(time_periods, function(time_period){
                    temp_future <- rast(
                        sprintf("data/variables/OtherEnvMean/%s_%s.tif", 
                                ssp, time_period), lyrs = driver)
                    temp_future - val_current
                })); names(val_changes) <- time_periods
                
                num_species <- do.call(c, lapply(time_periods, function(time_period){
                    rast(sprintf("results/climate_change/num_species_%s_%s_%s.tif", 
                                 driver, ssp, time_period))[[1]]
                })); names(num_species) <- time_periods
                
                ## Convert values to bivariates
                ## Calculate statistics
                data <- c(cellSize(shap_changes[[time_period]], unit = "km"), 
                          val_changes[[time_period]], shap_changes[[time_period]],
                          num_species[[time_period]])
                data <- as.data.frame(data, xy = TRUE) %>% na.omit()
                names(data) <- c("x", "y", "area", "variable", "SHAP", "num_species")
                
                # Remove non-zero values
                data <- data %>% filter(variable != 0) %>% filter(SHAP != 0)
                
                data$SHAP_br <- ifelse(data$SHAP < 0, 1, 2)
                data$SHAP_br <- factor(data$SHAP_br, levels = c(1, 2))
                data$variable_br <- ifelse(data$variable < 0, 1, 2)
                data$variable_br <- factor(data$variable_br, levels = c(1, 2))
                data <- bi_class(data, x = variable_br, y = SHAP_br, dim = 2)
                total_area <- sum(data$area)
                statistics <- data %>% 
                    dplyr::group_by(bi_class) %>% 
                    summarise(variable_1q = weighted.quantile(variable, area, 0.25),
                              variable_3q = weighted.quantile(variable, area, 0.75),
                              variable_median = weighted.quantile(variable, area, 0.5),
                              SHAP_1q = weighted.quantile(SHAP, area, 0.25),
                              SHAP_3q = weighted.quantile(SHAP, area, 0.75),
                              SHAP_median = weighted.quantile(SHAP, area, 0.5),
                              area = sum(area) / total_area * 100) %>% 
                    filter(area >= 0.01)
                statistics$area <- sprintf("%.2f%%", statistics$area)
                
                digit_var <- log10(1 / max(abs(statistics$variable_median)))
                digit_var <- ceiling(digit_var)
                digit_var <- ifelse(
                    max(abs(statistics$variable_median)) * 10^digit_var >= 10, 
                    digit_var - 1, digit_var)
                
                digit_shap <- log10(1 / max(abs(statistics$SHAP_median)))
                digit_shap <- ceiling(digit_shap)
                digit_shap <- ifelse(
                    max(abs(statistics$SHAP_median)) * 10^digit_shap >= 10, 
                    digit_shap - 1, digit_shap)
                
                statistics$SHAP <- sprintf(
                    "%.2f (%.2f-%.2f)", statistics$SHAP_median * 10^digit_shap, 
                    statistics$SHAP_1q * 10^digit_shap, statistics$SHAP_3q * 10^digit_shap)
                statistics$variable <- sprintf(
                    "%.2f (%.2f-%.2f)", statistics$variable_median * 10^digit_var, 
                    statistics$variable_1q * 10^digit_var, 
                    statistics$variable_3q * 10^digit_var)
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
                    ifelse(digit_var == 0, sprintf("%s change", nm), 
                           ifelse(digit_var > 0, 
                                  sprintf("%s~change~(10^-%s)", nm, digit_var),
                                  sprintf("%s~change~(10^%s)", nm, abs(digit_var)))), 
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
                          val_changes[[time_period]], shap_changes[[time_period]],
                          num_species[[time_period]])
                data <- project(data, plot_crs)
                data <- as.data.frame(data, xy = TRUE) %>% na.omit()
                names(data) <- c("x", "y", "area", "variable", "SHAP", "num_species")
                
                # Remove non-zero values
                data <- data %>% filter(variable != 0) %>% filter(SHAP != 0)
                data <- data %>% select(x, y, area, num_species, variable, SHAP)
                
                if (all(data$SHAP < 0)) {
                    data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, 1)
                    data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, 1)
                }
                if (all(data$SHAP > 0)) {
                    data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, -1)
                    data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, -1)
                }
                if (all(data$variable < 0)) {
                    data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, 1)
                    data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, -1)
                }
                if (all(data$variable > 0)) {
                    data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, 1)
                    data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, -1)
                }
                
                data$SHAP <- ifelse(data$SHAP < 0, 1, 2)
                data$SHAP <- as.factor(data$SHAP)
                
                data$variable <- ifelse(data$variable < 0, 1, 2)
                data$variable <- as.factor(data$variable)
                
                # species
                nums <- quantile(data$num_species, c(0.33, 0.66, 1.0), na.rm = TRUE)
                data$num_species <- ifelse(data$num_species <= nums[1], 1,
                                           ifelse(data$num_species <= nums[2], 2, 3))
                data$num_species <- as.factor(data$num_species)
                
                data_classified <- bi_class(data, x = variable, y = SHAP, dim = 2)
                breaks <- bi_class_breaks(data, x = variable, y = SHAP, 
                                          dim = 2, dig_lab = 2, split = TRUE)
                breaks$bi_x <- c("Drop", "Rise")
                breaks$bi_y <- c("Drop", "Rise")
                
                leg <- biscale:::bi_pal_pull(
                    pal = "DkBlue2", dim = 2, 
                    flip_axes = FALSE, 
                    rotate_pal = FALSE)
                leg <- data.frame(
                    bi_class = names(leg),
                    bi_fill = leg)
                
                leg <- lapply(c(0.7, 0.8, 1.0), function(val){
                    leg %>% mutate(value = val) %>% mutate(group = val)
                }) %>% bind_rows() %>% 
                    mutate(group = as.factor(group)) %>% 
                    mutate(bi_class = factor(
                        bi_class, levels = c("2-2", "2-1", "1-1","1-2"),
                        labels = c("2-2", "2-1", "1-1","1-2")))
                
                pol_lgd <- ggplotGrob(
                    ggplot(leg,
                           aes(x = bi_class,
                               y = value, alpha = group,
                               fill = bi_class)) +
                        geom_col(width = 1, color = "white", linewidth = 0.2) +
                        geom_text(x = 0, y = 0, label = 0, 
                                  size = 2, family = "Merriweather") +
                        geom_text(x = 0.5, y = 1, label = round(nums[1], 0), 
                                  size = 2, family = "Merriweather") +
                        geom_text(x = 0.5, y = 1.8, label = round(nums[2], 0), 
                                  size = 2, family = "Merriweather") +
                        geom_text(x = 0.5, y = 2.5, label = round(nums[3], 0), 
                                  size = 2, family = "Merriweather") +
                        scale_alpha_manual("", values = c(1.0, 0.7, 0.4)) +
                        ggtitle("Species richness") +
                        scale_fill_manual(
                            "", values = c("#434e87", "#52b6b6", "#d3d3d3", "#ad5b9c")) +
                        coord_polar() + theme_void() +
                        theme(legend.position = "none",
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              axis.line = element_blank(),
                              plot.title = element_text(
                                  hjust = 0.5, family = "Merriweather",
                                  size = 8, margin = margin(0, 0, 0, 0))))
                
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
                        pad_width = 0,
                        size = 10) + 
                        scale_fill_manual(values = rep("transparent", 4)) +
                        annotation_custom(grob = pol_lgd,
                                          xmin = 0.4, xmax = 2.8,
                                          ymin = 0.4, ymax = 2.8) +
                        theme(axis.text = element_text(
                            size = 8, color = "black", family = 'Merriweather'),
                            axis.title = element_text(
                                size = 8, color = "black", family = 'Merriweather'),
                            panel.background = element_rect(fill='transparent'),
                            plot.background = element_rect(fill='transparent', 
                                                           color = NA),
                            axis.title.x = element_text(
                                margin = margin(t = -10, unit = "mm")),
                            axis.title.y = element_text(
                                margin = margin(r = -40, unit = "mm")),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            axis.line = element_blank(),
                            legend.position = "none"))
                
                # Create the bivariate map using ggplot2
                ggplot() +
                    geom_sf(data = bry, fill = "transparent", 
                            color = "#d3d3d3", linewidth = 0.1) + 
                    geom_tile(data = data_classified, 
                              mapping = aes(x = x, y = y, fill = bi_class, 
                                            alpha = num_species), 
                              show.legend = FALSE) +
                    bi_scale_fill(
                        pal = "DkBlue2", dim = 2,
                        flip_axes = FALSE, rotate_pal = FALSE) +
                    scale_alpha_manual("", values = c(0.4, 0.7, 1.0)) +
                    ggtitle(plot_title) +
                    annotation_custom(grob = lgd, 
                                      xmin = -18108181, xmax = -6908181,
                                      ymin = -11142702, ymax = 1457298) +
                    annotation_custom(grob = ggplotGrob(tbl), 
                                      xmin = -10108181, xmax = 18208181,
                                      ymin = -13542702, ymax = 107298) +
                    theme_void() + 
                    theme(plot.margin = unit(c(-1, -0.4, 0.1, -0.4), "cm"),
                          plot.title = element_text(
                              family = "Merriweather", size = 10,
                              face = "bold", hjust = 0.5,
                              margin = margin(0, 0, -5, 0)))
            })
            
            ggarrange(plotlist = figs, ncol = 2)
            
            name <- file.path(
                fig_dir, sprintf("Figure_s_shift_%s_%s_%s.png", 
                                 nm, ssp, time_period))
            
            ggsave(name, width = 12, height = 3.3, dpi = 500)
        }
    }
}

# Clean up
rm(list = ls()); gc()

