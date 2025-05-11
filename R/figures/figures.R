#### Model evaluation ####
species_list <- read.csv(
    file.path(root_dir, "data/occurrences", "species_qualified_sdm.csv")) %>% 
    arrange(num_occ)
species_list <- species_list$species

evals <- lapply(species_list, function(sp){
    dr <- file.path(sdm_dir, sp)
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

write.csv(evals, file.path(tbl_dir, "model_evaluation.csv"), row.names = FALSE)

smr_evals <- evals %>% filter(type == "testing") %>% 
    group_by(metric.eval) %>% 
    summarise(mean = mean(value), sd = sd(value)) %>% 
    ungroup() %>% mutate(value = sprintf("%.2f\u00B1%.2f", mean, sd)) %>% 
    arrange(-mean) %>% 
    select(metric.eval, value)
colnames(smr_evals) <- NULL

tbl <- ggplotGrob(ggtexttable(
    smr_evals, rows = NULL, 
    theme = ttheme(
        padding = unit(c(1, 2), "mm"),
        colnames.style = colnames_style(
            size = 8, fill = "transparent", parse = TRUE, 
            face = "plain", font = "Merriweather"),
        tbody.style = tbody_style(
            size = 8, color = "black", font = "Merriweather"))))

g <- ggplot(data = evals %>% filter(type == "testing"), 
            aes(x = value, after_stat(density), fill = metric.eval)) +
    geom_density(alpha = 0.8) +
    xlab("Metric value(0 - 1)") +
    ylab("Density") +
    scale_fill_brewer(name = "Metric", palette = "Dark2") +
    theme_pubclean(base_size = 11, base_family = "Merriweather") +
    theme(axis.text = element_text(color = "black", family = 'Merriweather'))

ggarrange(g, tbl, nrow = 1, widths = c(0.7, 0.3))

ggsave(file.path(fig_dir, "Figure_s4_model_eval.png"),
       width = 5, height = 3, dpi = 500, bg = "white")

#### Variable selection ####
# Subset the species
species_list <- unique(evals$species)

# Variables
vars <- data.frame(
    var = names(rast(file.path(root_dir, "data/variables/Env", "AllEnv.tif"))),
    num = 0)
for (sp in species_list) {
    dr <- file.path(root_dir, "data/variables/variable_list")
    vars_sl <- read.csv(file.path(dr, sprintf("%s.csv", sp)))
    vars_sl <- na.omit(vars_sl$var_uncorrelated)
    vars[vars$var %in% vars_sl, 'num'] <- 
        vars[vars$var %in% vars_sl, 'num'] + 1
}

vars <- vars %>% mutate(ratio = num / length(species_list)) %>% 
    mutate(var = ifelse(var == "forest", "Forest coverage",
                        ifelse(var == "human_impact", "Human land use",
                               gsub("bio", "BIO", var)))) %>% 
    mutate(var = ifelse(var == "grassland", "Grassland", var)) %>% 
    mutate(selected = ifelse(ratio > 0.1, "yes", "no")) %>% 
    mutate(var = fct_reorder(var, num))

write.csv(vars, file.path(tbl_dir, "variable_selection.csv"), row.names = FALSE)

ggplot(data = vars, 
       aes(x = num, xend = 0,
           y = var, yend = var,
           color = selected)) +
    geom_segment() +
    geom_point() +
    scale_color_manual("", values = c("grey", "black")) +
    xlab("No. of species") + ylab("") +
    theme_pubclean(base_size = 11, base_family = "Merriweather") +
    theme(axis.text.x = element_text(
        color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(color = "black", size = 11),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(vjust = -1),
        legend.position = "none")

ggsave(file.path(fig_dir, "Figure_s2_vars_selected.png"),
       width = 4, height = 4.2, dpi = 500, bg = "white")

#### SHAP: number of Monte Carlo iterations ####
shaps <- lapply(species_list, function(sp){
    read.csv(
        file.path(sdm_dir, sp, sprintf("shap_cor_%s.csv", sp))) %>% 
        select(nshap, cor) %>% mutate(species = sp)
}) %>% bind_rows() %>% filter(nshap > 10)

shaps_mean <- shaps %>% group_by(nshap) %>% 
    summarise(sd = sd(cor), cor = mean(cor))

ggplot(shaps_mean) +
    geom_line(data = shaps, aes(x = nshap, y = cor, group = species), 
              color = "lightgrey", linewidth = 0.2) +
    geom_line(aes(x = nshap, y = cor)) +
    geom_vline(xintercept = 1000, color = "darkgrey", linetype = "dashed") +
    geom_point(aes(x = nshap, y = cor), color = "#EB5B00") +
    geom_errorbar(aes(x = nshap, ymin = cor - sd, ymax = cor + sd),
                  color = "#EB5B00", width = 300) +
    xlab("Number of repetitions (nsim)") + 
    ylab("Correlation with\nSHAP values(sim = 10,000)") +
    theme_pubclean(base_size = 11, base_family = "Merriweather") +
    theme(axis.text.x = element_text(
        color = "black"),
        axis.text.y = element_text(color = "black"))

ggsave(file.path(fig_dir, "Figure_s3_nshap.png"),
       width = 3, height = 4, dpi = 500, bg = "white")

#### SHAP magnitude shifts ####
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

for (ssp in c("SSP126", "SSP370", "SSP585")){
    # Manipulate data for plotting
    pts <- features_periods %>% 
        filter(scenario == ssp) %>% 
        mutate(driver = ifelse(driver == "forest", "Forest coverage",
                               ifelse(driver == "human_impact", "Human land use",
                                      gsub("bio", "BIO", driver)))) %>% 
        mutate(driver = ifelse(driver == "grassland", "Grassland coverage", driver))
    
    drivers <- pts %>% filter(time_period == "2011-2040") %>% 
        arrange(-suitable_change) %>% pull(driver)
    
    pts <- pts %>% 
        select(driver, suitable_change, unsuitable_change, time_period) %>% 
        pivot_longer(2:3, names_to = "type", values_to = "mean") %>% 
        mutate(type = factor(
            type, levels = c("suitable_change", "unsuitable_change"),
            labels = c("Baseline-favorable area", "Baseline-unfavorable area"))) %>% 
        mutate(driver = factor(driver, levels = drivers, labels = drivers)) %>% 
        mutate(driver = as.integer(driver)) %>% 
        left_join(
            pts %>% 
                select(driver, suitable_sd, unsuitable_sd, time_period) %>% 
                pivot_longer(2:3, names_to = "type", values_to = "sd") %>% 
                mutate(type = factor(
                    type, levels = c("suitable_sd", "unsuitable_sd"),
                    labels = c("Baseline-favorable area", "Baseline-unfavorable area"))) %>% 
                mutate(driver = factor(driver, levels = drivers, labels = drivers)) %>% 
                mutate(driver = as.integer(driver)),
            by = c("driver", "time_period", "type"))
    
    segments_pos <- pts %>% 
        filter(type == "Baseline-favorable area") %>% select(-sd) %>% 
        pivot_wider(names_from = time_period, values_from = mean) %>% 
        mutate(range1 = `2041-2070` - `2011-2040`,
               range2 = `2071-2100` - `2041-2070`)
    
    segments_pos_ins_1 <- segments_pos %>% 
        mutate(`2041-2070` = ifelse(range1 >= -0.0006, NA, `2041-2070`)) %>% 
        select(driver, type, `2011-2040`, `2041-2070`) %>% na.omit()
    
    segments_pos_ins_2 <- segments_pos %>% 
        mutate(`2071-2100` = ifelse(range2 >= -0.0006, NA, `2071-2100`)) %>% 
        select(driver, type, `2041-2070`, `2071-2100`) %>% na.omit()
    
    segments_pos_dec_1 <- segments_pos %>% 
        mutate(`2041-2070` = ifelse(range1 <= 0.0006, NA, `2041-2070`)) %>% 
        select(driver, type, `2011-2040`, `2041-2070`) %>% na.omit()
    
    segments_pos_dec_2 <- segments_pos %>% 
        mutate(`2071-2100` = ifelse(range2 <= 0.0006, NA, `2071-2100`)) %>% 
        select(driver, type, `2041-2070`, `2071-2100`) %>% na.omit()
    
    segments_neg <- pts %>% 
        filter(type == "Baseline-unfavorable area") %>% select(-sd) %>% 
        pivot_wider(names_from = time_period, values_from = mean) %>% 
        mutate(range1 = `2041-2070` - `2011-2040`,
               range2 = `2071-2100` - `2041-2070`)
    
    segments_neg_ins_1 <- segments_neg %>% 
        mutate(`2041-2070` = ifelse(range1 <= 0.0006, NA, `2041-2070`)) %>% 
        select(driver, type, `2011-2040`, `2041-2070`) %>% na.omit()
    
    segments_neg_ins_2 <- segments_neg %>% 
        mutate(`2071-2100` = ifelse(range2 <= 0.0006, NA, `2071-2100`)) %>% 
        select(driver, type, `2041-2070`, `2071-2100`) %>% na.omit()
    
    segments_neg_dec_1 <- segments_neg %>% 
        mutate(`2041-2070` = ifelse(range1 >= -0.0006, NA, `2041-2070`)) %>% 
        select(driver, type, `2011-2040`, `2041-2070`) %>% na.omit()
    
    segments_neg_dec_2 <- segments_neg %>% 
        mutate(`2071-2100` = ifelse(range2 >= -0.0006, NA, `2071-2100`)) %>% 
        select(driver, type, `2041-2070`, `2071-2100`) %>% na.omit()
    
    # Plot
    ggplot() +
        geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
        geom_errorbar(data = pts,
                      aes(x = driver, ymin = mean - sd, ymax = mean + sd,
                          color = time_period, group = time_period),
                      width = 0, linewidth = 2.8, alpha = 0.5) +
        geom_point(data = pts,
                   aes(y = mean, x = driver, 
                       group = time_period, shape = type), color ="white",
                   size = 3.4) +
        geom_point(data = pts,
                   aes(y = mean, x = driver, color = time_period,
                       group = time_period, shape = type), 
                   size = 2.6) +
        # add curves for suitable area
        geom_curve(
            data = segments_pos_ins_1,
            aes(x = driver + 0.2, y = `2011-2040`, 
                xend = driver + 0.21, yend = `2041-2070`, 
                group = type), col = '#757B82', 
            arrow = arrow(length = unit(0.01, "npc")), curvature = 0.2) +
        geom_curve(
            data = segments_pos_ins_2,
            aes(x = driver + 0.2, y = `2041-2070`, 
                xend = driver + 0.21, yend = `2071-2100`, 
                group = type), col = '#757B82', 
            arrow = arrow(length = unit(0.01, "npc")), curvature = 0.2) +
        geom_curve(
            data = segments_pos_dec_1,
            aes(x = driver - 0.2, y = `2011-2040`, 
                xend = driver - 0.21, yend = `2041-2070`, 
                group = type), col = '#757B82', 
            arrow = arrow(length = unit(0.01, "npc")), curvature = 0.2) +
        geom_curve(
            data = segments_pos_dec_2,
            aes(x = driver - 0.2, y = `2041-2070`, 
                xend = driver - 0.21, yend = `2071-2100`, 
                group = type), col = '#757B82', 
            arrow = arrow(length = unit(0.01, "npc")), curvature = 0.2) +
        # add curves for unsuitable area
        geom_curve(
            data = segments_neg_ins_1,
            aes(x = driver + 0.2, y = `2011-2040`, 
                xend = driver + 0.21, yend = `2041-2070`, 
                group = type), col = '#757B82', 
            arrow = arrow(length = unit(0.01, "npc")), curvature = -0.2) +
        geom_curve(
            data = segments_neg_ins_2,
            aes(x = driver + 0.2, y = `2041-2070`, 
                xend = driver + 0.21, yend = `2071-2100`, 
                group = type), col = '#757B82', 
            arrow = arrow(length = unit(0.01, "npc")), curvature = -0.2) +
        geom_curve(
            data = segments_neg_dec_1,
            aes(x = driver - 0.2, y = `2011-2040`, 
                xend = driver - 0.21, yend = `2041-2070`, 
                group = type), col = '#757B82', 
            arrow = arrow(length = unit(0.01, "npc")), curvature = -0.2) +
        geom_curve(
            data = segments_neg_dec_2,
            aes(x = driver - 0.2, y = `2041-2070`, 
                xend = driver - 0.21, yend = `2071-2100`, 
                group = type), col = '#757B82', 
            arrow = arrow(length = unit(0.01, "npc")), curvature = -0.2) +
        coord_flip() +
        labs(y = "SHAP value change",
             x = element_blank()) +
        scale_color_brewer(name = "  Time period", palette = "Dark2") +
        scale_shape_manual(name = "", values = c(20, 18)) +
        scale_y_continuous(expand = expansion(mult = 0.01)) +
        scale_x_continuous(breaks = c(1:17), labels = c(drivers, "")) +
        theme_pubclean(base_size = 10, base_family = 'Merriweather') +
        theme(axis.text = element_text(
            color = "black", size = 10),
            legend.text = element_text(size = 10),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.justification = c(-0.5, 0),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            axis.title.x = element_text(vjust = -1, hjust = 0.32),
            plot.subtitle = element_text(size = 10, vjust = -25, hjust = 2.7),
            legend.box = 'vertical', legend.box.just = "center") +
        guides(color = guide_legend("  Time period"),
               fill = guide_legend("  Time period"))
    
    fname <- ifelse(ssp == "SSP370", "docs/figures/Figure5_magnitude_shift.png",
                    sprintf("docs/figures/Figure_s_%s_magnitude_shift.png", ssp))
    ggsave(fname, width = 6.5, height = 5, dpi = 500)
}

#### Spatial magnitude shifts ####
## Set parameters
driver <- "bio1"
ssp <- "ssp370"
s <- 2 # 1 or 2
nm <- ifelse(driver == "forest", "Forest coverage",
             ifelse(driver == "human_impact", "Human land use",
                    gsub("bio", "BIO", driver)))
nm <- ifelse(nm == "grassland", "Grassland", nm)
sus <- ifelse(s == 1, "s", "us")

for (time_period in time_periods){
    # Load layers
    # Suitable area
    shap_changes <- do.call(c, lapply(time_periods, function(time_period){
        fnames <- list.files(
            data_dir, pattern = sprintf("%s_%s.tif", ssp, time_period), 
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
              val_changes[[time_period]], shap_changes[[time_period]],
              num_species[[time_period]])
    data <- project(data, plot_crs)
    data <- as.data.frame(data, xy = TRUE) %>% na.omit()
    names(data) <- c("x", "y", "area", "variable", "SHAP", "num_species")
    
    data <- data %>% select(x, y, area, num_species, variable, SHAP)
    
    # data <- data %>% 
    #     filter(SHAP >= quantile(SHAP, 0.01) & 
    #                SHAP <= quantile(SHAP, 0.99))
    
    if (all(data$SHAP <= 0)) {
        data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, 1)
        data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, 1)
    }
    if (all(data$SHAP > 0)) {
        data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, -1)
        data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, -1)
    }
    if (all(data$variable <= 0)) {
        data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, 1)
        data[nrow(data) + 1, ] <- c(rep(NA, 4), 1, -1)
    }
    if (all(data$variable > 0)) {
        data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, 1)
        data[nrow(data) + 1, ] <- c(rep(NA, 4), -1, -1)
    }
    
    data$SHAP <- ifelse(data$SHAP <= 0, 1, 2)
    data$SHAP <- as.factor(data$SHAP)
    
    data$variable <- ifelse(data$variable <= 0, 1, 2)
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
                      size = 8, margin = margin(0, 0, -2, 0))))
    
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
                              xmin = 0.4, xmax = 2.6,
                              ymin = 0.4, ymax = 2.6) +
            theme(axis.text = element_text(
                size = 8, color = "black", family = 'Merriweather'),
                axis.title = element_text(
                    size = 8, color = "black", family = 'Merriweather'),
                panel.background = element_rect(fill='transparent'),
                plot.background = element_rect(fill='transparent', 
                                               color = NA),
                axis.title.x = element_text(margin = margin(t = -10, unit = "mm")),
                axis.title.y = element_text(margin = margin(r = -40, unit = "mm")),
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
                  mapping = aes(x = x, y = y, fill = bi_class, alpha = num_species), 
                  show.legend = FALSE) +
        bi_scale_fill(
            pal = "DkBlue2", dim = 2,
            flip_axes = FALSE, rotate_pal = FALSE) +
        scale_alpha_manual("", values = c(0.4, 0.7, 1.0)) +
        annotation_custom(grob = lgd, 
                          xmin = -18108181, xmax = -5908181,
                          ymin = -10142702, ymax = 1557298) +
        annotation_custom(grob = ggplotGrob(tbl), 
                          xmin = -10108181, xmax = 18208181,
                          ymin = -13542702, ymax = 107298) +
        theme_void() + 
        theme(plot.margin = unit(c(-2.5, -0.4, -1, -0.6), "cm"))
    
    name <- ifelse(
        time_period == "2041-2070", 
        sprintf("docs/figures/Figure67_%s_%s_shap_%s.png", 
                sus, driver, time_period),
        sprintf("docs/figures/Figure_s_%s_%s_shap_%s.png", 
                sus, driver, time_period))
    
    ggsave(name, width = 6, height = 3.1, dpi = 500)
}
