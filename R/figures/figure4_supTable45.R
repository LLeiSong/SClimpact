# Load settings
source("R/figures/setting.R")

#### Analyze the magnitude shifts for each species ####
sp_analysis <- lapply(species_list, function(sp){
    fname <- file.path(sp_analysis_dir, sprintf("changes_%s.csv", sp))
    read.csv(fname) %>% 
        separate(scenario, c('scenario', 'year'), sep = "_") %>% 
        filter(class == "magnitude change") %>%
        select(-c(class)) %>% 
        pivot_wider(names_from = c(metrics), values_from = value) %>% 
        na.omit()
}) %>% bind_rows() %>% 
    mutate(sp = gsub("_", " ", sp))

# Attach IUCN status to the species analysis results
sp_analysis <- left_join(sp_analysis, sp_names, by = c("sp" = "species"))

# Standardize the names 
sp_analysis <- sp_analysis %>% filter(feature %in% features) %>%
    mutate(feature = case_when(
        feature == "forest" ~ "FOR",
        feature == "human_impact" ~ "HLU",
        feature == "grassland" ~ "GRA",
        str_detect(feature, "bio") ~ toupper(feature)))

##### Supplementary Table 4 ####

sp_analysis_tosave <- sp_analysis %>%
    pivot_wider(names_from = type, values_from = c(mean, sd)) %>% 
    select(scenario, year, sp, feature, mean_P, sd_P, 
           mean_N, sd_N, category, status) %>% 
    rename("Species" = sp, "IUCN category" = category, "Variable" = feature, 
           "Year" = year, "Scenario" = scenario, "Status" = status,
           "Baseline-suitable area\n(mean SHAP value change)" = mean_P,
           "Baseline-suitable area\n(SD of SHAP value change)" = sd_P,
           "Baseline-unsuitable area\n(mean SHAP value change)" = mean_N, 
           "Baseline-unsuitable area\n(SD of SHAP value change)" = sd_N)
write.csv(
    sp_analysis_tosave, file.path(tbl_dir, "supplementary_table4.csv"),
    row.names = FALSE, na = "")
rm(sp_analysis_tosave)

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
                rast(val_fnames[i])[[1]] / rast(num_fnames[i])[[2]]
            })); names(ss) <- c("SSP126", "SSP370", "SSP580")
            ss <- trim(ss)
            
            # Unsuitable
            uss <- do.call(c, lapply(1:length(val_fnames), function(i){
                rast(val_fnames[i])[[2]] / rast(num_fnames[i])[[3]]
            })); names(uss) <- c("SSP126", "SSP370", "SSP580")
            uss <- trim(uss)
            
            ss_val <- c(cellSize(ss, unit = "km"), ss)
            
            # Mask it if not global
            if (lat != "Global"){
                # Load area
                msk <- areas %>% filter(area == lat)
                
                ss_val <- trim(mask(ss_val, msk))
            }
            
            ss_val <- values(ss_val) %>% na.omit() %>% as.data.frame()
            
            ssp126 <- ss_val %>% select(area, SSP126) %>%
                filter(SSP126 >= quantile(SSP126, 0.01) & 
                           SSP126 <= quantile(SSP126, 0.99))
            ssp370 <- ss_val %>% select(area, SSP370) %>%
                filter(SSP370 >= quantile(SSP370, 0.01) & 
                           SSP370 <= quantile(SSP370, 0.99))
            ssp580 <- ss_val %>% select(area, SSP580) %>%
                filter(SSP580 >= quantile(SSP580, 0.01) & 
                           SSP580 <= quantile(SSP580, 0.99))
            
            ss_sds <- c(wtd.var(ssp126$SSP126, ssp126$area),
                        wtd.var(ssp370$SSP370, ssp370$area),
                        wtd.var(ssp580$SSP580, ssp580$area))
            
            ss_means <- c(wtd.mean(ssp126$SSP126, ssp126$area),
                          wtd.mean(ssp370$SSP370, ssp370$area),
                          wtd.mean(ssp580$SSP580, ssp580$area))
            
            uss_val <- c(cellSize(uss, unit = "km"), uss)
            
            # Mask it if not global
            if (lat != "Global"){
                # Load area
                msk <- areas %>% filter(area == lat)
                
                uss_val <- trim(mask(uss_val, msk))
            }
            
            uss_val <- values(uss_val) %>% na.omit() %>% as.data.frame()
            
            ssp126 <- uss_val %>% select(area, SSP126) %>%
                filter(SSP126 >= quantile(SSP126, 0.01) & 
                           SSP126 <= quantile(SSP126, 0.99))
            ssp370 <- uss_val %>% select(area, SSP370) %>%
                filter(SSP370 >= quantile(SSP370, 0.01) & 
                           SSP370 <= quantile(SSP370, 0.99))
            ssp580 <- uss_val %>% select(area, SSP580) %>%
                filter(SSP580 >= quantile(SSP580, 0.01) & 
                           SSP580 <= quantile(SSP580, 0.99))
            
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
    }) %>% bind_rows() %>% mutate(area = lat)
}) %>% bind_rows()

write.csv(
    magnitude_periods, 
    file.path(tbl_dir, "supporting_table_magnitude_shift.csv"), 
    row.names = FALSE)

##### Supplementary Table 5 ####

magnitude_periods_tosave <- magnitude_periods %>% 
    select(driver, area, scenario, time_period, suitable_change, 
           suitable_sd, unsuitable_change, unsuitable_sd) %>%
    mutate(driver = case_when(
        driver == "forest" ~ "FOR",
        driver == "human_impact" ~ "HLU",
        driver == "grassland" ~ "GRA",
        str_detect(driver, "bio") ~ toupper(driver))) %>% 
    mutate(area = ifelse(
        area == "Global", area, paste(area, "altitude"))) %>% as_tibble() %>% 
    rename("Variable" = driver, "Region" = area, "Scenario" = scenario,
           "Time period" = time_period, 
           "Baseline-suitable area\n(mean SHAP value change)" = suitable_change,
           "Baseline-suitable area\n(SD of SHAP value change)" = suitable_sd,
           "Baseline-unsuitable area\n(mean SHAP value change)" = unsuitable_change, 
           "Baseline-unsuitable area\n(SD of SHAP value change)" = unsuitable_sd)
write.csv(magnitude_periods_tosave, file.path(tbl_dir, "supplementary_table5.csv"), 
          row.names = FALSE)
rm(magnitude_periods_tosave)

#### Figure 4 ####
# Parameter setting
ssp <- "ssp370"
time_period <- "2041-2070"

# panel c
sp_analysis_fig <- sp_analysis %>% 
    filter(scenario == ssp & year == time_period)

# Add a common order based on the order of ALL status
feature_orders <- sp_analysis_fig %>% 
    group_by(feature, type, scenario, year) %>%
    summarise(mean = median(mean), .groups = "drop") %>%
    mutate(plot_order = paste(feature, type, sep = "_")) %>% 
    arrange(mean) %>% pull(plot_order)

sp_analysis_fig <- sp_analysis_fig %>% 
    mutate(plot_order = paste(feature, type, sep = "_")) %>% 
    mutate(plot_order = factor(plot_order, levels = feature_orders,
                               labels = feature_orders))

# Check if the distribution of EN group is different from all
mw_test <- lapply(unique(sp_analysis_fig$feature), function(driver){
    lapply(unique(sp_analysis_fig$type), function(tp){
        dat <- sp_analysis_fig %>% filter(feature == driver & type == tp)
        dat_en <- dat %>% filter(status == "EN")
        wt <- wilcox.test(dat$mean, dat_en$mean)
        data.frame(type = tp,
                   feature = driver,
                   p.value = wt$p.value)
    }) %>% bind_rows()
}) %>% bind_rows() %>% 
    mutate(significance = ifelse(
        p.value < 0.05, "Significant",
        ifelse(p.value >= 0.05 & p.value < 0.1, 
               "Marginally significant",
               "Not significant"))) %>% 
    mutate(significance = factor(
        significance,
        levels = c("Significant", "Marginally significant", 
                   "Not significant"),
        labels = c("p < 0.05", " 0.05 \U2264 P < 0.1", 
                   "P \u2265 0.1")))

en_fig <- sp_analysis_fig %>% 
    filter(status == "EN") %>% 
    group_by(feature, type, scenario, year, plot_order) %>% 
    summarise(mean = median(mean), .groups = "drop") %>% 
    left_join(mw_test, by = join_by(feature, type))

cols <- data.frame(
    labels = c("p < 0.05", " 0.05 \U2264 P < 0.1", "P \u2265 0.1"), 
    colors = c("#313695", "#8e0152", "#F99379"))
cols <- cols %>% filter(labels %in% unique(en_fig$significance)) %>% 
    pull(colors)

p1 <- ggplot(data = sp_analysis_fig) +
    geom_boxplot(aes(x = plot_order, 
                     y = mean, fill = type), outliers = FALSE, 
                 fatten = 1, show.legend = FALSE) +
    scale_fill_manual(values = c("#018571", '#a6611a')) +
    new_scale_fill() +
    geom_point(data = en_fig, 
               aes(x = plot_order, y = mean, fill = significance), 
               size = 1.6, shape = 21, color = "white", stroke = 0.2) +
    scale_fill_manual("Median of endangered species", values = cols) +
    coord_flip() +
    scale_x_discrete(labels = function(x){gsub("_P|_N", "", x)}) +
    theme_pubclean(base_family = "Merriweather", base_size = 11) + 
    xlab("") + ylab("SHAP value change") + 
    facet_wrap(
        ~factor(type, levels = c("P", "N"),
                labels = c("Baseline-favorable area", 
                           "Baseline-unfavorable area")),
        scales = "free") +
    ggtitle("(c)") +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black", size = 10),
          panel.grid.major.y = element_line(color = "white"),
          panel.grid.major.x = element_line(
              linetype = "dotted", color = "lightgrey"),
          strip.background = element_blank(),
          strip.text.x = element_text(hjust = 0.3),
          strip.text = element_text(face = "bold", size = 11),
          plot.margin = unit(c(-0.5, 0, 0, 0), "cm"),
          legend.direction = "horizontal",
          legend.position = "top",
          legend.title = element_text(size = 8), 
          legend.box.spacing = unit(rep(0, 4), "cm"),
          plot.title = element_text(
              family = "Merriweather", size = 10,
              face = "bold", vjust = -9)) +
    guides(fill = guide_legend(override.aes = list(size = 3)))

p1
ggsave('/Users/leisong/downloads/test.png', 
       width = 6.5, height = 3.8, dpi = 500, bg = "white")

# panel b
pts <- magnitude_periods %>% 
    filter(scenario == toupper(ssp)) %>% 
    mutate(driver = case_when(
        driver == "forest" ~ "FOR",
        driver == "human_impact" ~ "HLU",
        driver == "grassland" ~ "GRA",
        str_detect(driver, "bio") ~ toupper(driver)))

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

img <- image_read(file.path(fig_dir, "fig1_flow.png"))
g <- image_ggplot(img, interpolate = TRUE)

ggarrange(g, p, nrow = 2, heights = c(1, 3))

ggsave(file.path(fig_dir, "Figure4_mag_species.png"), 
       width = 6.5, height = 4.5, dpi = 500, bg = "white")
ggsave(fname, width = 6.5, height = 5, dpi = 500)
