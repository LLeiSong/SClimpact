# Load settings
source("R/figures/setting.R")

#### Analyze the areas of turnovers for each species ####
sp_turnovers <- lapply(species_list, function(sp){
    fname <- file.path(sp_analysis_dir, sprintf("changes_%s.csv", sp))
    numbers <- read.csv(fname) %>% 
        separate(scenario, c('scenario', 'year'), sep = "_") %>% 
        filter(class == "direction change") %>%
        select(-c(class, metrics))
    sp_sum <- numbers %>% 
        group_by(scenario, year, sp, feature) %>% 
        summarise(all = sum(value), .groups = "drop")
    left_join(numbers, sp_sum, by = join_by(scenario, year, sp, feature)) %>% 
        mutate(perct = value / all * 100) %>% select(-c(value, all))
}) %>% bind_rows() %>% 
    mutate(sp = gsub("_", " ", sp))

# Attach IUCN status to the species analysis results
sp_turnovers <- left_join(sp_turnovers, sp_names, by = c("sp" = "species"))

# Standardize the names 
sp_turnovers <- sp_turnovers %>% filter(feature %in% features) %>%
    mutate(feature = case_when(
        feature == "forest" ~ "FOR",
        feature == "human_impact" ~ "HLU",
        feature == "grassland" ~ "GRA",
        str_detect(feature, "bio") ~ toupper(feature)))

#### Analyze the magnitude shifts for each species ####
sp_shifts <- lapply(species_list, function(sp){
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
sp_shifts <- left_join(sp_shifts, sp_names, by = c("sp" = "species"))

# Standardize the names 
sp_shifts <- sp_shifts %>% filter(feature %in% features) %>%
    mutate(feature = case_when(
        feature == "forest" ~ "FOR",
        feature == "human_impact" ~ "HLU",
        feature == "grassland" ~ "GRA",
        str_detect(feature, "bio") ~ toupper(feature)))

#### Make figures and tables ####
##### Supplementary Table 1 ####

sp_turnovers_tosave <- sp_turnovers %>% 
    mutate(type = gsub(" to ", "2", type)) %>% 
    rename("Species" = sp, "Turnover" = type, 
           "Area (% of habitat)" = perct, "IUCN category" = category,
           "Variable" = feature, "Year" = year, "Scenario" = scenario,
           "Status" = status)
write.csv(
    sp_turnovers_tosave, file.path(tbl_dir, "supplementary_table1.csv"),
    row.names = FALSE)
rm(sp_turnovers_tosave)

##### Supplementary Table 2 ####
sp_shifts_tosave <- sp_shifts %>%
    mutate(type = case_when(
        type == "P" ~ "Favorable",
        type == "N" ~ "Unfavorable")) %>% 
    rename(
        "Species" = sp, "IUCN category" = category, "Variable" = feature, 
        "Baseline" = type, "Year" = year, "Scenario" = scenario, 
        "Status" = status, "SHAP value change (median)" = median,
        "SHAP value change (Q1)" = `1q`, "SHAP value change (Q3)" = `3q`)
write.csv(
    sp_shifts_tosave, file.path(tbl_dir, "supplementary_table2.csv"),
    row.names = FALSE, na = "")
rm(sp_shifts_tosave)

##### Figure 2 ####
sp_turnovers <- sp_turnovers %>% 
    # Only focus on these two types for figures
    filter(type %in% c("P to N", "N to P"))

# Parameter setting
ssp <- "ssp370"
time_period <- "2041-2070"

# Turnover
sp_turnovers_fig <- sp_turnovers %>% 
    filter(scenario == ssp & year == time_period)

# Add a common order based on the order of ALL status
feature_orders <- sp_turnovers_fig %>% 
    group_by(feature, type, scenario, year) %>%
    summarise(perct = median(perct), .groups = "drop") %>%
    mutate(plot_order = paste(feature, type, sep = "_")) %>% 
    arrange(type, perct) %>% pull(plot_order)

sp_turnovers_fig <- sp_turnovers_fig %>% 
    mutate(plot_order = paste(feature, type, sep = "_")) %>% 
    mutate(plot_order = factor(plot_order, levels = feature_orders,
                               labels = feature_orders))

# Check if the distribution of EN group is different from all
mw_test <- lapply(unique(sp_turnovers_fig$feature), function(driver){
    lapply(unique(sp_turnovers_fig$type), function(tp){
        dat <- sp_turnovers_fig %>% filter(feature == driver & type == tp)
        dat_en <- dat %>% filter(status == "EN")
        wt <- wilcox.test(dat$perct, dat_en$perct)
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

en_fig <- sp_turnovers_fig %>% 
    filter(status == "EN") %>% 
    group_by(feature, type, scenario, year, plot_order) %>% 
    summarise(perct = median(perct), .groups = "drop") %>% 
    left_join(mw_test, by = join_by(feature, type))

# Label some special values
maxs <- sp_turnovers_fig %>% 
    group_by(feature, type, plot_order) %>%
    summarise(max = boxplot(perct)$stats[5],
              .groups = "drop") %>% ungroup()

xlims <- data.frame(
    type = c("P to N", "N to P"),
    ymax = c(45, 25))

maxs <- left_join(maxs, xlims, by = "type") %>% 
    filter(max > ymax) %>% 
    mutate(label = sprintf("\U2192%.0f%s", max, "%"),
           plot_y = ymax)

# Define color table
cols <- data.frame(
    labels = c("p < 0.05", " 0.05 \U2264 P < 0.1", "P \u2265 0.1"), 
    colors = c("#313695", "#8e0152", "#F99379"))

# Label some special values
setting <- theme_pubclean(base_size = 12) + 
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black", size = 11),
          panel.grid.major.y = element_line(color = "white"),
          panel.grid.major.x = element_line(
              linetype = "dotted", color = "lightgrey"),
          strip.background = element_blank(),
          strip.text.x = element_text(hjust = 0.3),
          strip.text = element_text(face = "bold", size = 12),
          plot.margin = unit(rep(0.2, 4), "cm"),
          legend.direction = "horizontal",
          legend.position = "top",
          legend.title = element_text(size = 8), 
          legend.box.spacing = unit(rep(0, 4), "cm"),
          plot.title = element_text(
              size = 12,
              face = "bold", hjust = 0.5))

cols_p <- cols %>% 
    filter(
        labels %in% 
            unique(en_fig[en_fig$type == "P to N", ]$significance)) %>% 
    pull(colors)

drivers <- sp_turnovers_fig %>% filter(type == "P to N") %>% 
    arrange(plot_order) %>% pull(feature) %>% unique()
a_cols <- grps %>% arrange(match(var, drivers)) %>% pull(color)

p1 <- ggplot(data = sp_turnovers_fig %>% filter(type == "P to N")) +
    geom_boxplot(aes(x = plot_order, y = perct), 
                 fill = "#a6611a", outliers = FALSE, 
                 fatten = 1, show.legend = FALSE) +
    new_scale_fill() + ggtitle("Disfavoring turnover") +
    geom_point(data = en_fig %>% filter(type == "P to N"), 
               aes(x = plot_order, y = perct, fill = significance), 
               size = 1.6, shape = 21, color = "white", stroke = 0.2) +
    scale_fill_manual("Median of endangered species  ", values = cols_p) +
    scale_x_discrete(labels = function(x){gsub("_P to N|_N to P", "", x)}) +
    scale_y_continuous(labels = function(x){paste0(x, "%")}) +
    coord_flip(ylim = c(0, unlist(xlims[xlims$type == "P to N", c("ymax")]))) +
    xlab("") + ylab("Area (% of dispersable habitat)") + 
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    annotate("text", 
             x = (maxs %>% filter(type == "P to N") %>% 
                      pull(plot_order) %>% as.numeric()) - 16 + 0.4, 
             y = maxs %>% filter(type == "P to N") %>% pull(plot_y), 
             label = maxs %>% filter(type == "P to N") %>% pull(label),
             color = "black", size = 2) + setting +
    theme(axis.text.y = element_text(color = a_cols))

cols_n <- cols %>% 
    filter(labels %in% unique(en_fig[en_fig$type == "N to P", ]$significance)) %>% 
    pull(colors)

drivers <- sp_turnovers_fig %>% filter(type == "N to P") %>% 
    arrange(plot_order) %>% pull(feature) %>% unique()
a_cols <- grps %>% arrange(match(var, drivers)) %>% pull(color)

p2 <- ggplot(data = sp_turnovers_fig %>% filter(type == "N to P")) +
    geom_boxplot(aes(x = plot_order, y = perct), 
                 fill = "#018571", outliers = FALSE, 
                 fatten = 1, show.legend = FALSE) +
    new_scale_fill() + ggtitle("Favoring turnover") +
    geom_point(data = en_fig %>% filter(type == "N to P"), 
               aes(x = plot_order, y = perct, fill = significance), 
               size = 1.6, shape = 21, color = "white", stroke = 0.2) +
    scale_fill_manual("Median of endangered species  ", values = cols_p) +
    scale_x_discrete(labels = function(x){gsub("_P to N|_N to P", "", x)}) +
    scale_y_continuous(labels = function(x){paste0(x, "%")}) +
    coord_flip(ylim = c(0, unlist(xlims[xlims$type == "N to P", c("ymax")]))) +
    xlab("") + ylab("Area (% of dispersable habitat)") + 
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    annotate("text", 
             x = (maxs %>% filter(type == "N to P") %>% 
                      pull(plot_order) %>% as.numeric()) + 0.4, 
             y = maxs %>% filter(type == "N to P") %>% pull(plot_y), 
             label = maxs %>% filter(type == "N to P") %>% pull(label),
             color = "black", size = 2) + setting +
    theme(axis.text.y = element_text(color = a_cols))

p1 <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "none")

# shifts
sp_shifts_fig <- sp_shifts %>% 
    filter(scenario == ssp & year == time_period)

# Add a common order based on the order of ALL status
feature_orders <- sp_shifts_fig %>% 
    group_by(feature, type, scenario, year) %>%
    summarise(median = median(median), .groups = "drop") %>%
    mutate(plot_order = paste(feature, type, sep = "_")) %>% 
    mutate(median = ifelse(type == "P", -median, median)) %>% 
    arrange(type, median) %>% pull(plot_order)

sp_shifts_fig <- sp_shifts_fig %>% 
    mutate(plot_order = paste(feature, type, sep = "_")) %>% 
    mutate(plot_order = factor(plot_order, levels = feature_orders,
                               labels = feature_orders))

# Check if the distribution of EN group is different from all
mw_test <- lapply(unique(sp_shifts_fig$feature), function(driver){
    lapply(unique(sp_shifts_fig$type), function(tp){
        dat <- sp_shifts_fig %>% filter(feature == driver & type == tp)
        dat_en <- dat %>% filter(status == "EN")
        wt <- wilcox.test(dat$median, dat_en$median)
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

en_fig <- sp_shifts_fig %>% 
    filter(status == "EN") %>% 
    group_by(feature, type, scenario, year, plot_order) %>% 
    summarise(median = median(median), .groups = "drop") %>% 
    left_join(mw_test, by = join_by(feature, type))

# Label some special values
minmaxs <- sp_shifts_fig %>% 
    group_by(feature, type, scenario, year, plot_order) %>%
    summarise(max = boxplot(median)$stats[5], 
              min = boxplot(median)$stats[1], 
              .groups = "drop") %>% ungroup()

xlims <- data.frame(
    type = c("P", "N"),
    ymin = c(-0.1, -0.025),
    ymax = c(0.05, 0.05))

minmaxs <- left_join(minmaxs, xlims, by = "type") %>% 
    filter(max > ymax | min < ymin)

minmaxs <- minmaxs %>% filter(max > ymax) %>% 
    mutate(label = sprintf("\U2192%.2f", max),
           plot_y = ymax) %>% 
    rbind(minmaxs %>% filter(min < ymin) %>% 
              mutate(label = sprintf("%.2f\U2190", min),
                     plot_y = ymin))

cols_p <- cols %>% 
    filter(labels %in% unique(en_fig[en_fig$type == "P", ]$significance)) %>% 
    pull(colors)

drivers <- sp_shifts_fig %>% filter(type == "P") %>% 
    arrange(plot_order) %>% pull(feature) %>% unique()
a_cols <- grps %>% arrange(match(var, drivers)) %>% pull(color)

p2 <- ggplot(data = sp_shifts_fig %>% filter(type == "P")) +
    geom_hline(yintercept = 0, color = 'white') +
    geom_hline(yintercept = 0, color = '#e74c3c', linetype = "dotted") +
    geom_boxplot(aes(x = plot_order, y = median), 
                 fill = "#DAC0A2", outliers = FALSE, 
                 fatten = 1, show.legend = FALSE) +
    new_scale_fill() + ggtitle("Baseline-favorable area") +
    geom_point(data = en_fig %>% filter(type == "P"), 
               aes(x = plot_order, y = median, fill = significance), 
               size = 1.6, shape = 21, color = "white", stroke = 0.2) +
    scale_fill_manual("Median of endangered species  ", values = cols_p) +
    scale_x_discrete(labels = function(x){gsub("_P|_N", "", x)}) +
    coord_flip(ylim = unlist(xlims[xlims$type == "P", c('ymin', "ymax")])) +
    xlab("") + ylab("\u2206 SHAP") + 
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    annotate("text", 
             x = (minmaxs %>% filter(type == "P") %>% 
                      pull(plot_order) %>% as.numeric()) - 16 + 0.4, 
             y = minmaxs %>% filter(type == "P") %>% pull(plot_y), 
             label = minmaxs %>% filter(type == "P") %>% pull(label),
             color = "black", size = 2) + setting +
    theme(axis.text.y = element_text(color = a_cols))

cols_n <- cols %>% 
    filter(labels %in% unique(en_fig[en_fig$type == "N", ]$significance)) %>% 
    pull(colors)

drivers <- sp_shifts_fig %>% filter(type == "N") %>% 
    arrange(plot_order) %>% pull(feature) %>% unique()
a_cols <- grps %>% arrange(match(var, drivers)) %>% pull(color)

p3 <- ggplot(data = sp_shifts_fig %>% filter(type == "N")) +
    geom_hline(yintercept = 0, color = 'white') +
    geom_hline(yintercept = 0, color = '#e74c3c', linetype = "dotted") +
    geom_boxplot(aes(x = plot_order, y = median), 
                 fill = '#99CEC6', outliers = FALSE, 
                 fatten = 1, show.legend = FALSE) +
    new_scale_fill() + ggtitle("Baseline-unfavorable area") +
    geom_point(data = en_fig %>% filter(type == "N"), 
               aes(x = plot_order, y = median, fill = significance), 
               size = 1.6, shape = 21, color = "white", stroke = 0.2) +
    scale_fill_manual("Median of endangered species  ", values = cols_n) +
    scale_x_discrete(labels = function(x){gsub("_P|_N", "", x)}) +
    coord_flip(ylim = unlist(xlims[xlims$type == "N", c('ymin', "ymax")])) +
    xlab("") + ylab("\u2206 SHAP") + 
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    annotate("text", 
             x = (minmaxs %>% filter(type == "N") %>% 
                      pull(plot_order) %>% as.numeric()) + 0.4, 
             y = minmaxs %>% filter(type == "N") %>% pull(plot_y), 
             label = minmaxs %>% filter(type == "N") %>% pull(label),
             color = "black", size = 2) + setting +
    theme(axis.text.y = element_text(color = a_cols))

p2 <- ggarrange(p2, p3, ncol = 2, common.legend = TRUE, legend = "none")

# Make a label
lgd_orig <- ggplot() +
    geom_point(data = en_fig, 
               aes(x = plot_order, y = median, color = significance), 
               size = 1.6) +
    scale_color_manual(name = " Median of endangered species ", 
                       values = cols$colors) +
    theme_pubclean(base_size = 12) + 
    theme(legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10))

lgd <- ggplot() + 
    geom_point(aes(x = 0:2, y = c(0, 0, 1)), 
               color = "transparent") +
    annotation_custom(
        ggplotGrob(as_ggplot(get_legend(lgd_orig))), 
        xmin = 0, xmax = 2, ymin = 0, ymax = 0.5) +
    annotate(
        geom = "text", y = 0.8, 
        x = 0.5, label = "Temperature", color = "#e66101", 
        size = 3.5) +
    annotate(
        geom = "text", y = 0.8, x = 1, 
        label = "Precipitation", color = "#072ac8", 
        size = 3.5) +
    annotate(
        geom = "text", y = 0.8, x = 1.5, 
        label = "Landcover", color = "#1c2541", 
        size = 3.5) + theme_void()

img <- image_read(file.path(fig_dir, "fig2_flow.jpg"))
g <- image_ggplot(img, interpolate = TRUE)

ggarrange(g, p1, lgd, p2, nrow = 4, ncol = 1,
          labels = c("(a)", "(b)", "", "(c)"),
          font.label = list(
              size = 12, color = "black", face = "bold"),
          heights = c(2.8, 5, 0.8, 5))

ggsave(file.path(fig_dir, "Figure2_driver_species.png"), 
       width = 6.5, height = 7.5, dpi = 500, bg = "white")

##### Extended Data Table 1 ####
en_max_fig <- sp_turnovers_fig %>% 
    filter(status == "EN") %>% 
    mutate(type = case_when(
        type == "P to N" ~ "Disfavoring",
        type == "N to P" ~ "Favoring")) %>% 
    group_by(feature, type, scenario, year) %>% 
    arrange(-perct) %>% slice_head(n = 1) %>% 
    arrange(desc(type), desc(plot_order)) %>% 
    ungroup() %>% select(type, feature, sp, perct, category) %>% 
    mutate(perct = round(perct, 2)) %>% 
    rename("Variable" = feature, "Turnover" = type, "Species" = sp, 
           "Area (% of \ndispersable-habitat)" = perct, "IUCN category" = category)

en_max_fig %>% mutate(Species = sprintf(" %s ", Species)) %>% 
    flextable() %>% autofit() %>% 
    bg(bg = "white", part = "all") %>% 
    align(align = "left", part = "all") %>% 
    bold(part = "header") %>% 
    merge_at(i = 1:16, j = 1) %>% 
    merge_at(i = 17:32, j = 1) %>% 
    # Add some lines
    border(i = 16, j = 1:5, 
           border.bottom = fp_border(color = "gray")) %>% 
    save_as_image(file.path(fig_dir, "extended_data_table1.png"), res = 500)

##### Extended Data Table 2 ####
en_max_fig <- sp_shifts_fig %>% 
    mutate(median_order = ifelse(type == "P", -median, median)) %>% 
    filter(status == "EN") %>% 
    mutate(type = case_when(
        type == "P" ~ "Favorable area",
        type == "N" ~ "Unfavorable area")) %>% 
    group_by(feature, type, scenario, year) %>% 
    arrange(-median_order) %>% slice_head(n = 1) %>% 
    arrange(type, desc(type), desc(plot_order)) %>% 
    ungroup() %>% 
    mutate(label = sprintf(" %.2f (%.2f-%.2f) ", median, `1q`, `3q`)) %>% 
    select(type, feature, sp, label, category) %>% 
    rename("Variable" = feature, "Baseline" = type, "Species" = sp, 
           "SHAP value change\nMedian (Q1 - Q3)" = label, 
           "IUCN category" = category)

en_max_fig %>% mutate(Species = sprintf(" %s ", Species)) %>% 
    flextable() %>% autofit() %>% 
    bg(bg = "white", part = "all") %>% 
    align(align = "left", part = "all") %>% 
    align(i = 1, j = 4, align = "center", part = "header") %>% 
    bold(part = "header") %>% 
    merge_at(i = 1:16, j = 1) %>% 
    merge_at(i = 17:32, j = 1) %>% 
    # Add some lines
    border(i = 16, j = 1:5, 
           border.bottom = fp_border(color = "gray")) %>% 
    save_as_image(file.path(fig_dir, "extended_data_table2.png"), res = 500)

##### Extended Data Fig.1 and Fig.2 ####

# Fig. 1
inds <- c(letters[1:5], letters[5:8])
figs <- lapply(c("ssp126", "ssp370", "ssp585"), function(ssp){
    figs <- lapply(c("2011-2040", "2041-2070", "2071-2100"), 
                   function(time_period){
        sp_analysis_fig <- sp_turnovers %>% 
            filter(scenario == ssp & year == time_period)
        
        id1 <- which(c("2011-2040", "2041-2070", "2071-2100") == time_period)
        id2 <- which(c("ssp126", "ssp370", "ssp585") == ssp)
        ind <- inds[(id2 - 1) * 3 + id1]
        
        # Add a common order based on the order of ALL status
        feature_orders <- sp_analysis_fig %>% 
            group_by(feature, type, scenario, year) %>%
            summarise(perct = median(perct), .groups = "drop") %>%
            mutate(plot_order = paste(feature, type, sep = "_")) %>% 
            arrange(perct) %>% pull(plot_order)
        
        sp_analysis_fig <- sp_analysis_fig %>% 
            mutate(plot_order = paste(feature, type, sep = "_")) %>% 
            mutate(plot_order = factor(plot_order, levels = feature_orders,
                                       labels = feature_orders))
        
        # Check if the distribution of EN group is different from all
        mw_test <- lapply(unique(sp_analysis_fig$feature), function(driver){
            lapply(unique(sp_analysis_fig$type), function(tp){
                dat <- sp_analysis_fig %>% 
                    filter(feature == driver & type == tp)
                dat_en <- dat %>% filter(status == "EN")
                wt <- wilcox.test(dat$perct, dat_en$perct)
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
            summarise(perct = median(perct), .groups = "drop") %>% 
            left_join(mw_test, by = join_by(feature, type))
        
        cols <- data.frame(
            labels = c("p < 0.05", " 0.05 \U2264 P < 0.1", "P \u2265 0.1"), 
            colors = c("#313695", "#8e0152", "#F99379"))
        cols <- cols %>% filter(labels %in% unique(en_fig$significance)) %>% 
            pull(colors)
        
        drivers <- sp_shifts_fig %>% filter(type == "P") %>% 
            arrange(plot_order) %>% pull(feature) %>% unique()
        a_cols <- grps %>% arrange(match(var, drivers)) %>% pull(color)
        
        if (ssp == "ssp370" & time_period == "2041-2070"){
            NULL
        } else {
            ggplot(data = sp_analysis_fig) +
                geom_boxplot(aes(x = plot_order, y = perct, fill = type), 
                             outliers = FALSE, fatten = 1, coef = 0,
                             show.legend = FALSE) +
                scale_fill_manual(values = c("#018571", '#a6611a')) +
                new_scale_fill() +
                geom_point(data = en_fig, 
                           aes(x = plot_order, y = perct, fill = significance), 
                           size = 2, shape = 21, color = "white", stroke = 0.2,
                           show.legend = FALSE) +
                scale_fill_manual("", values = cols) +
                coord_flip() +
                scale_x_discrete(labels = function(x){
                    lb <- gsub("_N to P|_P to N", "", x)
                    cols <- grps[match(lb, grps$var), "color"]
                    sapply(1:length(lb), function(i) {
                        paste0("<span style = 'color: ", 
                               cols[i], ";'>", lb[i], "</span>")
                    })}) +
                scale_y_continuous(labels = function(x){paste0(x, "%")}) +
                theme_pubclean(base_size = 12) + 
                xlab("") + ylab("Area (% of dispersable habitat)") + 
                ggtitle(sprintf("%s. %s, %s", ind, ssp, time_period)) +
                facet_wrap(
                    ~factor(type, levels = c("P to N", "N to P"),
                            labels = c("Disfavoring turnover", 
                                       "Favoring turnover")),
                    scales = "free") +
                theme(axis.text = element_text(color = "black"),
                      axis.text.y = element_markdown(),
                      axis.title = element_text(size = 12),
                      panel.grid.major.y = element_line(color = "white"),
                      panel.grid.major.x = element_line(
                          linetype = "dotted", color = "lightgrey"),
                      strip.background = element_blank(),
                      strip.text.x = element_text(hjust = 0.5),
                      strip.text = element_text(face = "italic", size = 12),
                      plot.title = element_text(
                          face = "bold", size = 12, hjust = 0.5, 
                          margin = margin(b = -5, t = 10)),
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
        }
    })
})

figs <- do.call(c, figs)
figs <- figs[!sapply(figs, is.null)]

lgd <- ggplot() + 
    geom_point(aes(x = 0:2, y = c(0, 0, 1)), 
               color = "transparent") +
    annotation_custom(
        ggplotGrob(as_ggplot(get_legend(lgd_orig))), 
        xmin = 0, xmax = 2, ymin = 0, ymax = 0.5) +
    annotate(
        geom = "text", y = 0.7, 
        x = 0.5, label = "Temperature", color = "#e66101", 
        size = 3.5) +
    annotate(
        geom = "text", y = 0.7, x = 1, 
        label = "Precipitation", color = "#072ac8", 
        size = 3.5) +
    annotate(
        geom = "text", y = 0.7, x = 1.5, 
        label = "Landcover", color = "#1c2541", 
        size = 3.5) + theme_void()

ggarrange(ggarrange(plotlist = figs, nrow = 4, ncol = 2), lgd, nrow = 2,
          heights = c(12, 0.4))

ggsave(file.path(fig_dir, "extended_data_fig1.png"), 
       width = 9, height = 12, dpi = 500, bg = "white")

# Fig.2
inds <- c(letters[1:5], letters[5:8])
figs <- lapply(c("ssp126", "ssp370", "ssp585"), function(ssp){
    figs <- lapply(c("2011-2040", "2041-2070", "2071-2100"), 
                   function(time_period){
        sp_analysis_fig <- sp_shifts %>% 
            filter(scenario == ssp & year == time_period)
        
        id1 <- which(c("2011-2040", "2041-2070", "2071-2100") == time_period)
        id2 <- which(c("ssp126", "ssp370", "ssp585") == ssp)
        ind <- inds[(id2 - 1) * 3 + id1]
        
        # Add a common order based on the order of ALL status
        feature_orders <- sp_analysis_fig %>% 
            group_by(feature, type, scenario, year) %>%
            summarise(median = median(median), .groups = "drop") %>%
            mutate(plot_order = paste(feature, type, sep = "_")) %>% 
            mutate(median = ifelse(type == "P", -median, median)) %>% 
            arrange(type, median) %>% pull(plot_order)
        
        sp_analysis_fig <- sp_analysis_fig %>% 
            mutate(plot_order = paste(feature, type, sep = "_")) %>% 
            mutate(plot_order = factor(plot_order, levels = feature_orders,
                                       labels = feature_orders))
        
        # Check if the distribution of EN group is different from all
        mw_test <- lapply(unique(sp_analysis_fig$feature), function(driver){
            lapply(unique(sp_analysis_fig$type), function(tp){
                dat <- sp_analysis_fig %>% 
                    filter(feature == driver & type == tp)
                dat_en <- dat %>% filter(status == "EN")
                wt <- wilcox.test(dat$median, dat_en$median)
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
            summarise(median = median(median), .groups = "drop") %>% 
            left_join(mw_test, by = join_by(feature, type))
        
        cols <- data.frame(
            labels = c("p < 0.05", " 0.05 \U2264 P < 0.1", "P \u2265 0.1"), 
            colors = c("#313695", "#8e0152", "#F99379"))
        cols <- cols %>% filter(labels %in% unique(en_fig$significance)) %>% 
            pull(colors)
        
        if (ssp == "ssp370" & time_period == "2041-2070"){
            NULL
        } else {
            ggplot(data = sp_analysis_fig) +
                geom_hline(yintercept = 0, color = 'white') +
                geom_hline(yintercept = 0, color = '#e74c3c', linetype = "dotted") +
                geom_boxplot(aes(x = plot_order, y = median, fill = type), 
                             outliers = FALSE, fatten = 1, coef = 0,
                             show.legend = FALSE) +
                scale_fill_manual(values = c("#99CEC6", '#DAC0A2')) +
                new_scale_fill() +
                geom_point(data = en_fig, 
                           aes(x = plot_order, y = median, fill = significance), 
                           size = 2, shape = 21, color = "white", stroke = 0.2,
                           show.legend = FALSE) +
                scale_fill_manual("", values = cols) +
                coord_flip() +
                scale_x_discrete(labels = function(x){
                    lb <- gsub("_N|_P", "", x)
                    cols <- grps[match(lb, grps$var), "color"]
                    sapply(1:length(lb), function(i) {
                        paste0("<span style = 'color: ", 
                               cols[i], ";'>", lb[i], "</span>")
                    })}) +
                theme_pubclean(base_size = 12) + 
                xlab("") + ylab("\u2206 SHAP") + 
                ggtitle(sprintf("%s. %s, %s", ind, ssp, time_period)) +
                facet_wrap(
                    ~factor(type, levels = c("P", "N"),
                            labels = c("Baseline-favorable area", 
                                       "Baseline-unfavorable area")),
                    scales = "free") +
                theme(axis.text = element_text(color = "black"),
                      axis.text.y = element_markdown(),
                      axis.title = element_text(size = 12),
                      panel.grid.major.y = element_line(color = "white"),
                      panel.grid.major.x = element_line(
                          linetype = "dotted", color = "lightgrey"),
                      strip.background = element_blank(),
                      strip.text.x = element_text(hjust = 0.5),
                      strip.text = element_text(face = "italic", size = 10),
                      plot.title = element_text(
                          face = "bold", size = 12, hjust = 0.5, 
                          margin = margin(b = -5, t = 10)),
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
        }
    })
})

figs <- do.call(c, figs)
figs <- figs[!sapply(figs, is.null)]

ggarrange(ggarrange(plotlist = figs, nrow = 4, ncol = 2), lgd, nrow = 2,
          heights = c(12, 0.4))

ggsave(file.path(fig_dir, "extended_data_fig2.png"), 
       width = 9, height = 12, dpi = 500, bg = "white")

# Clean up
rm(list = ls()); gc()
