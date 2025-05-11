# Load settings
source("R/figures/setting.R")

#### Analyze the areas of turnovers for each species ####
sp_analysis <- lapply(species_list, function(sp){
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
sp_analysis <- left_join(sp_analysis, sp_names, by = c("sp" = "species"))

# Standardize the names 
sp_analysis <- sp_analysis %>% filter(feature %in% features) %>%
    mutate(feature = case_when(
        feature == "forest" ~ "FOR",
        feature == "human_impact" ~ "HLU",
        feature == "grassland" ~ "GRA",
        str_detect(feature, "bio") ~ toupper(feature)))

#### Make figures and tables ####
##### Supplementary Table 1 ####

sp_analysis_tosave <- sp_analysis %>% 
    rename("Species" = sp, "Turnover" = type, 
           "Area (% of habitat)" = perct, "IUCN category" = category,
           "Variable" = feature, "Year" = year, "Scenario" = scenario,
           "Status" = status)
write.csv(
    sp_analysis_tosave, file.path(tbl_dir, "supplementary_table1.csv"),
    row.names = FALSE)
rm(sp_analysis_tosave)

##### Figure 1 ####
sp_analysis <- sp_analysis %>% 
    # Only focus on these two types for figures
    filter(type %in% c("P to N", "N to P"))

# Parameter setting
ssp <- "ssp370"
time_period <- "2041-2070"

sp_analysis_fig <- sp_analysis %>% 
    filter(scenario == ssp & year == time_period)

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
        dat <- sp_analysis_fig %>% filter(feature == driver & type == tp)
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

p <- ggplot(data = sp_analysis_fig) +
    geom_boxplot(aes(x = plot_order, 
                     y = perct, fill = type), outliers = FALSE, 
                 fatten = 1, show.legend = FALSE) +
    scale_fill_manual(values = c("#018571", '#a6611a')) +
    new_scale_fill() +
    geom_point(data = en_fig, 
               aes(x = plot_order, y = perct, fill = significance), 
               size = 1.6, shape = 21, color = "white", stroke = 0.2) +
    scale_fill_manual("Median of endangered species", 
                      values = c("#313695", "#8e0152", "#F99379")) +
    coord_flip() +
    scale_x_discrete(labels = function(x){gsub("_N to P|_P to N", "", x)}) +
    scale_y_continuous(labels = function(x){paste0(x, "%")}) +
    theme_pubclean(base_family = "Merriweather", base_size = 11) + 
    xlab("") + ylab("") + 
    facet_wrap(
        ~factor(type, levels = c("P to N", "N to P"),
                labels = c("(b) P2N turnover", "(c) N2P turnover")),
        scales = "free") +
    theme(axis.text = element_text(color = "black"),
          panel.grid.major.y = element_line(color = "white"),
          panel.grid.major.x = element_line(
              linetype = "dotted", color = "lightgrey"),
          strip.background = element_blank(),
          strip.text.x = element_text(hjust = 0.3),
          strip.text = element_text(face = "bold", size = 11),
          plot.margin = unit(c(0, 0, -0.3, 0), "cm"),
          legend.direction = "vertical",
          legend.position = "inside",
          legend.title = element_text(size = 8), 
          legend.position.inside = c(0.32, 0.3),
          legend.spacing = unit(rep(0, 4), "cm")) +
    guides(fill = guide_legend(override.aes = list(size = 3)))

img <- image_read(file.path(fig_dir, "fig1_flow.png"))
g <- image_ggplot(img, interpolate = TRUE)

ggarrange(g, NULL, p, nrow = 3, heights = c(1, 0.1, 3))

ggsave(file.path(fig_dir, "Figure1_driver_species.png"), 
       width = 6.5, height = 4.6, dpi = 500, bg = "white")

##### Extended Data Table 1 ####
en_max_fig <- sp_analysis_fig %>% 
    filter(status == "EN") %>% 
    mutate(type = case_when(
        type == "P to N" ~ "P2N",
        type == "N to P" ~ "N2P")) %>% 
    group_by(feature, type, scenario, year) %>% 
    arrange(-perct) %>% slice_head(n = 1) %>% 
    arrange(desc(type), desc(plot_order)) %>% 
    ungroup() %>% select(type, feature, sp, perct, category) %>% 
    mutate(perct = round(perct, 2)) %>% 
    rename("Variable" = feature, "Turnover" = type, "Species" = sp, 
           "Area (% of habitat)" = perct, "IUCN category" = category)

en_max_fig %>% mutate(Species = sprintf(" %s ", Species)) %>% 
    flextable() %>% autofit() %>% 
    align(align = "left", part = "all") %>% 
    bold(part = "header") %>% 
    font(fontname = "Merriweather", part = "all") %>% 
    merge_at(i = 1:16, j = 1) %>% 
    merge_at(i = 17:32, j = 1) %>% 
    # Add some lines
    border(i = 16, j = 1:5, 
           border.bottom = fp_border(color = "gray")) %>% 
    save_as_image(file.path(fig_dir, "extended_data_table1.png"), res = 500)

##### Extended Data Fig.1 ####

figs <- lapply(c("ssp126", "ssp370", "ssp585"), function(ssp){
    figs <- lapply(c("2011-2040", "2041-2070", "2071-2100"), 
                   function(time_period){
        sp_analysis_fig <- sp_analysis %>% 
            filter(scenario == ssp & year == time_period)
        
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
        
        if (ssp == "ssp370" & time_period == "2041-2070"){
            g <- ggplot(data = sp_analysis_fig) +
                geom_boxplot(aes(x = plot_order, y = perct, fill = type), 
                             outliers = FALSE, fatten = 1) +
                scale_fill_manual("Turnover", values = c("#018571", '#a6611a')) +
                new_scale_fill() +
                geom_point(data = en_fig, 
                           aes(x = plot_order, y = perct, fill = significance), 
                           size = 2, shape = 21, color = "white", stroke = 0.2) +
                scale_fill_manual("Median of endangered species", values = cols) +
                coord_flip() +
                scale_x_discrete(labels = function(x){
                    gsub("_N to P|_P to N", "", x)}) +
                scale_y_continuous(labels = function(x){paste0(x, "%")}) +
                theme_pubclean(base_family = "Merriweather", base_size = 11) + 
                xlab("") + ylab("") + 
                facet_wrap(
                    ~factor(type, levels = c("P to N", "N to P"),
                            labels = c(sprintf("(%s %s)\nP2N turnover", 
                                               ssp, time_period), 
                                       sprintf("(%s %s)\nN2P turnover", 
                                               ssp, time_period))),
                    scales = "free") +
                theme(axis.text = element_text(color = "black"),
                      panel.grid.major.y = element_line(color = "white"),
                      panel.grid.major.x = element_line(
                          linetype = "dotted", color = "lightgrey"),
                      strip.background = element_blank(),
                      strip.text.x = element_text(hjust = 0.3),
                      strip.text = element_text(face = "bold", size = 10),
                      legend.direction = "vertical",
                      legend.position = "inside",
                      legend.position.inside = c(0.5, 0.5),
                      legend.text = element_text(size = 13),
                      legend.title = element_text(size = 13)) +
                guides(fill = guide_legend(override.aes = list(size = 6)))
            
            as_ggplot(get_legend(g))
        } else {
            ggplot(data = sp_analysis_fig) +
                geom_boxplot(aes(x = plot_order, y = perct, fill = type), 
                             outliers = FALSE, fatten = 1, show.legend = FALSE) +
                scale_fill_manual(values = c("#018571", '#a6611a')) +
                new_scale_fill() +
                geom_point(data = en_fig, 
                           aes(x = plot_order, y = perct, fill = significance), 
                           size = 2, shape = 21, color = "white", stroke = 0.2,
                           show.legend = FALSE) +
                scale_fill_manual("", values = cols) +
                coord_flip() +
                scale_x_discrete(labels = function(x){
                    gsub("_N to P|_P to N", "", x)}) +
                scale_y_continuous(labels = function(x){paste0(x, "%")}) +
                theme_pubclean(base_family = "Merriweather", base_size = 11) + 
                xlab("") + ylab("") + 
                facet_wrap(
                    ~factor(type, levels = c("P to N", "N to P"),
                            labels = c(sprintf("(%s %s)\nP2N turnover", 
                                               ssp, time_period), 
                                       sprintf("(%s %s)\nN2P turnover", 
                                               ssp, time_period))),
                    scales = "free") +
                theme(axis.text = element_text(color = "black"),
                      panel.grid.major.y = element_line(color = "white"),
                      panel.grid.major.x = element_line(
                          linetype = "dotted", color = "lightgrey"),
                      strip.background = element_blank(),
                      strip.text.x = element_text(hjust = 0.3),
                      strip.text = element_text(face = "bold", size = 10),
                      plot.margin = unit(c(0, 0, -0.3, 0), "cm"))
        }
    })
})

figs <- do.call(c, figs)

ggarrange(plotlist = figs, nrow = 3, ncol = 3)

ggsave(file.path(fig_dir, "extended_data_fig1.png"), 
       width = 18, height = 12, dpi = 500, bg = "white")

# Clean up
rm(list = ls()); gc()
