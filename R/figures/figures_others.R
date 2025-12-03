# Load settings
source("R/figures/setting.R")

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
            size = 10, fill = "transparent", parse = TRUE, 
            face = "plain"),
        tbody.style = tbody_style(
            size = 10, color = "black"))))

g <- ggplot(data = evals %>% filter(type == "testing"), 
            aes(x = value, after_stat(density), fill = metric.eval)) +
    geom_density(alpha = 0.8) +
    xlab("Metric value(0 - 1)") +
    ylab("Density") +
    scale_fill_brewer(name = "Metric", palette = "Dark2") +
    theme_pubclean(base_size = 12) +
    theme(axis.text = element_text(color = "black"))

ggarrange(g, tbl, nrow = 1, widths = c(0.7, 0.3))

ggsave(file.path(fig_dir, "Figure_s_model_eval.png"),
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
    mutate(var = ifelse(var == "forest", "Forest coverage (FOR)",
                        ifelse(var == "human_impact", "Human land use (HLU)",
                               gsub("bio", "BIO", var)))) %>% 
    mutate(var = ifelse(var == "grassland", "Grassland (GRA)", var)) %>% 
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
    theme_pubclean(base_size = 12) +
    theme(axis.text.x = element_text(
        color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(color = "black", size = 10),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(vjust = -1),
        legend.position = "none")

ggsave(file.path(fig_dir, "Figure_s_vars_selected.png"),
       width = 4.5, height = 4, dpi = 500, bg = "white")

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
    theme_pubclean(base_size = 12) +
    theme(axis.text.x = element_text(
        color = "black"),
        axis.text.y = element_text(color = "black"))

ggsave(file.path(fig_dir, "Figure_s_nshap.png"),
       width = 3, height = 4, dpi = 500, bg = "white")

#### Justification for SHAP ####
# Load all statistics
species_list <- read.csv(
    file.path(root_dir, "results/species_reliable_final.csv"))

evals_omni <- read.csv(file.path(root_dir, "results/validate/evals_raw.csv")) %>% 
    left_join(., species_list %>% select(species, range), 
              by = c("sp" = "species"))

# Rename the item
evals_omni <- evals_omni %>% 
    mutate(item = factor(
        item, 
        levels = c("dat", "dat_less", "bg_bigger", "bg_smaller"),
        labels = c("Matched domain", "Matched domain with less samples",
                   "Bigger domain", "Smaller domain"))) %>% 
    rename(Background = item)

evals_omni %>% group_by(Background) %>% 
    summarise(r = median(r, na.rm = TRUE), iou = median(iou, na.rm = TRUE))

# Filter out outliers
evals_filtered <- lapply(unique(evals_omni$Background), function(x){
    evals_omni[evals_omni$Background == x, ][
        !is_outlier(evals_omni[evals_omni$Background == x, ]$r), ]
}) %>% bind_rows()

# Get the bad example, the name is pseudo without any meaning
evals_bad <- evals_omni %>% filter(sp == "Ovis_gmelini")

# Figure
ggplot(evals_omni, aes(x = r, y = iou, color = Background)) +
    geom_smooth(
        data = evals_filtered, aes(x = r, y = iou, color = Background),
        method = lm, se = FALSE, linewidth = 1) +
    geom_point(size = 1) + scale_color_npg() +
    geom_point(data = evals_bad, aes(x = r, y = iou, color = Background), 
               shape = 18, size = 4, show.legend = FALSE) +
    labs(x = "Correlation", y = "IoU") +
    theme_pubclean(base_size = 11) + 
    theme(axis.text = element_text(color = "black", size = 11),
          axis.title = element_text(color = "black", size = 11),
          legend.text = element_text(color = "black", size = 11),
          panel.grid.major.x = element_line(
              linetype = "dotted", color = "lightgrey")) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(fig_dir, "Figure_s_shap.png"), 
       width = 6.5, height = 5.5, dpi = 300, bg = "white")

# Load the species
load(file.path(
    root_dir, "results/validate/examples", 
    "Ovis_gmelini", "v_sp.rda"))

shaps_dat <- rast(file.path(
    root_dir, "results/validate/examples", "Ovis_gmelini", "shaps_dat.tif"))

suits_dat <- rast(file.path(
    root_dir, "results/validate/examples", "Ovis_gmelini", "suits_dat.tif"))

load("results/validate/examples/Ovis_gmelini/acg_bg_dat.rda")

map_theme <- theme_pubclean(base_size = 10) + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(color = "black", size = 10),
          legend.text = element_text(color = "black", size = 10),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.key.height = unit(0.3, "cm"))

g_shaps <- purrr::map(
    seq_len(nlyr(shaps_dat)),
    function(x) {
        ggplot() +
            geom_spatraster(data = shaps_dat[[x]]) +
            scale_fill_whitebox_c(
                name = "Shapley\nvalue", palette = "viridi", 
                na.value = "transparent",
                labels = scales::label_number(accuracy = 0.1)) +
            labs(title = toupper(names(shaps_dat[[x]]))) +
            map_theme
    }
)

g_suits <- purrr::map(
    seq_len(nlyr(suits_dat)),
    function(x) {
        ggplot() +
            geom_spatraster(data = suits_dat[[x]]) +
            scale_fill_whitebox_c(
                name = "Suitability\n(0-1)", palette = "viridi", 
                na.value = "transparent",
                labels = scales::label_number(accuracy = 0.1)) +
            labs(title = toupper(names(shaps_dat[[x]]))) +
            map_theme
    }
)

cowplot::plot_grid(plotlist = c(g_shaps, g_suits), nrow = 2)

ggsave(file.path(fig_dir, "Figure_s_exp.png"), 
       width = 6.5, height = 4.5, dpi = 300, bg = "white")

shaps_dat <- shaps_dat >= 0
g_shaps <- purrr::map(
    seq_len(nlyr(shaps_dat)),
    function(x) {
        ggplot() +
            geom_spatraster(data = shaps_dat[[x]]) +
            scale_fill_manual(
                name = "Binary (Shapley)", values = c("#bde0fe", "#003049"),
                labels = c(0, 1), na.translate = FALSE) +
            labs(title = toupper(names(shaps_dat[[x]]))) +
            map_theme
    }
)

suits_dat <- suits_dat >= avgs_bg
g_suits <- purrr::map(
    seq_len(nlyr(suits_dat)),
    function(x) {
        ggplot() +
            geom_spatraster(data = suits_dat[[x]]) +
            scale_fill_manual(
                name = "Binary (Suitability)", values = c("#bde0fe", "#003049"),
                labels = c(0, 1), na.translate = FALSE) +
            labs(title = toupper(names(shaps_dat[[x]]))) +
            map_theme
    }
)

cowplot::plot_grid(plotlist = c(g_shaps, g_suits), nrow = 2)

ggsave(file.path(fig_dir, "Figure_s_exp_binary.png"), 
       width = 6.5, height = 4.5, dpi = 300, bg = "white")

#### Justification for SHAP and RF ####
evals_rf <- read.csv(file.path(root_dir, "results/validate/evals_rf_raw.csv")) %>% 
    left_join(., species_list %>% select(species, range), 
              by = c("sp" = "species"))

evals_rf <- evals_rf %>% 
    mutate(item = factor(
        item, 
        levels = c("dat", "dat_less", "bg_bigger", "bg_smaller"),
        labels = c("Matched domain", "Matched domain with less samples",
                   "Bigger domain", "Smaller domain"))) %>% 
    rename(Background = item)

evals_rf %>% group_by(Background) %>% 
    summarise(r = median(r, na.rm = TRUE), iou = median(iou, na.rm = TRUE))

evals_bad <- evals_rf %>% filter(sp == "Platyrrhinus_dorsalis")

# Figure
ggplot(evals_rf, aes(x = r, y = iou, color = Background)) +
    geom_smooth(
        data = evals_rf, aes(x = r, y = iou, color = Background),
        method = lm, se = FALSE, linewidth = 1) +
    geom_point(size = 1) + scale_color_npg() +
    geom_point(data = evals_bad, aes(x = r, y = iou, color = Background),
               shape = 18, size = 4, show.legend = FALSE) +
    labs(x = "Correlation", y = "IoU") +
    theme_pubclean(base_size = 11) + 
    theme(axis.text = element_text(color = "black", size = 11),
          axis.title = element_text(color = "black", size = 11),
          legend.text = element_text(color = "black", size = 11),
          panel.grid.major.x = element_line(
              linetype = "dotted", color = "lightgrey")) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(fig_dir, "Figure_s_shap_rf.png"), 
       width = 6.5, height = 5.2, dpi = 300, bg = "white")

sp <- "Platyrrhinus_dorsalis"

load(file.path(
    root_dir, "results/validate/examples", sp, "v_sp.rda"))

shaps_dat <- rast(file.path(
    root_dir, "results/validate/examples", sp, "shaps_dat.tif"))

suits_dat <- rast(file.path(
    root_dir, "results/validate/examples", sp, "suits_dat.tif"))

load(file.path("results/validate/examples", sp, "acg_bg_dat.rda"))

map_theme <- theme_pubclean(base_size = 10) + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(color = "black", size = 10),
          legend.text = element_text(color = "black", size = 10),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.key.height = unit(0.3, "cm"))

g_shaps <- purrr::map(
    seq_len(nlyr(shaps_dat)),
    function(x) {
        ggplot() +
            geom_spatraster(data = shaps_dat[[x]]) +
            scale_fill_whitebox_c(
                name = "Shapley\nvalue", palette = "viridi", 
                na.value = "transparent",
                labels = scales::label_number(accuracy = 0.1)) +
            labs(title = toupper(names(shaps_dat[[x]]))) +
            map_theme
    }
)

g_suits <- purrr::map(
    seq_len(nlyr(suits_dat)),
    function(x) {
        ggplot() +
            geom_spatraster(data = suits_dat[[x]]) +
            scale_fill_whitebox_c(
                name = "Suitability\n(0-1)", palette = "viridi", 
                na.value = "transparent",
                labels = scales::label_number(accuracy = 0.1)) +
            labs(title = toupper(names(shaps_dat[[x]]))) +
            map_theme
    }
)

cowplot::plot_grid(plotlist = c(g_shaps, g_suits), nrow = 2)

ggsave(file.path(fig_dir, "Figure_s_exp_rf.png"), 
       width = 6.5, height = 4.5, dpi = 300, bg = "white")

shaps_dat <- shaps_dat >= 0
g_shaps <- purrr::map(
    seq_len(nlyr(shaps_dat)),
    function(x) {
        ggplot() +
            geom_spatraster(data = shaps_dat[[x]]) +
            scale_fill_manual(
                name = "Binary (Shapley)", values = c("#bde0fe", "#003049"),
                labels = c(0, 1), na.translate = FALSE) +
            labs(title = toupper(names(shaps_dat[[x]]))) +
            map_theme
    }
)

suits_dat <- suits_dat >= avgs_bg
g_suits <- purrr::map(
    seq_len(nlyr(suits_dat)),
    function(x) {
        ggplot() +
            geom_spatraster(data = suits_dat[[x]]) +
            scale_fill_manual(
                name = "Binary (Suitability)", values = c("#bde0fe", "#003049"),
                labels = c(0, 1), na.translate = FALSE) +
            labs(title = toupper(names(shaps_dat[[x]]))) +
            map_theme
    }
)

cowplot::plot_grid(plotlist = c(g_shaps, g_suits), nrow = 2)

ggsave(file.path(fig_dir, "Figure_s_exp_binary_rf.png"), 
       width = 6.5, height = 4.5, dpi = 300, bg = "white")
