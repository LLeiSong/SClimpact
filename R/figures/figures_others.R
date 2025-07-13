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
