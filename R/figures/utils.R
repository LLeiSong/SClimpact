# Make some pseudo figures for the workflow

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(flextable)
library(ggpubr)
library(showtext)
font_add_google('Merriweather')
showtext_auto()
showtext_opts(dpi = 500)

# Use one example
sp_dir <- "results/sdm"
sp <- "Myrmecophaga_tridactyla"

shap_base <- rast(
    file.path(sp_dir, sp, sprintf("shap_base_%s_1000.tif", sp)))

shap_fut <- rast(
    file.path(sp_dir, sp, sprintf("shap_ssp370_2041-2070_%s_1000.tif", sp)))

shap_base <- aggregate(shap_base, fact = 29)
shap_fut <- aggregate(shap_fut, fact = 29)

shap_change <- shap_fut - shap_base
# 1 2
# 3 4
turnovers <- ((shap_base >= 0) + 1) * ((shap_fut >= 0) + 3)
turnovers[turnovers %in% c(3, 8)] <- 0
types <- shap_base >= 0

shap_base <- as.data.frame(shap_base, xy = TRUE)
shap_fut <- as.data.frame(shap_fut, xy = TRUE)
shap_change <- as.data.frame(shap_change, xy = TRUE)
turnovers <- as.data.frame(turnovers, xy = TRUE)
types <- as.data.frame(types, xy = TRUE)

for (var in names(shap_base)[3:6]){
    dat <- shap_base %>% select(x, y, all_of(var))
    names(dat)[3] <- "value"
    ggplot() +
        geom_tile(
            data = dat, aes(x = x, y = y, fill = value),
            color = "white", show.legend = FALSE, linewidth = 0.4) +
        colorspace::scale_fill_continuous_diverging(
            palette = "Blue-Red 2") +
        coord_equal() +
        theme_void()
    ggsave(sprintf("docs/workflow/baseline_%s.svg", var), height = 4, width = 4)
    
    dat <- shap_fut %>% select(x, y, all_of(var))
    names(dat)[3] <- "value"
    ggplot() +
        geom_tile(
            data = dat, aes(x = x, y = y, fill = value),
            color = "white", show.legend = TRUE, linewidth = 0.4) +
        colorspace::scale_fill_continuous_diverging(
            palette = "Blue-Red 2") +
        coord_equal() +
        theme_void()
    ggsave(sprintf("docs/workflow/future_%s.svg", var), height = 4, width = 4)
    
    dat <- shap_change %>% select(x, y, all_of(var))
    names(dat)[3] <- "value"
    ggplot() +
        geom_tile(
            data = dat, aes(x = x, y = y, fill = value),
            color = "white", show.legend = FALSE, linewidth = 0.4) +
        colorspace::scale_fill_continuous_divergingx(
            palette = "BrBG", rev = TRUE) +
        coord_equal() +
        theme_void()
    ggsave(sprintf("docs/workflow/shift_%s.svg", var), height = 4, width = 4)
    
    dat <- turnovers %>% select(x, y, all_of(var))
    names(dat)[3] <- "value"
    dat <- dat %>% 
        mutate(value = factor(value, levels = c(0, 4, 6), labels = c(0, 4, 6)))
    ggplot() +
        geom_tile(
            data = dat, aes(x = x, y = y, fill = value),
            color = "white", show.legend = FALSE, linewidth = 0.4) +
        scale_fill_manual(values = c("#c9c9c9", "#a6611a", "#018571")) +
        coord_equal() +
        theme_void()
    ggsave(sprintf("docs/workflow/turnover_%s.svg", var), height = 4, width = 4)
    
    dat <- types %>% select(x, y, all_of(var))
    names(dat)[3] <- "value"
    dat <- dat %>% 
        mutate(value = factor(value, levels = c(0, 1), labels = c(0, 1)))
    ggplot() +
        geom_tile(
            data = dat, aes(x = x, y = y, fill = value),
            color = "white", show.legend = FALSE, linewidth = 0.4) +
        scale_fill_manual(values = c("#99CEC6", "#dadbe2")) +
        coord_equal() +
        theme_void()
    ggsave(sprintf("docs/workflow/unfavorable_%s.svg", var), height = 4, width = 4)
    
    ggplot() +
        geom_tile(
            data = dat, aes(x = x, y = y, fill = value),
            color = "white", show.legend = FALSE, linewidth = 0.4) +
        scale_fill_manual(values = c("#dadbe2", "#DAC0A2")) +
        coord_equal() +
        theme_void()
    ggsave(sprintf("docs/workflow/favorable_%s.svg", var), height = 4, width = 4)
    
    dat <- shap_fut %>% select(x, y, all_of(var))
    names(dat)[3] <- "value"
    dat <- dat %>% mutate(value = ifelse(value >= 0, value, 0))
    ggplot() +
        geom_tile(
            data = dat, aes(x = x, y = y, fill = value),
            color = "white", show.legend = FALSE, linewidth = 0.4) +
        colorspace::scale_fill_continuous_divergingx(
            palette = "BrBG") +
        coord_equal() +
        theme_void()
    ggsave(sprintf("docs/workflow/favoring_sp_%s.svg", var), height = 4, width = 4)
    
    dat <- shap_fut %>% select(x, y, all_of(var))
    names(dat)[3] <- "value"
    dat <- dat %>% mutate(value = ifelse(value <= 0, value, 0))
    ggplot() +
        geom_tile(
            data = dat, aes(x = x, y = y, fill = value),
            color = "white", show.legend = FALSE, linewidth = 0.4) +
        colorspace::scale_fill_continuous_divergingx(
            palette = "BrBG") +
        coord_equal() +
        theme_void()
    ggsave(sprintf("docs/workflow/unfavoring_sp_%s.svg", var), height = 4, width = 4)
}

definitions <- data.frame(
    Term = c("Favoring turnover", "Unfavoring turnover", "Magnitude shifts"),
    Definition = c(
        "For a species at a location, the SHAP value of a given variable changes from ≤ 0 at baseline to > 0 in the future.",
        "For a species at a location, the SHAP value of a given variable changes from > 0 at baseline to ≤ 0 in the future.",
        "For a species at a location, the change in the SHAP value of a given variable from baseline to future."),
    Interpretation = c(
        "A location where an environmental variable changes from a negative (or neutral) to a positive contribution to species suitability, indicating emerging local benefits from that variable.",
        "A location where an environmental variable changes from a positive to a negative (or neutral) contribution to species suitability, indicating a local loss of suitability associated with that variable.",
        "A change in the magnitude of a variable’s contribution to species suitability, indicating increasing or decreasing importance of that variable in shaping local suitability.")
)

ggtexttable(
    definitions, rows = NULL, 
    theme = ttheme(
        padding = unit(c(1, 2), "mm"),
        colnames.style = colnames_style(
            size = 8, fill = "transparent", parse = TRUE, 
            face = "plain", font = "Merriweather"),
        tbody.style = tbody_style(
            size = 8, color = "black", font = "Merriweather")))
ggsave('/Users/leisong/downloads/tbl.svg', height = 4, width = 8)
