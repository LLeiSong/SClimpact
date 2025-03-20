shaps <- lapply(species_list, function(sp){
    read.csv(
        file.path(root_dir, "results/sdm", sp, 
                  sprintf("shap_cor_%s.csv", sp))) %>% 
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

ggsave("docs/figures/nshap.png",
       width = 3, height = 4, dpi = 500)
