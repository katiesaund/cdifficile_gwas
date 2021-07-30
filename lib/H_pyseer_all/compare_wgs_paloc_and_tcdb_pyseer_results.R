library(tidyverse)
library(forcats)

tcdb_df <- read_tsv("../../data/15_pyseer/tcdb_model_results/best_model_results.tsv")
paloc_df <- read_tsv("../../data/15_pyseer/paloc_model_results/best_model_results.tsv")
wgs_df <- read_tsv("../../data/15_pyseer/wgs_model_results/best_model_results.tsv")

wgs_df <- wgs_df %>% 
  filter(phenotype == "log_toxin") %>% 
  mutate(genotype = "wgs")

df <- rbind(wgs_df, paloc_df, tcdb_df)

set.seed(1) # to get a reproducible jitter distribution where all observed values aren't covered by negative controls
pyseer_plot <- df %>% 
  mutate(plot_order = ifelse(genotype == "paloc", 2, 1)) %>% 
  mutate(plot_order = ifelse(genotype == "wgs", 3, plot_order)) %>% 
  ggplot(mapping = aes(y = Rsquared,
                       x = fct_reorder(genotype, plot_order), 
                       col = negative_control)) + 
  geom_jitter(width = 0.15, height = 0) +  
  scale_color_manual(values = c("grey", "black")) + 
  theme_bw() + 
  # ggtitle(unique(best_results_df$longhand)) + 
  xlab("") + 
  scale_x_discrete(labels = c(expression(italic("tcdB")), "PaLoc", "WG")) + 
  ylim(c(0, 1)) + 
  ylab(expression(italic("R"^2))) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.2, 0.9), 
        legend.background = element_rect(fill = "white", colour = "black")) 

ggsave(plot = pyseer_plot,
       filename = paste0("../../figures/pyseer_results/all_three_EN_models_jitter_best_alpha_log_toxin.png"), 
       units = "in",
       height = 4,
       width = 4)
save(pyseer_plot, 
     file = paste0("../../figures/pyseer_results/all_three_EN_models_jitter_best_alpha_log_toxin.RData"))
