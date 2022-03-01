# Make figure 2: 
# Phenotype phylogenetic signal and genomic model.
# (A) The phylogenetic signal of each phenotype (black) and its negative controls (grey). WN = white noise. BM = Brownian motion.
# (B) Elastic nets modeling each phenotype were generated with SNPs, indels, and accessory genes as input variables. 

# Plot for 2A was originally made in: ../summarize_all_pheno_data/5B_plot_phylogenetic_signal.R
# Plot for 2B war originally made in: ../pyseer_analysis/4_plot_and_compare_model_performance.R
# library(aplot)
library(tidyverse)
library(patchwork)

# Recall: 
# A continuous phenotype that is modelled well by Brownian motion has a λ near 1
#    while a white noise phenotype has a λ near 0 [19]
# a binary phenotype that is modelled well by Brownian motion has a d near 0 
#    while a white noise phenotype has a D near 1 [20]"
# Now with negative controls
phylo_signal_w_neg_ctrl_df <- read_tsv("../../data/4_phenotypes/phenotypic_analysis_results/phylo_signal_with_neg_ctrls.tsv", 
                                       col_names = TRUE)
phylo_signal_w_neg_ctrl_df <- pivot_longer(phylo_signal_w_neg_ctrl_df, 
                                           cols = c(phylo_signal, neg_ctrl_phylo_signal_1, 
                                                    neg_ctrl_phylo_signal_2, 
                                                    neg_ctrl_phylo_signal_3), 
                                           names_to = "phylo_signal_type")
phylo_signal_w_neg_ctrl_df <- phylo_signal_w_neg_ctrl_df %>% 
  mutate("negative_ctrl" = grepl("neg", phylo_signal_type)) %>% 
  mutate(pretty_name = gsub("ln[(]", "", longhand) %>% gsub("[)]", "", .)) %>% 
  mutate(pretty_name = gsub(" and ", "&", pretty_name)) %>% 
  mutate(plot_order = 1)

phylo_signal_w_neg_ctrl_df$negative_ctrl[phylo_signal_w_neg_ctrl_df$negative_ctrl == TRUE] <- "Negative control"
phylo_signal_w_neg_ctrl_df$negative_ctrl[phylo_signal_w_neg_ctrl_df$negative_ctrl == FALSE] <- "Observed value"

phylo_signal_w_neg_ctrl_df <- 
  phylo_signal_w_neg_ctrl_df %>% 
  mutate(plot_name = case_when(pretty_name == "germination in Tc" ~ "Germ. in TCA", 
                               pretty_name == "germination in Tc&Gly" ~ "Germ. in TCA&Gly", 
                               pretty_name == "fqR" ~ "Fluoroquinolone Resistance", 
                               pretty_name == "spore viability" ~ "Spore viability", 
                               pretty_name == "growth rate" ~ "Max growth rate", 
                               pretty_name == "# spores" ~ "Total spores",
                               pretty_name == "toxin activity" ~ "Cytotoxicity", 
                               pretty_name == "Severity" ~ "Severe infection"))

# Want to change plot order to: 
# 1. germ tc
# 2. germ TG
# 3. growth
# 4. # spores
# 5. spore viability
# 6. toxin 
phylo_signal_w_neg_ctrl_df$plot_order[phylo_signal_w_neg_ctrl_df$pretty_name == "germination in Tc"] <- 1
phylo_signal_w_neg_ctrl_df$plot_order[phylo_signal_w_neg_ctrl_df$pretty_name == "germination in Tc&Gly"] <- 2
phylo_signal_w_neg_ctrl_df$plot_order[phylo_signal_w_neg_ctrl_df$pretty_name == "growth rate"] <- 3
phylo_signal_w_neg_ctrl_df$plot_order[phylo_signal_w_neg_ctrl_df$pretty_name == "# spores"] <- 4
phylo_signal_w_neg_ctrl_df$plot_order[phylo_signal_w_neg_ctrl_df$pretty_name == "spore viability"] <- 5
phylo_signal_w_neg_ctrl_df$plot_order[phylo_signal_w_neg_ctrl_df$pretty_name == "toxin activity"] <- 6

set.seed(6) # to get a reproducible jitter distribution where all observed values aren't covered by negative controls
dot_plot_wo_leg <- phylo_signal_w_neg_ctrl_df %>% 
  filter(bin_cont == "continuous") %>% 
  ggplot(aes(x = fct_reorder(plot_name, plot_order), y = value, col = negative_ctrl)) + 
  geom_jitter(width = 0.25, height = 0) + 
  theme_bw() + 
  xlab("") + 
  ylab(expression(paste("Phylogenetic Signal (", lambda, ")"))) + 
  scale_color_manual(values = c("grey", "black"), guide = "none") + 
  theme(legend.title = element_blank()) + 
  ylim(c(0, 1.0)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(tag = "A")

# 2B
best_results_df <- read_tsv("../../data/15_pyseer/wgs_model_results/best_model_results.tsv", col_names = TRUE)
best_results_df <- best_results_df %>% 
  filter(bin_cont == "continuous") %>% 
  mutate(pretty_name = gsub("ln[(]", "", longhand) %>% gsub("[)]", "", .)) %>% 
  mutate(pretty_name = gsub(" and ", "&", pretty_name))  %>% 
  mutate(plot_order = 1)
best_results_df$plot_order[best_results_df$pretty_name == "germination in Tc"] <- 1
best_results_df$plot_order[best_results_df$pretty_name == "germination in Tc&Gly"] <- 2
best_results_df$plot_order[best_results_df$pretty_name == "growth rate"] <- 3
best_results_df$plot_order[best_results_df$pretty_name == "# spores"] <- 4
best_results_df$plot_order[best_results_df$pretty_name == "spore viability"] <- 5
best_results_df$plot_order[best_results_df$pretty_name == "toxin activity"] <- 6

best_results_df <- 
  best_results_df %>% 
  mutate(plot_name = case_when(pretty_name == "germination in Tc" ~ "Germ. in TCA", 
                               pretty_name == "germination in Tc&Gly" ~ "Germ. in TCA&Gly", 
                               pretty_name == "fqR" ~ "Fluoroquinolone Resistance", 
                               pretty_name == "spore viability" ~ "Spore viability", 
                               pretty_name == "growth rate" ~ "Max growth rate", 
                               pretty_name == "# spores" ~ "Total spores",
                               pretty_name == "toxin activity" ~ "Cytotoxicity", 
                               pretty_name == "Severity" ~ "Severe infection"))

set.seed(2) # to get a reproducible jitter distribution where all observed values aren't covered by negative controls
pyseer_plot <- best_results_df %>% 
  filter(bin_cont == "continuous") %>% 
  ggplot(mapping = aes(y = Rsquared,
                       x = fct_reorder(plot_name, plot_order), 
                       col = negative_control)) + 
  geom_jitter(width = 0.25, height = 0) +  
  scale_color_manual(values = c("grey", "black")) + 
  theme_bw() + 
  xlab("") + 
  ylab(expression(italic("R"^2))) + 
  theme(legend.position = "none") + 
  labs(tag = "B") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dot_plot_wo_leg + pyseer_plot + 
  plot_layout(ncol = 2, widths = c(7, 7)) & 
  theme(plot.tag = element_text(face = "bold"))
ggsave("../../figures/Fig2.png", width = 8, height = 4, units = "in")
ggsave("../../figures/Fig2.pdf", width = 8, height = 4, units = "in")
