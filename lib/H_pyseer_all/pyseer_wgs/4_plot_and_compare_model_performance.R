# Goal -- compare model performance

# Read in model results
library(tidyverse)
phenotypes <- read_tsv("../../../data/4_phenotypes/pheno_names.tsv")

results_df <- read_tsv("../../../data/15_pyseer/wgs_model_results/EN_model_summary.tsv", col_names = TRUE)

results_df <- results_df %>% 
  filter(!rsq %in% c("error","TBD")) # for now, these are jobs that ran out of memory/time or are only the indel matrix

phenos <- c("fqR", "log_cfe", "log_germ_tc_and_gly", 
              "log_germ_tc", "log_growth", "log_sporulation", 
              "log_toxin", "severity")

for (i in 1:length(phenos)) {
  current_pheno <- phenos[i]

  results_df %>% 
    filter(grepl(pattern = current_pheno, x = phenotype)) %>% 
    mutate(true_pheno_or_neg_ctrl = case_when(grepl("neg_ctrl_", phenotype) ~ "neg_ctrl_pheno", 
                                              !grepl("neg_ctrl_", phenotype) ~ current_pheno), 
           Rsquared = as.numeric(as.character(rsq))) %>% 
    ggplot(mapping = aes(y = Rsquared,
                         x = n_fold, 
                         col = alpha)) + 
    geom_jitter(width = 0.05, height = 0) +  
    facet_grid(genotype ~ true_pheno_or_neg_ctrl + lineage) + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = 0, col = "grey") + 
    scale_color_gradient(low = "blue", high = "red")
  
  ggsave(paste0("../../../figures/pyseer_results/wgs_EN_model_performance_jitter_",  current_pheno, ".png"), 
         units = "in", height = 8.5, width = 7)
}

# Now plot all together 
# Need to get phenotype to be just one of the 8 and add a column for neg ctrl

results_df <- results_df %>% 
  mutate(negative_control = grepl("neg", phenotype)) %>% 
  mutate(phenotype = gsub("neg_ctrl_", "", phenotype) %>% gsub("_[123]", "", .)) %>% 
  mutate( Rsquared = as.numeric(as.character(rsq)))


results_df %>% 
  ggplot(mapping = aes(y = Rsquared,
                       x = phenotype, 
                       shape = negative_control, 
                       col = alpha)) + 
  geom_jitter(width = 0.05, height = 0) +  
  geom_hline(yintercept = 0, col = "grey") + 
  scale_color_gradient(low = "blue", high = "red")
ggsave(paste0("../../../figures/pyseer_results/wgs_EN_model_performance_jitter_all_pheno.png"), 
       units = "in", height = 8.5, width = 7)

# Which alpha is the best for the majority of phenotypes?
best_alpha_df <- results_df %>% 
  filter(negative_control == FALSE) %>% 
  group_by(phenotype) %>% 
  filter(rsq == max(rsq)) %>%
  select(alpha, phenotype, rsq)
write_tsv(best_alpha_df, "../../../data/15_pyseer/wgs_model_results/best_alpha_values.tsv")
# We observe that:
# fqR is 1 for all 5 alpha values
# alpha phenotype             rsq
# <dbl> <chr>               <dbl>
#   1 0.01  log_toxin           0.899
# 2 0.01  severity            1    
# 3 0.01  log_germ_tc_and_gly 0.029
# 4 0.01  log_germ_tc         0.997
# 5 0.255 fqR                 1    
# 6 0.255 log_cfe             0.407
# 7 0.255 log_growth          0.4  
# 8 0.255 log_sporulation     0.064
# 9 0.5   fqR                 1    
# 10 0.745 fqR                 1    
# 11 0.99  fqR                 1    

# Let's plot just when the Rsq is the best - specific value for each phenotype

alpha_01_df <- results_df %>% 
  filter(alpha == 0.01, 
         phenotype %in% c("fqr",
                          "log_toxin",
                          "severity", 
                          "log_germ_tc_and_gly",
                          "log_germ_tc"))
alpha_255_df <- results_df %>% 
  filter(alpha == 0.255, 
         phenotype %in% c("log_cfe", 
                          "log_growth", 
                          "log_sporulation"))

best_results_df <- rbind(alpha_01_df, alpha_255_df)

# add long hand name for plotting
best_results_df <- left_join(best_results_df, phenotypes, by = c("phenotype" =  "shorthand"))
best_results_df$negative_control[best_results_df$negative_control == TRUE] <- "Negative control"
best_results_df$negative_control[best_results_df$negative_control == FALSE] <- "Observed value"

set.seed(4) # to get a reproducible jitter distribution where all observed values aren't covered by negative controls
best_results_df %>% 
  ggplot(mapping = aes(y = Rsquared,
                       x = longhand, 
                       col = negative_control)) + 
  geom_jitter(width = 0.25, height = 0) +  
  scale_color_manual(values = c("grey", "black")) + 
  theme_bw()+ 
  xlab("") + 
  ylab(expression("R"^2)) + 
  theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("../../../figures/pyseer_results/wgs_EN_model_performance_jitter_best_alpha_all_pheno.png"), 
       units = "in", height = 8.5, width = 7)

set.seed(3) # to get a reproducible jitter distribution where all observed values aren't covered by negative controls
best_results_df %>% 
  filter(bin_cont == "continuous") %>% 
  ggplot(mapping = aes(y = Rsquared,
                       x = longhand, 
                       col = negative_control)) + 
  geom_jitter(width = 0.25, height = 0) +  
  scale_color_manual(values = c("grey", "black")) + 
  theme_bw() + 
  xlab("") + 
  ylab(expression(italic("R"^2))) + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.2, 0.9), 
        legend.background = element_rect(fill = "white", colour = "black")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0("../../../figures/pyseer_results/wgs_EN_model_performance_jitter_best_alpha_continuous_pheno.png"), 
       units = "in", height = 8.5, width = 7)

write_tsv(best_results_df, "../../../data/15_pyseer/wgs_model_results/best_model_results.tsv", col_names = TRUE)

