# Goal -- compare model performance

# Read in model results
library(tidyverse)
phenotypes <- read_tsv("../../../data/4_phenotypes/pheno_names.tsv")

results_df <- read_tsv("../../../data/15_pyseer/tcdb_model_results//EN_model_summary.tsv", col_names = TRUE)

results_df <- results_df %>% 
  filter(!rsq %in% c("error","TBD")) # for now, these are jobs that ran out of memory/time or are only the indel matrix

phenos <- c("fqR", "log_cfe", "log_germ_tc_and_gly", 
            "log_germ_tc", "log_growth", "log_sporulation", 
            "log_toxin", "severity")


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
ggsave(paste0("../../../figures/pyseer_results/tcdb_EN_model_performance_jitter_all_pheno.png"), 
       units = "in", height = 8.5, width = 7)

# Which alpha is the best for the majority of phenotypes?
best_alpha_df <- results_df %>% 
  filter(negative_control == FALSE) %>% 
  group_by(phenotype) %>% 
  filter(rsq == max(rsq)) %>%
  select(alpha, phenotype, rsq)
write_tsv(best_alpha_df, "../../../data/15_pyseer/tcdb_model_results/best_alpha_values.tsv")
# We observe that the best alpha is 0.5

# Let's plot just when the Rsq is the best - specific value for each phenotype
best_results_df <- results_df %>%  filter(alpha == 0.5) 
      
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

ggsave(paste0("../../../figures/pyseer_results/tcdb_EN_model_performance_jitter_best_alpha_all_pheno.png"), 
       units = "in", height = 8.5, width = 7)

write_tsv(best_results_df, "../../../data/15_pyseer/tcdb_model_results/best_model_results.tsv", col_names = TRUE)

