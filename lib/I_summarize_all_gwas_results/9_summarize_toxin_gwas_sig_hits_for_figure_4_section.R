library(tidyverse)
data_dir <- "../../data/"

thresholds <- read_csv("../../data/7_thresholds/phenotype_thresholds.csv", col_names = TRUE)
toxin_df <- read_tsv("../../data/13_summaries/log_toxin_indiv_genotypes_ranked_by_epsilon_then_p.tsv")
pval_thresh <- thresholds %>% filter(shorthand == "log_toxin") %>% pull(indiv_pval)
ep_thresh <- thresholds %>% filter(shorthand == "log_toxin") %>% pull(indiv_epsilon)
hogwash_all_df <- read_tsv("../../data/13_summaries/hogwash_summary_data_by_phenotype.tsv")
hogwash_all_df %>% 
  filter(`Input or output` == "sig. output") %>% 
  filter(shorthand != "fqR") %>%  
  filter(`Indiv. or Group.` == "Individual") %>% 
  select(shorthand, longhand, `Variant type`, n, test) %>% 
  View()

num_toxin_sig_hits <- toxin_df %>% 
  filter(`-ln(FDR P-value)` > -log(pval_thresh)) %>% 
  nrow()

num_toxin_hi_ep <- toxin_df %>% 
  filter(Epsilon > ep_thresh) %>% 
  nrow()

num_toxin_sig_hits <- toxin_df %>% 
  filter(`-ln(FDR P-value)` > -log(pval_thresh)) %>% 
  nrow()

num_toxin_hi_ep_and_sig_hits <- toxin_df %>% 
  filter(Epsilon > ep_thresh) %>% 
  filter(`-ln(FDR P-value)` > -log(pval_thresh)) %>% 
  nrow()


summary_df <- matrix(NA, nrow = 5, ncol = 2)
summary_df[1, ] <- c("toxin high epsilon", num_toxin_hi_ep)
summary_df[2, ] <- c("toxin sig p-value", num_toxin_sig_hits)
summary_df[3, ] <- c("toxin high epsilon and sig p-value", num_toxin_hi_ep_and_sig_hits)
summary_df <- as.data.frame(summary_df)

write_tsv(summary_df, 
          path = paste0(data_dir, "13_summaries/hogwash_summary_toxin_invid_sig.tsv"), 
          col_names = TRUE)

