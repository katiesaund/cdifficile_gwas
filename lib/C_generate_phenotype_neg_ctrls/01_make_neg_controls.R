# Goal: 
# We'll create 3 negative controls for each specific phenotype: 
# A negative control is the phenotype randomly scrambled on the tree. 

library(tidyverse)
ctrl_dir_path <- "../../data/4_phenotypes/negative_controls/"
data_dir_path <-  "../../data/15_pyseer/phenotypes_for_pyseer/"
pheno_paths <- paste0(data_dir_path, 
                      c("cleaned_fqR.tsv", 
                        "cleaned_log_cfe.tsv",
                        "cleaned_log_germ_tc_and_gly.tsv", 
                        "cleaned_log_germ_tc.tsv",
                        "cleaned_log_growth.tsv",
                        "cleaned_log_sporulation.tsv",
                        "cleaned_log_toxin.tsv",
                        "cleaned_severity.tsv"))
names(pheno_paths) <- c("fqR", "log_cfe", "log_germ_tc_and_gly", "log_germ_tc", 
                        "log_growth", "log_sporulation", "log_toxin", "severity")
num_pheno <- length(pheno_paths)

for (i in 1:num_pheno) {
  current_df <- read_tsv(pheno_paths[i], col_names = TRUE)
  num_sample <- nrow(current_df)
  current_name <- names(pheno_paths)[i]
  set.seed(1)
  jumbled_pheno_df <- cbind(current_df, 
                            sample(x = unname(unlist(current_df[, 2])),
                                          size = num_sample,
                                          replace = FALSE), # jumbled 1
                            sample(x = unname(unlist(current_df[, 2])),
                                          size = num_sample, 
                                          replace = FALSE), # jumbled 2
                            sample(x = unname(unlist(current_df[, 2])),
                                          size = num_sample, 
                                          replace = FALSE)) # jumbled 3
  colnames(jumbled_pheno_df) <- c("sample", 
                                  "true_pheno",
                                  "neg_ctrl_1_jumble", 
                                  "neg_ctrl_2_jumble",
                                  "neg_ctrl_3_jumble")
  
  
  write_tsv(jumbled_pheno_df %>% select(sample, neg_ctrl_1_jumble),
            path = paste0(ctrl_dir_path, "neg_ctrl_", current_name, "_1.tsv"), 
            col_names = TRUE)
  write_tsv(jumbled_pheno_df %>% select(sample, neg_ctrl_2_jumble),
            path = paste0(ctrl_dir_path, "neg_ctrl_", current_name, "_2.tsv"), 
            col_names = TRUE)
  write_tsv(jumbled_pheno_df %>% select(sample, neg_ctrl_3_jumble),
            path = paste0(ctrl_dir_path, "neg_ctrl_", current_name, "_3.tsv"), 
            col_names = TRUE)
}

