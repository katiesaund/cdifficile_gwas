# Calculate phylogenetic signal for negative controls
library(tidyverse)
library(hogwash)
library(ape)
options(scipen = 999) # turn off scientific notation


phylo_signal_df <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")
phylo_signal_df <- cbind(phylo_signal_df, rep(NA, 8), 
                         rep(NA, 8), rep(NA, 8), rep(NA, 8))
colnames(phylo_signal_df) <- c("shorthand", "longhand", "bin_cont",
                               "phylo_signal",
                               "neg_ctrl_phylo_signal_1",
                               "neg_ctrl_phylo_signal_2",
                               "neg_ctrl_phylo_signal_3")

for (i in 1:nrow(phylo_signal_df)) {
  current_pheno <- read.table(paste0("../../data/4_phenotypes/cleaned_", 
                                     phylo_signal_df$shorthand[i], 
                                     ".tsv"))
  current_pheno <- as.matrix(current_pheno)
  
  # neg ctrls:   neg_ctrl_log_germ_tc_3.tsv 
  neg_p_1 <- read_tsv(paste0("../../data/4_phenotypes//negative_controls/neg_ctrl_", 
                               phylo_signal_df$shorthand[i], 
                               "_1.tsv"))
  neg_p_1_mat <- as.matrix(neg_p_1$neg_ctrl_1_jumble)
  row.names(neg_p_1_mat) <- neg_p_1$sample
  neg_p_2 <- read_tsv(paste0("../../data/4_phenotypes/negative_controls/neg_ctrl_", 
                               phylo_signal_df$shorthand[i], 
                               "_2.tsv"))
  neg_p_2_mat <- as.matrix(neg_p_2$neg_ctrl_2_jumble)
  row.names(neg_p_2_mat) <- neg_p_2$sample
  neg_p_3 <- read_tsv(paste0("../../data/4_phenotypes/negative_controls/neg_ctrl_", 
                               phylo_signal_df$shorthand[i], 
                               "_3.tsv"))
  neg_p_3_mat <- as.matrix(neg_p_3$neg_ctrl_3_jumble)
  row.names(neg_p_3_mat) <- neg_p_3$sample
  
  
  current_tree <- read.tree(paste0("../../data/3_trees/", 
                                   phylo_signal_df$shorthand[i], 
                                   ".tree"))
  
  current_phylo_signal <- hogwash::report_phylogenetic_signal(current_pheno, current_tree)
  neg_1_phylo_signal <- hogwash::report_phylogenetic_signal(neg_p_1_mat, current_tree)
  neg_2_phylo_signal <- hogwash::report_phylogenetic_signal(neg_p_2_mat, current_tree)
  neg_3_phylo_signal <- hogwash::report_phylogenetic_signal(neg_p_3_mat, current_tree)
  
  if (phylo_signal_df$bin_cont[i] == "binary") {
    current_phylo_signal <- as.numeric(gsub(".*D = ", "", current_phylo_signal))
    neg_1_phylo_signal <- as.numeric(gsub(".*D = ", "", neg_1_phylo_signal))
    neg_2_phylo_signal <- as.numeric(gsub(".*D = ", "", neg_2_phylo_signal))
    neg_3_phylo_signal <- as.numeric(gsub(".*D = ", "", neg_3_phylo_signal))
  } else {
    current_phylo_signal <- as.numeric(gsub(".*a = ", "", current_phylo_signal))
    neg_1_phylo_signal <- as.numeric(gsub(".*a = ", "", neg_1_phylo_signal))
    neg_2_phylo_signal <- as.numeric(gsub(".*a = ", "", neg_2_phylo_signal))
    neg_3_phylo_signal <- as.numeric(gsub(".*a = ", "", neg_3_phylo_signal))
  }

  phylo_signal_df$phylo_signal[i] <- current_phylo_signal
  phylo_signal_df$neg_ctrl_phylo_signal_1[i] <- neg_1_phylo_signal
  phylo_signal_df$neg_ctrl_phylo_signal_2[i] <- neg_2_phylo_signal
  phylo_signal_df$neg_ctrl_phylo_signal_3[i] <- neg_3_phylo_signal
  
}

write_tsv(phylo_signal_df, "../../data/4_phenotypes/phenotypic_analysis_results/phylo_signal_with_neg_ctrls.tsv", col_names = TRUE)
