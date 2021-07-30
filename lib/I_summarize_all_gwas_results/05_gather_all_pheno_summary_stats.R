# Goal round up summary stats for all of the phenotypes into one dataframe. 

# Phylogenetic Signal - already collected
# Number of significant p-values
# Number of high epsilon p-values
# Number of & type of variants: included as inputs & included as outputs (difference = number filtered out)
# Number of & type groups used in grouped analysis. 

source("gwas_analysis_lib.R")
data_dir <- "../../data/"
phylo_signal_df <- read_tsv(paste0(data_dir, "4_phenotypes/phenotypic_analysis_results/phylo_signal_with_neg_ctrls.tsv"), col_names = TRUE)
thresholds_df <- read_csv(paste0(data_dir, "7_thresholds/phenotype_thresholds.csv"), col_names = TRUE)
indiv_variant_df <- read_tsv(paste0(data_dir, "12_variant_types/variant_type_indiv.tsv"), col_names = TRUE)
group_variant_df <- read_tsv(paste0(data_dir, "12_variant_types/variant_type_group.tsv"), col_names = TRUE)
variant_df <- rbind(indiv_variant_df, group_variant_df)
phylo_tresh_df <- left_join(phylo_signal_df, thresholds_df, by = c("shorthand", "longhand", "bin_cont"))
num_pheno <- nrow(phylo_tresh_df)

summary_df <- left_join(phylo_tresh_df, variant_df, by = c("shorthand", "bin_cont"))
write_tsv(summary_df, path = paste0(data_dir, "13_summaries/hogwash_summary_data_by_phenotype.tsv"), col_names = TRUE)
