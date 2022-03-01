library(tidyverse)
data_dir <- "../../data/"

thresholds <- read_csv("../../data/7_thresholds/phenotype_thresholds.csv", col_names = TRUE)

cfe_df <- read_tsv("../../data/13_summaries/log_cfe_indiv_genotypes_ranked_by_epsilon_then_p.tsv")
germ_tc_gly_df <- read_tsv("../../data/13_summaries/log_germ_tc_and_gly_indiv_genotypes_ranked_by_epsilon_then_p.tsv")
germ_tc_df <- read_tsv("../../data/13_summaries/log_germ_tc_indiv_genotypes_ranked_by_epsilon_then_p.tsv")
grow_df <- read_tsv("../../data/13_summaries/log_growth_indiv_genotypes_ranked_by_epsilon_then_p.tsv")
spor_df <- read_tsv("../../data/13_summaries/log_sporulation_indiv_genotypes_ranked_by_epsilon_then_p.tsv")
toxin_df <- read_tsv("../../data/13_summaries/log_toxin_indiv_genotypes_ranked_by_epsilon_then_p.tsv")
sev_phyc_df <- read_tsv("../../data/13_summaries/severity_phyc_indiv_genotypes_ranked_by_epsilon_then_p.tsv")
sev_sync_df <- read_tsv("../../data/13_summaries/severity_synchronous_indiv_genotypes_ranked_by_epsilon_then_p.tsv")

pval_thresh <- thresholds %>% filter(shorthand == "log_toxin") %>% pull(indiv_pval)
ep_thresh <- thresholds %>% filter(shorthand == "log_toxin") %>% pull(indiv_epsilon)
hogwash_all_df <- read_tsv("../../data/13_summaries/hogwash_summary_data_by_phenotype.tsv")

# number of significant and high epsilon variants per phenotype
small_hogwash_summary <- hogwash_all_df %>% 
  filter(`Input or output` == "sig. output") %>% 
  filter(shorthand != "fqR") %>%  
  filter(`Indiv. or Group.` == "Individual") %>% 
  select(shorthand, longhand, `Variant type`, n, Test) 

# Just toxin below

summarize_hits <- function(df) {
  num_sig_hits <- df %>% 
    filter(`-ln(FDR P-value)` > -log(pval_thresh)) %>% 
    nrow()
  
  num_hi_ep <- df %>% 
    filter(Epsilon > ep_thresh) %>% 
    nrow()
  
  num_hi_ep_and_sig_hits <- df %>% 
    filter(Epsilon > ep_thresh) %>% 
    filter(`-ln(FDR P-value)` > -log(pval_thresh)) %>% 
    nrow()
  return(c(num_hi_ep_and_sig_hits, num_sig_hits, num_hi_ep))
}


summary_df <- matrix(NA, nrow = 3, ncol = 9)
colnames(summary_df) <- c(thresholds$longhand[1:7], "severity: phyc", "severity:sync")
rownames(summary_df) <- c("high epsilon and sig p-value", 
                          "sig p-value", 
                          "high epsilon")
summary_df <- as.data.frame(summary_df)
summary_df$`ln(spore viability)` <- summarize_hits(cfe_df)
summary_df$`ln(germination in Tc and Gly)` <- summarize_hits(germ_tc_gly_df)
summary_df$`ln(germination in Tc)` <- summarize_hits(germ_tc_df)
summary_df$`ln(growth rate)` <- summarize_hits(grow_df)
summary_df$`ln(# spores)` <- summarize_hits(spor_df)
summary_df$`ln(toxin activity)` <- summarize_hits(toxin_df)
summary_df$`severity: phyc` <- summarize_hits(sev_phyc_df)
summary_df$`severity:sync` <- summarize_hits(sev_sync_df)

write_tsv(summary_df, 
          file = paste0(data_dir, "13_summaries/hogwash_summary_for_figure_4.tsv"), 
          col_names = TRUE)

pos_1967_locus <- toxin_df %>% 
  filter(grepl("1967", Genotype)) %>% 
  filter(`-ln(FDR P-value)` > -log(pval_thresh))

write_tsv(pos_1967_locus, 
          file = paste0(data_dir, "13_summaries/hogwash_summary_for_figure_4_tcdB_pos_1967_locus.tsv"), 
          col_names = TRUE)


tcdr_tcdb_locus <- toxin_df %>% 
  filter(grepl("CD630_06590-CD630_06600", Genotype)) %>% 
  filter(`-ln(FDR P-value)` > -log(pval_thresh)) 

write_tsv(tcdr_tcdb_locus, 
          file = paste0(data_dir, "13_summaries/hogwash_summary_for_figure_4_tcdR_tcdB_locus.tsv"), 
          col_names = TRUE)

tcdC_117 <- toxin_df %>% 
  filter(grepl("tcdC", Genotype)) %>% 
  filter(grepl("117", Genotype)) 

write_tsv(tcdC_117, 
          file = paste0(data_dir, "13_summaries/hogwash_summary_for_figure_4_tcdC_117.tsv"), 
          col_names = TRUE)

top8_toxin <- toxin_df %>% 
  filter(Epsilon > ep_thresh) %>% 
  filter(`-ln(FDR P-value)` > -log(pval_thresh)) %>% 
  select(Genotype, roary_prokka_ID, `-ln(FDR P-value)`, Epsilon, Rank) %>% 
  mutate(P = exp(-`-ln(FDR P-value)`)
  )

write_tsv(top8_toxin, 
          file = paste0(data_dir, "13_summaries/hogwash_summary_for_figure_4_toxin_table3_part1.tsv"), 
          col_names = TRUE)
