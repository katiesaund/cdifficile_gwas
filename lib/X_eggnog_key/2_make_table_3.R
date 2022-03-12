# make table 3: gwas cytoxicity top hits

eggnog_tbl <- read_tsv( "../../data/18_eggnog/file_s6_pangenome_eggnog_key.tsv", col_names = TRUE)
toxin_tbl <- read_tsv("../../data/13_summaries/log_toxin_indiv_genotypes_ranked_by_epsilon_then_p.tsv", col_names = TRUE)

top8_tbl <- toxin_tbl %>% 
  filter(Epsilon > 0.15) %>% 
  filter(`-ln(FDR P-value)` > -log(.15)) %>% 
  left_join(., eggnog_tbl, by = c("Genotype" = "Genome_ID_plus_Prokka_ID")) %>% 
  select(c(-reference_allele, -variant_allele, -amino_acid_change)) 

write_tsv(top8_tbl,
          file = paste0("../../data/18_eggnog/table_3_cytotoxicity_with_eggnog_key.tsv"))
