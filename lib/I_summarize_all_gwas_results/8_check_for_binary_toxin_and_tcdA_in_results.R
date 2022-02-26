# Double checking for P-value and epsilon for binary toxin and toxin A in log 
# toxin GWAS results.
library(tidyverse)

toxin_tbl <- read_tsv("../../data/13_summaries/log_toxin_indiv_genotypes_ranked_by_epsilon_then_p.tsv", col_names = TRUE)

toxin_tbl %>% 
  filter(grepl("CdtAB", Genotype) | grepl("cdtR", Genotype) ) %>% 
  View()



toxin_tbl %>% 
  filter(grepl("Cdt", Genotype) | grepl("cdt", Genotype) ) %>% 
  select(Genotype) %>% 
  mutate(gene_info = sub(".*locus_tag=", "", Genotype)) %>% 
  pull(gene_info) 

# actually this isn't a good approach. Only a cdtAB pseudogene and cdtR is in 
# the core. cdtA/B is in the accessory (duh). 

toxin_tbl %>% 
  filter(grepl("tcdA", Genotype)) %>% 
  View()


toxin_tbl %>% 
  filter(grepl("tcdA", Genotype)) %>% 
  slice(1:1) %>%
  pull(Genotype)

toxin_tbl %>% 
  filter(grepl("[|]cdtR", Genotype)) %>% 
  slice(1:1) %>%
  pull(Genotype)

