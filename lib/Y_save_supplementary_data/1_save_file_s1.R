# Need to save file S1 (cytoxicity gwas results) as two files do to size 
# constraints
library(tidyverse)
df <- read_tsv("../../data/13_summaries/log_toxin_indiv_genotypes_ranked_by_epsilon_then_p.tsv")
dim(df)
halfway_num <- (nrow(df)/2)

# let's split it in half: S1a and S1b
dfa <- df %>% slice(1:halfway_num)
dfb <- df %>% slice((halfway_num + 1):nrow(df))
stopifnot(nrow(dfa) + nrow(dfb) == nrow(df))
stopifnot(sum(dfa$Genotype %in% dfb$Genotype) == 0)
stopifnot(sum(dfb$Genotype %in% dfa$Genotype) == 0)
dfa[1:5, ]
dfb[1:5, ]
dfb[9001:9005, ]
write_tsv(dfa,
          "../../data/13_summaries/supplementary_file_1a.tsv", 
          col_names = TRUE)
write_tsv(dfb, "../../data/13_summaries/supplementary_file_1b.tsv", 
col_names = TRUE)
