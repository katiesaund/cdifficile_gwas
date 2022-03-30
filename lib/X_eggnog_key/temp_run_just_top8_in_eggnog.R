library(tidyverse)
top8_tbl <- read_tsv("../../data/18_eggnog/table_3_cytotoxicity_with_eggnog_key.tsv")
temp <- paste(top8_tbl$roary_prokka_ID, top8_tbl$Genotype)
temp <- temp[2:8]

pangenome <- read.FASTA("../../data/1_combine_genotype_matrices/pangenomemat/pan_genome_reference.fa")
subset <- pangenome[temp]
write.FASTA(subset, file = "../../data/18_eggnog/temp_top8_toxin.fa")
