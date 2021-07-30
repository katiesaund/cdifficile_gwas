library(tidyverse)
snp_annotations <- read_tsv("../../data/1_combine_genotype_matrices/snpmat/snp_annots.tsv")
save(snp_annotations, file = "../../data/1_combine_genotype_matrices/snpmat/snp_annots.RData")


indel_annotations <- read_tsv("../../data/1_combine_genotype_matrices/indelmat/indel_annots.tsv")
save(indel_annotations, file = "../../data/1_combine_genotype_matrices/indelmat/indel_annots.RData")
