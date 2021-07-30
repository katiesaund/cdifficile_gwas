# Goal: 
#   Write SNP and indel to genes/intergenic regions keys.
source("data_prep_lib.R")

# Read in data -------------------------------------------------------------
load("../../data/1_combine_genotype_matrices/combinedmat/indep_annot.RData")
#indel_annot_mat

load("../../data/1_combine_genotype_matrices/combinedmat/snp_annot.RData")
#snp_annot_mat

# Gene Key -----------------------------------------------------------------
snp_annot_gene_subset_mat <- snp_annot_mat %>% 
  select(-raw_rownames) %>% 
  filter(intergenic == FALSE, snpeff_impact != "MODIFIER") %>% 
  select(bin_row_name, locus_tag)
indel_annot_gene_subset_mat <- indel_annot_mat %>% 
  select(-raw_rownames) %>% 
  filter(intergenic == FALSE, snpeff_impact != "MODIFIER") %>% 
  select(bin_row_name, locus_tag)
gene_key_mat <- 
  rbind(snp_annot_gene_subset_mat,
        indel_annot_gene_subset_mat)
colnames(gene_key_mat) <- c("Locus", "Group")

# Intergenic Key -----------------------------------------------------------
snp_annot_intergenic_subset_mat <- snp_annot_mat %>% 
  select(-raw_rownames) %>% 
  filter(intergenic == TRUE, snpeff_impact == "MODIFIER") %>% 
  select(bin_row_name, locus_tag)
indel_annot_intergenic_subset_mat <- indel_annot_mat %>% 
  select(-raw_rownames) %>% 
  filter(intergenic == TRUE, snpeff_impact == "MODIFIER") %>% 
  select(bin_row_name, locus_tag)
intergenic_key_mat <- 
  rbind(snp_annot_intergenic_subset_mat,
        indel_annot_intergenic_subset_mat)
colnames(intergenic_key_mat) <- c("Locus", "Group")

# Joint Key ----------------------------------------------------------------
gene_intergenic_key_mat <- rbind(gene_key_mat, intergenic_key_mat)

# Note we lose ~2000 genotypes this way: snp_annot_mat %>% filter(!bin_row_name %in% snp_annot_intergenic_subset_mat$bin_row_name & !bin_row_name %in% snp_annot_gene_subset_mat$bin_row_name) %>% View()
# These seem to be in exons (noncoding snps that are intragenic)

write.table(gene_key_mat, 
            file = "../../data/2_gene_keys/gene_key.tsv",
            quote = TRUE,
            sep = "\t", 
            col.names = TRUE, 
            row.names = TRUE)
write.table(intergenic_key_mat, 
            file = "../../data/2_gene_keys/intergenic_key.tsv", 
            quote = TRUE,
            sep = "\t", 
            col.names = TRUE, 
            row.names = TRUE)
write.table(gene_intergenic_key_mat, 
            file = "../../data/2_gene_keys/gene_and_intergenic_key.tsv", 
            quote = TRUE,
            sep = "\t", 
            col.names = TRUE, 
            row.names = TRUE)
