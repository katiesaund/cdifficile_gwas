# Filter and reference the SNPs. 
# Make sure to always use the original Hanna Lab IDs -- no more using the 
# sequencing IDs (sequencing IDs started with PH, Hanna Lab IDs start with DA, 
# Cd, or FOBT)

source("data_prep_lib.R")

# Inputs ---------------------------
# Read in trees
tree <- read.tree("../../data/3_trees/midpoint_rooted.tree")

# Read in SNP matrix
# Code
snpmat_code <- read.table("../../data/1_combine_genotype_matrices/snpmat/SNP_matrix_code.tsv",
                          sep = "\t", 
                          header = TRUE, 
                          row.names = 1, 
                          quote = "",
                          stringsAsFactors = FALSE)
# Allele
snpmat_allele <- read.table("../../data/1_combine_genotype_matrices/snpmat/SNP_matrix_allele_new.tsv",
                            sep = "\t", 
                            header = TRUE, 
                            row.names = 1, 
                            quote = "",
                            stringsAsFactors = FALSE)

# Read in name key
key <- read.table("../../data/0_sample_list/genomes_to_use_in_analysis.tsv",
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE)

# Clean up names in matrices and tree ------ 
snpmat_code <- clean_up_mat_names(snpmat_code, key)  
snpmat_allele <- clean_up_mat_names(snpmat_allele, key)

snpmat_code <- snpmat_code[, colnames(snpmat_code) %in% key$Hanna_ID, drop = FALSE]
snpmat_allele <- snpmat_allele[, colnames(snpmat_allele) %in% key$Hanna_ID, drop = FALSE]

# Remove samples in snpmats that are not included in trees
snpmat_code <- subset_mat_to_tips(tree, snpmat_code)
snpmat_allele <- subset_mat_to_tips(tree, snpmat_allele)

# Now reference the variants
# Rereference variants based on ancestral state at root of trees
snp_reref_out <- 
  snitkitr::parse_snps(varmat_code = snpmat_code,
                       varmat_allele = snpmat_allele, 
                       tree = tree, 
                       remove_multi_annots = FALSE,
                       return_binary_matrix = TRUE, 
                       ref_to_anc = TRUE)

write.table(snp_reref_out$bin$mat, 
            file = "../../data/1_combine_genotype_matrices/snpmat/ancestral_referenced_snpmat.tsv",
            sep = "\t", 
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE)

write.table(snp_reref_out$bin$annots, 
            file = "../../data/1_combine_genotype_matrices/snpmat/snp_annots.tsv",
            sep = "\t", 
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE)

save(snp_reref_out, file = "../../data/1_combine_genotype_matrices/snpmat/snp_ancestral_referenced.RData")
