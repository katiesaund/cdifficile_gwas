# Filter and reference the indels. 
# Make sure to always use the original Hanna Lab IDs -- no more using the 
# sequencing IDs (sequencing IDs started with PH, Hanna Lab IDs start with DA, 
# Cd, or FOBT)

source("data_prep_lib.R")

# Inputs ---------------------------
# Read in trees
tree <- read.tree("../../data/3_trees/midpoint_rooted.tree")

# Read in indel matrix
# Code
indelmat_code <- read.table("../../data/1_combine_genotype_matrices/indelmat/Indel_matrix_code.tsv",
                          sep = "\t", 
                          header = TRUE, 
                          row.names = 1, 
                          quote = "",
                          stringsAsFactors = FALSE)
# Allele
indelmat_allele <- read.table("../../data/1_combine_genotype_matrices/indelmat/Indel_matrix_allele.tsv",
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
indelmat_code <- clean_up_mat_names(indelmat_code, key)  
indelmat_allele <- clean_up_mat_names(indelmat_allele, key)

indelmat_code <- indelmat_code[, colnames(indelmat_code) %in% key$Hanna_ID, drop = FALSE]
indelmat_allele <- indelmat_allele[, colnames(indelmat_allele) %in% key$Hanna_ID, drop = FALSE]

# Remove samples in indelmats that are not included in trees
indelmat_code <- subset_mat_to_tips(tree, indelmat_code)
indelmat_allele <- subset_mat_to_tips(tree, indelmat_allele)

# Now reference the variants
# Rereference variants based on ancestral state at root of trees
indel_reref_out <- 
  snitkitr::parse_indels(varmat_code = indelmat_code,
                       varmat_allele = indelmat_allele, 
                       tree = tree, 
                       remove_multi_annots = FALSE,
                       return_binary_matrix = TRUE, 
                       ref_to_anc = TRUE)

write.table(indel_reref_out$bin$mat, 
            file = "../../data/1_combine_genotype_matrices/indelmat/ancestral_referenced_indelmat.tsv",
            sep = "\t", 
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE)

write.table(indel_reref_out$bin$annots, 
            file = "../../data/1_combine_genotype_matrices/indelmat/indel_annots.tsv",
            sep = "\t", 
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE)

save(indel_reref_out, file = "../../data/1_combine_genotype_matrices/indelmat/indel_ancestral_referenced.RData")
