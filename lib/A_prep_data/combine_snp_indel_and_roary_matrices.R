# Goals: 
#  1. Combine rereferenced SNP and Indel matrices with the roary pan genome. 
#  2. For SNPs keep only those annotated by SnpEFF to be of moderate, high, or 
#     modifier predicted functional impact
source("data_prep_lib.R")

# Read in data -------------------------------------------------------------
load("../../data/1_combine_genotype_matrices/indelmat/indel_ancestral_referenced.RData")
indel_mat <- indel_reref_out$bin$mat
indel_annot_mat <- indel_reref_out$bin$annots

load("../../data/1_combine_genotype_matrices/snpmat/snp_ancestral_referenced.RData")
snp_mat <- snp_reref_out$bin$mat
snp_annot_mat <- snp_reref_out$bin$annots

pangenome_mat <- 
  read.table("../../data/1_combine_genotype_matrices/pangenomemat/cleaned_names_pangenomemat.tsv",
             header = TRUE, 
             row.names = 1, 
             sep = "\t", 
             quote = "", 
             stringsAsFactors = FALSE)

if (ncol(indel_mat) != ncol(snp_mat) | ncol(snp_mat) != ncol(pangenome_mat)) {
  stop("Col number mismatch")
}

if (!identical(nrow(indel_mat), nrow(indel_annot_mat)) | 
    !identical(nrow(snp_mat), nrow(snp_annot_mat))) {
  stop("Size mismatch")
}

tree <- read.tree("../../data/3_trees/midpoint_rooted.tree")

# Subset SNPs to only modifier, moderate, or high snpeff annotations -----------    
snp_rows_to_keep <- snp_annot_mat$snpeff_impact %in% c("MODERATE", "HIGH", "MODIFIER")
snp_annot_mat <- snp_annot_mat[snp_rows_to_keep, ]
snp_mat <- snp_mat[snp_rows_to_keep, ]

# No need to subset indels, because they're by definition not LOW impact

# Add rownames from binary matrices to annotation matrix
if (!identical(nrow(indel_mat), nrow(indel_annot_mat)) | 
    !identical(nrow(snp_mat), nrow(snp_annot_mat))) {
  stop("Size mismatch")
}

snp_annot_mat$bin_row_name <- row.names(snp_mat)
indel_annot_mat$bin_row_name <- row.names(indel_mat)

# The trees were rooted and tip labels correctly ordered in 0_clean_up_tree.R

# Check that all of the names agree --------------------------------------------
sum(!tree$tip.label %in% colnames(snp_mat)) # 0
sum(!tree$tip.label %in% colnames(indel_mat)) # 0
sum(!tree$tip.label %in% colnames(pangenome_mat)) # 0

# Check that there aren't any NAs in the matrices ------------------------------
if (sum(sum(is.na(snp_mat)), 
        sum(is.na(indel_mat)), 
        sum(is.na(pangenome_mat))) != 0) {
  stop("Some NAs in one of the matrices")
}

# Reorder matrices by tree tips ------------------------------------------------
pangenome_mat <- 
  pangenome_mat[, match(tree$tip.label, colnames(pangenome_mat)), drop = FALSE]
snp_mat <- 
  snp_mat[, match(tree$tip.label, colnames(snp_mat)), drop = FALSE]
indel_mat <- 
  indel_mat[, match(tree$tip.label, colnames(indel_mat)), drop = FALSE]

if (sum(sum(tree$tip.label != colnames(pangenome_mat)), 
        sum(tree$tip.label != colnames(snp_mat)), 
        sum(tree$tip.label != colnames(indel_mat))) != 0) {
  stop("Matrices still out of order")
}

# Combine snp & indel & pan ----------------------------------------------------
joint <- rbind(pangenome_mat, snp_mat, indel_mat)

# QC ---------------------------------------------------------------------------
dim(joint) # 99499   107

# Check that the matrix is binary
range(joint) # [1] 0 1

# How much variation occurs in each variant (row)?
range(rowSums(joint)) # 0 107

# Drop rows without any variation
joint <- joint[rowSums(joint) > 0 & rowSums(joint) < ncol(joint), , drop = FALSE]

# How much variation occurs in each variant (row)?
range(rowSums(joint)) # 1 106
dim(joint) # 85487   107

# Get the samples in the rows and genotype in the columns as per hogwash 
# requirements
joint <- t(joint)

dim(joint) # 107 85487

# Subset annotations to only keep the genotypes in this final version of the joint matrix
genotypes_in_joint <- colnames(joint)

dim(indel_annot_mat) #  8471   33
dim(snp_annot_mat) #  80859    31
indel_annot_mat <- indel_annot_mat[indel_annot_mat$bin_row_name %in% genotypes_in_joint, , drop = FALSE]
snp_annot_mat <- snp_annot_mat[snp_annot_mat$bin_row_name %in% genotypes_in_joint, , drop = FALSE]
dim(indel_annot_mat) # 8351   33
dim(snp_annot_mat) # 769592    31
# Save the outputs as tsvs -----------------------------------------------------
write.table(x = joint, 
            quote = TRUE,
            file = "../../data/1_combine_genotype_matrices/combinedmat/pan_and_rereferenced_snp_indel.tsv",
            sep = "\t", 
            col.names = TRUE, 
            row.names = TRUE)

save(x = indel_annot_mat, 
     file = "../../data/1_combine_genotype_matrices/combinedmat/indep_annot.RData")

save(x = snp_annot_mat, 
     file = "../../data/1_combine_genotype_matrices/combinedmat/snp_annot.RData")
