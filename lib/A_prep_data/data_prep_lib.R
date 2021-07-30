library(ape)
library(tidyverse)
library(snitkitr)
library(stringr)
library(data.table)

print("Session Info")
print(sessionInfo())

clean_up_tree_names <- function(tree, key) {
  tree$tip.label <- gsub("_", "", tree$tip.label)
  
  for (i in 1:Ntip(tree)) {
    for (j in 1:nrow(key)) {
      tree_tip_name <- tree$tip.label[i]
      sequencing_ID <- key$Unique_ID[j]
      Hanna_lab_ID <- key$Hanna_ID[j]
      if (tree_tip_name == sequencing_ID) {
        tree$tip.label[i] <- Hanna_lab_ID
      }
    }
  }
  return(tree)
}

clean_up_mat_names <- function(mat, key) {
  colnames(mat) <- 
    gsub("R1.*", "", colnames(mat)) %>%
    sub(".fastq.gz", "", . ) %>% 
    gsub("_S[0-9].*", "", .) %>% 
    gsub("_", "", .)
  
  for (i in 1:ncol(mat)) {
    for (j in 1:nrow(key)) {
      mat_sample_name <- colnames(mat)[i]
      sequencing_ID <- key$Unique_ID[j]
      Hanna_lab_ID <- key$Hanna_ID[j]
      if (mat_sample_name == sequencing_ID) {
        colnames(mat)[i] <- Hanna_lab_ID
      }
    }
  }
  return(mat)
}


remove_og_and_ref <- function(mat, outgroup, ref) {
  mat <- 
    mat[, !(colnames(mat) %in% c(outgroup, ref)), drop = FALSE]
  return(mat)
}

subset_mat_to_tips <- function(tr, mat){
  not_in_tree <- colnames(mat)[!colnames(mat) %in% tr$tip.label]
  mat <- mat[, !(colnames(mat) %in% not_in_tree), drop = FALSE]
  return(mat)
}

keep_good_samples <- function(pheno_mat, key_mat) {
  pheno_mat <- pheno_mat[row.names(pheno_mat) %in% key$Hanna_ID, , drop = FALSE]
  return(pheno_mat)
}

make_pheno_specific_tree <- function(tree, pheno_mat) {
  pheno_drop <- tree$tip.label[!tree$tip.label %in% row.names(pheno_mat)]
  pheno_tree <- tree
  if (!purrr::is_empty(pheno_drop)) {
    pheno_tree <- drop.tip(tree, pheno_drop)
  } 
  return(pheno_tree)
}

order_pheno_by_tree <- function(tree, pheno) {
  # Order phenotype by tree tips
  pheno_mat <- pheno[match(tree$tip.label, row.names(pheno)), , drop = FALSE]
  return(pheno_mat)
}

check_tree_and_pheno_match <- function(tree, pheno) {
  if (sum(tree$tip.label != row.names(pheno)) != 0) {
    stop("pheno doesn't match tree")
  }
}

save_tree_and_pheno <- function(tree, pheno_mat, name) {
  write.tree(tree, file = paste0("../../data/3_trees/", name, ".tree"))
  write.table(pheno_mat, 
              sep = "\t", 
              quote = FALSE,
              row.names = TRUE, 
              col.names = TRUE, 
              file = paste0("../../data/4_phenotypes/cleaned_", name, ".tsv"))
}

prep_phenotype <- function(tree, pheno_mat, key, name) {
  pheno_mat <- keep_good_samples(pheno_mat, key)
  tree <- make_pheno_specific_tree(tree, pheno_mat)
  pheno_mat <- order_pheno_by_tree(tree, pheno_mat)
  check_tree_and_pheno_match(tree, pheno_mat)
  save_tree_and_pheno(tree, pheno_mat, name)
}

save_pheno_specific_geno <- function(pheno_name, pheno_mat, geno_mat, key_mat) {
  pheno_samples <- row.names(pheno_mat)
  
  print(pheno_name)
  print(length(pheno_samples))
  print(dim(geno_mat))
  
  # Subset genotype matrix to have only the samples found in this specific phenotype
  keep_geno_mat <- geno_mat[row.names(geno_mat) %in% pheno_samples, , drop = FALSE]
  
  # Subset genotype matrix to keep only those genotypes found in at least 1 sample and less than all samples
  keep_geno_mat <- keep_geno_mat[, colSums(keep_geno_mat) < nrow(keep_geno_mat) & 
                                   colSums(keep_geno_mat) > 0]
  
  # Subset the genotype key to keep only those genotypes found in the subset genotype matrix
  key_mat <- key_mat[, 2:3] # drop row names column
  print(dim(key_mat))
  
  keep_key_mat <- key_mat[key_mat[, 1] %in% colnames(keep_geno_mat), , drop = FALSE]
  
  # Reorder genotype matrix to match samples as in the phenotype matrix
  keep_geno_mat <- 
    keep_geno_mat[match(pheno_samples, row.names(keep_geno_mat)), , drop = FALSE]
  # No need to reorder the key, because it doesn't have sample names in it. 
  
  print("after subseting and reordering")
  print(dim(keep_geno_mat))
  print(dim(keep_key_mat))
  
  write.table(x = keep_geno_mat, 
              quote = TRUE,
              file = paste0("../../data/1_combine_genotype_matrices/combinedmat/", pheno_name, "_pan_and_rereferenced_snp_indel.tsv"),
              sep = "\t", 
              col.names = TRUE, 
              row.names = TRUE)
  
  write.table(x = keep_key_mat, 
              quote = TRUE,
              file = paste0("../../data/2_gene_keys/", pheno_name, "_gene_and_intergenic_key.tsv"),
              sep = "\t", 
              col.names = TRUE, 
              row.names = TRUE)
  
}
