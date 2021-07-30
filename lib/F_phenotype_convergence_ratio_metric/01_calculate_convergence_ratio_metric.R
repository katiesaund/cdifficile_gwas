# paper= https://www.sciencedirect.com/science/article/pii/S0022519308000209
# Journal of Theoretical Biology
# Volume 252, Issue 1, 7 May 2008, Pages 1-14
# Journal of Theoretical Biology
# Is convergence surprising? An examination of the frequency of convergence in simulated datasets
# C. TristanStayton

# "Another method, inspired by a classic paper by Raup and Gould (1974), involves taking the ratio of the patristic distance to the phenetic distance between all possible pairs of taxa and then averaging this ratio over the entire tree. If two taxa are highly convergent, then the ratio of their patristic distance to phenetic distance will be very large (they will be close in multivariate space but distantly related). In a dataset with a large amount of convergence, the average of all such ratios will be very large.
# Patristic and phenetic distances will depend on the data being used (tree length and traits). To compare values between trees, phenetic and patristic distances were divided by the maximum distances observed in each dataset in order to standardize for varying tree length and trait structure. This Ratio metric can only be used in a relative sense, to compare convergence between trees; patristic distance will depend on the type of branch lengths being used and there is no theoretical upper limit to this metric, as there is with the others. This metric measures convergence, as well as geometric parallelism (when patristic distance increases, but phenetic distance does not), but as instances of convergence will have a larger effect on this metric, it is also used in this exploratory analysis."

library(ape)
library(tidyverse)

# For each tree and phenotype calculate the average of the patristic distance / phenetic distance for each taxa pair
# The relative convergence of each phenotype is conveyed by comparing these averages. 

phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")
phenotypes <- phenotypes %>% filter(bin_cont == "continuous") 
phenotypes <- phenotypes %>% 
  mutate(mean_pheno_to_pat_distance = 99999)

ratio_vec_summary <- rep(list(), 6)


for (i in 1:nrow(phenotypes)) {
  print(phenotypes$shorthand[i])
  # Use ape to get distance matrix for each tree
  tree <- ape::read.tree(paste0("../../data/3_trees/", 
                                phenotypes$shorthand[i], 
                                ".tree"))
  pheno <- read.table(paste0("../../data/4_phenotypes/cleaned_", 
                      phenotypes$shorthand[i], 
                      ".tsv"))
  
  pairwise_edge_dist_mat <- ape::cophenetic.phylo(tree)
  max_pairwise_edge_dat <- max(pairwise_edge_dist_mat, na.rm = TRUE)
  pairwise_edge_dist_mat <- pairwise_edge_dist_mat / max_pairwise_edge_dat
  
  pheno_dist_mat <- matrix(NA, nrow = Ntip(tree), ncol = Ntip(tree))
  colnames(pheno_dist_mat) <- row.names(pheno_dist_mat) <- tree$tip.label
  # same layout: 
  sum(colnames(pairwise_edge_dist_mat) != colnames(pheno_dist_mat))
  sum(row.names(pairwise_edge_dist_mat) != row.names(pheno_dist_mat))
  # calculate phenetic distance 
  for (k in 1:Ntip(tree)) {
    for (j in 1:Ntip(tree)) {
      current_pheno_k <-  pheno[row.names(pheno) == tree$tip.label[k], 1]
      current_pheno_j <-  pheno[row.names(pheno) == tree$tip.label[j], 1]
      pheno_dist_mat[k, j] <- abs(current_pheno_k - current_pheno_j)
    }
  }
  
  pheno_dist_mat <- pheno_dist_mat / max(pheno_dist_mat, na.rm = TRUE)
  
  # calculate ratio for each pair
  ratio_mat <- pheno_dist_mat / pairwise_edge_dist_mat # This is the same as the paper convratio.m script
        # function y = convratio(x,z)
        # %convratio  - Measures convergence by calculating the ratio
        # %phenetic distance/phylogenetic distance for each pair of taxa, and taking
        # %an average of these ratios for the entire tree.  Convergence is found when
        # %phenetic distance is much less than phylogenetic distance.  Both phenetic
        # %distance and phylogenetic distance are scaled before the ratio is
        # %calculated.  Distances are scaled by dividing by the maximum
        # %distance among taxa.
  
  # subset so each ratio is only counted once and remove self-self pairs
  lower_tri_ratio_vec <- ratio_mat[lower.tri(ratio_mat)]
  ratio_vec_summary[[i]] <- lower_tri_ratio_vec
  names(ratio_vec_summary)[[i]] <- phenotypes$shorthand[i]
  ratio_avg <- mean(lower_tri_ratio_vec)
  # get average for each phenotype
  
  phenotypes$mean_pheno_to_pat_distance[i] <- ratio_avg
}

write_tsv(phenotypes, "../../data/4_phenotypes/phenotypic_analysis_results/phenotype_convergence.tsv")
