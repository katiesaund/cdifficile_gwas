# hogwash and pyseer require genotype matrices to be flipped relative to each other
# pyseer needs samples in the column names and variants in the row names

pheno_names <- c("fqR", "log_cfe", "log_germ_tc_and_gly", "log_germ_tc",
                 "log_growth", "log_sporulation", "log_toxin", "severity")

num_pheno <- length(pheno_names)

for (i in 1:num_pheno) {
  geno_mat_name <- paste0("../../../data/1_combine_genotype_matrices/combinedmat/",
                          "cleaned_",
                          pheno_names[i],
                          "_pan_and_rereferenced_snp_indel.tsv")
  
  genotype <- data.table::fread(file = geno_mat_name,
                                sep = "\t", 
                                data.table = FALSE, 
                                stringsAsFactors = FALSE, 
                                check.names = FALSE)
  genotype <- t(genotype)
  # genotype <- as.matrix(genotype)
  # row.names(genotype) <- genotype[, 1]
  # genotype <- genotype[, 2:ncol(genotype)]
  # storage.mode(genotype) <- "numeric"
  write.table(x = genotype, 
              quote = FALSE,
              file = paste0("../../../data/15_pyseer/genotypes_for_pyseer/", 
                            "pyseer_cleaned_",
                            pheno_names[i],
                            "_pan_and_rereferenced_snp_indel.tsv"), 
              sep = "\t", 
              col.names = FALSE, 
              row.names = TRUE)
}
