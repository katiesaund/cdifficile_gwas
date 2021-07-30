library(readr)
phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")
num_pheno <- nrow(phenotypes)
thresholds <- cbind(phenotypes, 
                    rep(NA, num_pheno),  
                    rep(NA, num_pheno), 
                    rep(NA, num_pheno), 
                    rep(NA, num_pheno))

colnames(thresholds) <- c(colnames(phenotypes),
                          "indiv_pval", 
                          "indiv_epsilon", 
                          "group_pval", 
                          "group_epsilon")

thresholds$indiv_pval <- thresholds$group_pval  <- rep(0.15, num_pheno)
thresholds$indiv_epsilon <- thresholds$group_epsilon  <- rep(0.15, num_pheno)


write_csv(thresholds, "../../data/7_thresholds/phenotype_thresholds.csv", col_names = TRUE)
