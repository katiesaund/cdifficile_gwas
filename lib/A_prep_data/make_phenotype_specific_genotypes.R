# this script will generate a genotype specific to just the samples in each phenotype (phenotypes have non-identical sets of genotypes)
source("data_prep_lib.R")

print("start reading phenotypes")
# Read in phenotypes
fqR <-  read.table(file = "../../data/4_phenotypes/cleaned_fqR.tsv",
                   sep = "\t", 
                   header = TRUE,
                   row.names = 1, 
                   stringsAsFactors = FALSE)

cfe <-  read.table(file = "../../data/4_phenotypes/cleaned_log_cfe.tsv",
                   sep = "\t", 
                   header = TRUE,
                   row.names = 1, 
                   stringsAsFactors = FALSE)

germ_tc_and_gly <-  read.table(file = "../../data/4_phenotypes/cleaned_log_germ_tc_and_gly.tsv",
                   sep = "\t", 
                   header = TRUE,
                   row.names = 1, 
                   stringsAsFactors = FALSE)

germ_tc <-  read.table(file = "../../data/4_phenotypes/cleaned_log_germ_tc.tsv",
                   sep = "\t", 
                   header = TRUE,
                   row.names = 1, 
                   stringsAsFactors = FALSE)

growth <-  read.table(file = "../../data/4_phenotypes/cleaned_log_growth.tsv",
                   sep = "\t", 
                   header = TRUE,
                   row.names = 1, 
                   stringsAsFactors = FALSE)

sporulation <-  read.table(file = "../../data/4_phenotypes/cleaned_log_sporulation.tsv",
                   sep = "\t", 
                   header = TRUE,
                   row.names = 1, 
                   stringsAsFactors = FALSE)

toxin <-  read.table(file = "../../data/4_phenotypes/cleaned_log_toxin.tsv",
                   sep = "\t", 
                   header = TRUE,
                   row.names = 1, 
                   stringsAsFactors = FALSE)

severity <-  read.table(file = "../../data/4_phenotypes/cleaned_severity.tsv",
                   sep = "\t", 
                   header = TRUE,
                   row.names = 1, 
                   stringsAsFactors = FALSE)

print("end pheno")
all_geno_mat <- data.table::fread(file = "../../data/1_combine_genotype_matrices/combinedmat/pan_and_rereferenced_snp_indel.tsv", 
                  sep = "\t", data.table = FALSE, 
                  stringsAsFactors = FALSE, 
                  check.names = FALSE)
all_geno_mat <- as.matrix(all_geno_mat)
row.names(all_geno_mat) <- all_geno_mat[, 1]
all_geno_mat <- all_geno_mat[, 2:ncol(all_geno_mat)]
storage.mode(all_geno_mat) <- "numeric"

print("end geno mat")
all_geno_key <-  read.table(file = "../../data/2_gene_keys/gene_and_intergenic_key.tsv",
                            sep = "\t", 
                            header = TRUE,
                            row.names = NULL, 
                            stringsAsFactors = FALSE)


print("end key")
# Save phenotype-specific genotypes
save_pheno_specific_geno("cleaned_fqR", fqR, all_geno_mat, all_geno_key)
save_pheno_specific_geno("cleaned_log_cfe",cfe, all_geno_mat, all_geno_key)
save_pheno_specific_geno("cleaned_log_germ_tc_and_gly",germ_tc_and_gly, all_geno_mat, all_geno_key)
save_pheno_specific_geno("cleaned_log_germ_tc",germ_tc, all_geno_mat, all_geno_key)
save_pheno_specific_geno("cleaned_log_growth",growth, all_geno_mat, all_geno_key)
save_pheno_specific_geno("cleaned_log_sporulation",sporulation, all_geno_mat, all_geno_key)
save_pheno_specific_geno("cleaned_log_toxin",toxin, all_geno_mat, all_geno_key)
save_pheno_specific_geno("cleaned_severity",severity, all_geno_mat, all_geno_key)

print("end script")
