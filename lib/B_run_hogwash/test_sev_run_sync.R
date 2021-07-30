library(hogwash)
library(data.table)
library(ape)
print("session Info")
print(sessionInfo())
args <- commandArgs(trailingOnly = TRUE) # Grab arguments from the .sbat file
phenotype <- read.table("../../data//4_phenotypes/cleaned_severity.tsv",
                        sep = "\t",
                        row.names = 1,
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        check.names = FALSE)
phenotype <- as.matrix(phenotype)

# Because this genotype is VERY big, we're using fread instead of read.table
genotype <- data.table::fread(file = "../../data//1_combine_genotype_matrices/combinedmat/cleaned_severity_pan_and_rereferenced_snp_indel.tsv",
                              sep = "\t", 
                              data.table = FALSE, 
                              stringsAsFactors = FALSE, 
                              check.names = FALSE)
genotype <- as.matrix(genotype)
row.names(genotype) <- genotype[, 1]
genotype <- genotype[, 2:300]
storage.mode(genotype) <- "numeric"

tr  <- ape::read.tree("../../data//3_trees/severity.tree")
f_name <- "test_sync"
out_dir <- "../../data//5_hogwash/severity/"
perm_num <- 1000
fdr_value <- as.numeric(0.3)
bootstrap_threshold <- as.numeric(0.90)
test_type <- as.character("synchronous")

group_method <- as.character("post-ar")

strain_key <- read.table("../../data//6_ribotype/severity_ribotype_no_colors.tsv",
                         sep = "\t",
                         row.names = 1,
                         header = TRUE,
                         stringsAsFactors = FALSE,
                         check.names = FALSE)
strain_key <- as.matrix(strain_key)

key <- NULL

print("finish reading inputs, start hogwash")
hogwash(pheno = phenotype,
        geno = genotype,
        tree = tr,
        file_name = f_name,
        dir = out_dir,
        perm = perm_num,
        fdr = fdr_value,
        bootstrap = bootstrap_threshold,
        test = test_type, 
        grouping_method = group_method, 
        group_genotype_key = key, 
        strain_key = strain_key)
print("finish hogwash")

