load("../../data/1_combine_genotype_matrices/cleaned_log_toxin_pan_and_rereferenced_snp_indel.RData") # genotype
load("../../data/5_hogwash/log_toxin/hogwash_continuous_grouped_continuous_log_toxin_grouped_post-ar.rda") # hogwash_continuous
key <- read.table( "../../data/2_gene_keys/cleaned_log_toxin_gene_and_intergenic_key.tsv",
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  check.names = FALSE)
key <- as.matrix(key)

thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv",
                       stringsAsFactors = FALSE)
pval_thresh <- thresholds$group_pval[thresholds$shorthand == "log_toxin"]
epsilon_thresh <- thresholds$group_epsilon[thresholds$shorthand == "log_toxin"]

# copied from hogwash::build_gene_genotype_from_snps_post_ar()

unique_genes <- unique(key[, 2])
samples_by_genes <- matrix(0, nrow = nrow(genotype), ncol = length(unique_genes))
colnames(samples_by_genes) <- unique_genes
row.names(samples_by_genes) <- row.names(genotype)

for (j in 1:length(unique_genes)) {
  temp_mat <-
    genotype[, colnames(genotype) %in% key[ , 1][key[, 2] == unique_genes[j]],
         drop = FALSE]
  class(temp_mat) <- "numeric"
  temp_column <- rowSums(temp_mat)
  samples_by_genes[, j] <- temp_column
}

samples_by_genes <- samples_by_genes > 0
class(samples_by_genes) <- "numeric"

sig_index <- hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(pval_thresh)

samples_by_genes <- samples_by_genes[, colnames(samples_by_genes) %in% row.names(hogwash_continuous$hit_pvals), drop = FALSE]

samples_by_genes <- samples_by_genes[, match(row.names(hogwash_continuous$hit_pvals), colnames(samples_by_genes)), drop = FALSE]

if(!identical(colnames(samples_by_genes), row.names(hogwash_continuous$hit_pvals))){stop("ERROR")}
if(!identical(length(sig_index), ncol(samples_by_genes))){stop("ERROR")}

sig_genotype <- samples_by_genes[, sig_index, drop = FALSE]

print(sum(sig_index))
print(dim(samples_by_genes))
print(dim(sig_genotype))
write.table(x = sig_genotype, 
            quote = TRUE,
            file = "../../data/1_combine_genotype_matrices/combinedmat/cleaned_log_toxin_sig_group_pan_and_rereferenced_snp_indel.tsv",
            sep = "\t", 
            col.names = TRUE, 
            row.names = TRUE)

save(sig_genotype, 
     file = "../../data/1_combine_genotype_matrices/combinedmat/cleaned_toxin_sig_group_pan_and_rereferenced_snp_indel.RData")

