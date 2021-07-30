load("../../data/1_combine_genotype_matrices/cleaned_log_toxin_pan_and_rereferenced_snp_indel.RData") # genotype
load("../../data/5_hogwash/log_toxin/hogwash_continuous_continuous_log_toxin.rda") #hogwash_continuous

thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv",
                       stringsAsFactors = FALSE)
pval_thresh <- thresholds$indiv_pval[thresholds$shorthand == "log_toxin"]
epsilon_thresh <- thresholds$indiv_epsilon[thresholds$shorthand == "log_toxin"]

# stats::p.adjust(hit_values, method = "fdr")
fdr_corrected_pvals <- -log(p.adjust(exp(-hogwash_continuous$raw_pvals$neg_log_unadjusted_pvals), method = "fdr"))
foo <- fdr_corrected_pvals > -log(pval_thresh)

sig_index <- hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(pval_thresh)
if(!identical(foo, sig_index)){stop("ERROR")}
genotype <- genotype[, colnames(genotype) %in% row.names(hogwash_continuous$hit_pvals), drop = FALSE]

genotype <- genotype[, match(row.names(hogwash_continuous$hit_pvals), colnames(genotype)), drop = FALSE]

if(!identical(colnames(genotype), row.names(hogwash_continuous$hit_pvals))){stop("ERROR")}
if(!identical(ncol(genotype), length(sig_index))){stop("ERROR")}

sig_genotype <- genotype[, sig_index, drop = FALSE]

print(sum(sig_index))
print(dim(genotype))
print(dim(sig_genotype))
write.table(x = sig_genotype, 
            quote = TRUE,
            file = "../../data/1_combine_genotype_matrices/combinedmat/cleaned_toxin_sig_indiv_pan_and_rereferenced_snp_indel.tsv",
            sep = "\t", 
            col.names = TRUE, 
            row.names = TRUE)

save(sig_genotype, 
     file = "../../data/1_combine_genotype_matrices/combinedmat/cleaned_toxin_sig_indiv_pan_and_rereferenced_snp_indel.RData")
