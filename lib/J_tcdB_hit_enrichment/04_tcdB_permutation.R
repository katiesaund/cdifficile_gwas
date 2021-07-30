options(scipen = 999) # turn off scientific notation
library(tidyverse)
load("../../data/5_hogwash/log_toxin/hogwash_continuous_continuous_log_toxin.rda")

thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv",
  stringsAsFactors = FALSE
)
pval_thresh <- thresholds$indiv_pval[thresholds$shorthand == "log_toxin"]

load("../../data/1_combine_genotype_matrices/cleaned_log_toxin_paloc_only_pan_and_rereferenced_snp_indel.RData")

### Are PaLoc genes enriched in the high significance set of hogwash results?
# * Calculate the number of significant loci, e.g. those above the P-value threshold, called X.
# * Calculate the number of PaLoc variants that are significant. This is the "observed" number.
# * Create a vector that classifies each genotype assigned a P-value by hogwash as either a PaLoc variant or non-PaLoc variant.
# * For 10,000 permutations: sample X entries from the vector and count the number of PaLoc genes to generate the null distribution.
# * Calculate P-value: (1 + (number of null values >= observed)) / (1 + num permutations)

# Set up for permutation test
locus_tags_for_tcd_genes <- c(
  "CD630_06590",
  "CD630_06600",
  "CD630_06610",
  "CD630_06630",
  "CD630_06640"
)

paloc_key_df <- as.data.frame(
  cbind(
    c("CD630_06590", "CD630_06600", "CD630_06610", "CD630_06630", "CD630_06640"),
    c("tcdR", "tcdB", "tcdE", "tcdA", "tcdC")
  )
)
colnames(paloc_key_df) <- c("locus_tag", "gene_name")
# "tcdR", "CD630_06590"
# "tcdB" "CD630_06600"
# "tcdE" "CD630_06610"
# "tcdA" "CD630_06630"
# "tcdC" "CD630_06640"

paloc_gene_log <- as.logical(
  grepl(locus_tags_for_tcd_genes[1], row.names(hogwash_continuous$hit_pvals)) +
    grepl(locus_tags_for_tcd_genes[2], row.names(hogwash_continuous$hit_pvals)) +
    grepl(locus_tags_for_tcd_genes[3], row.names(hogwash_continuous$hit_pvals)) +
    grepl(locus_tags_for_tcd_genes[4], row.names(hogwash_continuous$hit_pvals)) +
    grepl(locus_tags_for_tcd_genes[5], row.names(hogwash_continuous$hit_pvals))
)

tcdB_var_log <-
  as.logical(grepl("CD630_06600", row.names(hogwash_continuous$hit_pvals)))
# Can include indels, snps, and intergenic stuff

# sum(tcdB_var_log) # 302
num_sig <-
  sum(hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(pval_thresh))
num_sig_tcdB <-
  sum(tcdB_var_log & hogwash_continuous$hit_pvals$fdr_corrected_pvals >
    -log(pval_thresh))
num_perm <- 10000
null_num_sig_tcdB <- rep(NA, num_perm)

# Run permutation test

set.seed(1)
for (i in 1:num_perm) {
  null_num_sig_tcdB[i] <-
    sum(sample(x = tcdB_var_log, size = num_sig, replace = FALSE))
}

tcdB_pval <- (1 + sum(null_num_sig_tcdB >= num_sig_tcdB)) / (1 + num_perm)


num_sig_intergenic_tcdB <-
  sum(
    tcdB_var_log &
      hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(pval_thresh) &
      as.logical(grepl("intergenic", row.names(hogwash_continuous$hit_pvals)))
  )

num_sig_genic_tcdB <- sum(
  tcdB_var_log &
    hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(pval_thresh) &
    as.logical(!grepl("intergenic", row.names(hogwash_continuous$hit_pvals)))
)

tcdB_summary_mat <- matrix(NA, nrow = 10, ncol = 2)
tcdB_summary_mat[, 1] <- c(
  "Observed_sig_tcdB_variant",
  "Pval",
  "Null median",
  "Null mean",
  "Null min",
  "Null max",
  "Num tcdB SNPs genes in indiv results",
  "Num tcdB Indels in indiv results",
  "Num significant tcdB intergenic variants",
  "Num significant variants within tcdB"
)
tcdB_summary_mat[1, 2] <- num_sig_tcdB
tcdB_summary_mat[2, 2] <- round(tcdB_pval, 5)
tcdB_summary_mat[3, 2] <- median(null_num_sig_tcdB)
tcdB_summary_mat[4, 2] <- mean(null_num_sig_tcdB)
tcdB_summary_mat[5, 2] <- min(null_num_sig_tcdB)
tcdB_summary_mat[6, 2] <- max(null_num_sig_tcdB)
tcdB_summary_mat[7, 2] <-
  sum(grepl("SNP", row.names(hogwash_continuous$hit_pvals)[tcdB_var_log]))
tcdB_summary_mat[8, 2] <-
  sum(grepl("Indel", row.names(hogwash_continuous$hit_pvals)[tcdB_var_log]))
tcdB_summary_mat[9, 2] <- num_sig_intergenic_tcdB
tcdB_summary_mat[10, 2] <- num_sig_genic_tcdB

write.csv(tcdB_summary_mat,
  "../../data/13_summaries/tcdB_summary.csv",
  quote = FALSE,
  row.names = FALSE
)

png("../../figures/tcdb_permutation_test/tcdb_permutation_test.png")
hist(null_num_sig_tcdB,
  xlim = c(0, num_sig_tcdB + 5),
  col = "grey",
  breaks = 10,
  ylab = "Count",
  xlab = "Significant tcdB loci",
  main = paste0(
    "Null distribution of significant loci\n in individual toxin results; P-value=",
    round(tcdB_pval, 6)
  )
)
abline(v = num_sig_tcdB, col = "red")
dev.off()
