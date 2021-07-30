library(tidyverse)
# Read in the data
load("../../data/5_hogwash/log_toxin/hogwash_continuous_grouped_continuous_log_toxin_grouped_post-ar.rda")
group_obj <- hogwash_continuous
hogwash_continuous <- NULL

thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv",
  stringsAsFactors = FALSE
)
group_pval_thresh <- thresholds$group_pval[thresholds$shorthand == "log_toxin"]

load(file = "../../data/1_combine_genotype_matrices/snpmat/snp_annots.RData")
load(file = "../../data/1_combine_genotype_matrices/indelmat/indel_annots.RData")

# Grouped analysis
# What are the flagellar genes in CD630?
# CD630_02260- CD630_02720 - or 293,002 - 333,302
# So let's read in the SNP annotations and get all of the gene IDs with positions between 293,002 - 333,302
flagellar_genes <- snp_annotations %>%
  filter(
    pos < 333303 & pos > 293001,
    intergenic == FALSE
  ) %>%
  select(locus_tag, annotation_1) %>%
  unique()

temp_genes <- indel_annotations %>%
  filter(
    pos < 333303 & pos > 293001,
    intergenic == FALSE
  ) %>%
  select(locus_tag, annotation_1) %>%
  unique()

flagellar_genes <- rbind(flagellar_genes, temp_genes)

uni_geno <- unique(flagellar_genes$locus_tag)
fixed_df <- as.data.frame(matrix(NA, nrow = 0, ncol = 2))

for (i in 1:length(uni_geno)) {
  temp_df <- flagellar_genes %>% filter(locus_tag == uni_geno[i])

  new_annotations <- temp_df %>%
    pull(annotation_1) %>%
    unique() %>%
    str_flatten(., collapse = "; ")

  temp_df$annotation_1 <- new_annotations

  temp_df <- temp_df %>% unique()
  fixed_df <- rbind(fixed_df, temp_df)
}
flagellar_genes <- fixed_df


# Are flagellar genes enriched in the high significance set of hogwash results?
# To calculate this, let's do a permutation test:
#
# * Calculate the number of significant loci, e.g. those above the P-value threshold, called X.
# * Calculate the number of flagellar genes that are significant. This is the "observed" number.
# * Create a vector that classifies each genotype assigned a P-value by hogwash as either a flageller gene or not.
# * For 10,000 permutations: sample X entries from the vector and count the number of flagellar genes to generate the null distribution.
# * Calculate P-value: (1 + (number of null values >= observed)) / (1 + num permutations)

# To account also for intergenic regions we need to add some nuance here
grouped_genotype_names <- row.names(group_obj$hit_pvals)

grouped_genotype_genes_and_intergenic <-
  as.data.frame(grouped_genotype_names) %>%
  mutate(
    locus_tag1 = gsub("-.*", "", grouped_genotype_names),
    locus_tag2 = gsub(".*-", "", grouped_genotype_names)
  )

grouped_flagellar_genes_and_intergenic <-
  grouped_genotype_genes_and_intergenic %>%
  filter(locus_tag1 %in% flagellar_genes$locus_tag |
    locus_tag2 %in% flagellar_genes$locus_tag)

flagellar_gene_log <-
  row.names(group_obj$hit_pvals) %in%
  grouped_flagellar_genes_and_intergenic$grouped_genotype_names

num_sig <- sum(group_obj$hit_pvals$fdr_corrected_pvals > -log(group_pval_thresh))

num_sig_flagellar <-
  sum(flagellar_gene_log &
    group_obj$hit_pvals$fdr_corrected_pvals > -log(group_pval_thresh))

num_perm <- 10000

null_num_sig_flagellar <- rep(NA, num_perm)

# Run permutation test
set.seed(1)
for (i in 1:num_perm) {
  null_num_sig_flagellar[i] <-
    sum(sample(x = flagellar_gene_log, size = num_sig, replace = FALSE))
}
flagellar_pval <-
  (1 + sum(null_num_sig_flagellar >= num_sig_flagellar)) / (1 + num_perm)

flagella_summary_mat <- matrix(NA, nrow = 8, ncol = 2)
flagella_summary_mat[, 1] <- c(
  "Observed_sig_flag_gene_or_intergenic",
  "Pval",
  "Null median",
  "Null mean",
  "Null min",
  "Null max",
  "Num flagellar genes in grouped results",
  "Num flagellar intergenic regions in grouped results"
)
flagella_summary_mat[1, 2] <- num_sig_flagellar
flagella_summary_mat[2, 2] <- round(flagellar_pval, 5)
flagella_summary_mat[3, 2] <- median(null_num_sig_flagellar)
flagella_summary_mat[4, 2] <- mean(null_num_sig_flagellar)
flagella_summary_mat[5, 2] <- min(null_num_sig_flagellar)
flagella_summary_mat[6, 2] <- max(null_num_sig_flagellar)
flagella_summary_mat[7, 2] <-
  grouped_flagellar_genes_and_intergenic %>%
  select(grouped_genotype_names) %>%
  filter(!grepl("-", grouped_genotype_names)) %>%
  nrow()
flagella_summary_mat[8, 2] <-
  grouped_flagellar_genes_and_intergenic %>%
  select(grouped_genotype_names) %>%
  filter(grepl("-", grouped_genotype_names)) %>%
  nrow()

write.csv(flagella_summary_mat,
  "../../data/9_flagella/flagella_summary.csv",
  quote = FALSE,
  row.names = FALSE
)

png("../../figures/flagella_permutation_test/flagella_permutation_test.png")
hist(null_num_sig_flagellar,
     xlim = c(0, num_sig_flagellar + 5),
     col = "grey",
     breaks = 10,
     ylab = "Count",
     xlab = "Significant flagellar genes",
     main = paste0("Null distribution of significant flagellar genes\n in grouped toxin results; P-value=",
                        round(flagellar_pval, 6)))
abline(v = num_sig_flagellar, col = "red")
dev.off()

# names of significant flagellar genes
sig_gene_names <-
  row.names(group_obj$hit_pvals)[
    group_obj$hit_pvals$fdr_corrected_pvals > -log(group_pval_thresh)
  ]

sig_flagellar_genes <-
  flagellar_genes %>% filter(locus_tag %in% sig_gene_names)

write.csv(
  flagellar_genes,
  "../../data/9_flagella/flagella_gene_names_and_annotation.csv",
  quote = FALSE,
  row.names = FALSE
)