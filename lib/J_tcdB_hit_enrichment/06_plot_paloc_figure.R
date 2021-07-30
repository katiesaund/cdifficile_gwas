# Code copied from notebook/Annotate_PaLoc_figure.Rmd
library(tidyverse)
source("color_palettes.R")

load("../../data/1_combine_genotype_matrices/cleaned_log_toxin_paloc_only_pan_and_rereferenced_snp_indel.RData")
snp_geno <- paloc_genotype[, grepl("SNP", colnames(paloc_genotype)), drop = FALSE]
indel_geno <- paloc_genotype[, grepl("Indel", colnames(paloc_genotype)), drop = FALSE]

tested_snps <- gsub(".*at ", "", colnames(snp_geno)) %>% gsub(" >.*", "", .) %>% as.numeric()
tested_indels <- gsub(".*at ", "", colnames(indel_geno)) %>% gsub(" >.*", "", .) %>% as.numeric()

thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv",
                       stringsAsFactors = FALSE)
indiv_pval_thresh <- thresholds$indiv_pval[thresholds$shorthand == "log_toxin"]
indiv_epsilon_thresh <- thresholds$indiv_epsilon[thresholds$shorthand == "log_toxin"]

load("../../data/5_hogwash/log_toxin/hogwash_continuous_continuous_log_toxin.rda")

significant_loci <- row.names(hogwash_continuous$hit_pvals)[hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(indiv_pval_thresh)]
significant_loci <- significant_loci[significant_loci %in% colnames(paloc_genotype)]
significant_snp <- significant_loci[grepl("SNP", significant_loci)]
significant_indel <- significant_loci[grepl("Indel", significant_loci)]

significant_snp <- gsub(".*at ", "", significant_snp) %>% gsub(" >.*", "", .) %>% as.numeric()
significant_indel <- gsub(".*at ", "", significant_indel) %>% gsub(" >.*", "", .) %>% as.numeric()

# Start / stop genomic positions
# https://www.ncbi.nlm.nih.gov/nuccore/AM180355.1?&withparts=on&expand-gaps=on

cdu1_start <- 785233
cdu1_stop <- 785613
tcdR_start <- 786521
tcdR_stop <- 787075
tcdB_start <- 787393
tcdB_stop <- 794493
tcdE_start <- 794615
tcdE_stop <- 795115
uviB_pseudo_start <- 795224
uviB_pseudo_stop <- 795355
tcdA_start <- 795843
tcdA_stop <- 803975
tcdC_comp_start <- 804310
tcdC_comp_stop <- 805008 
hyp_protein_start <- 805434
hyp_protein_stop <- 805679

# To get the promoter start and stop locations I blasted the promoter sequences from Hundsberger et al 1997 vs a PaLoc_plus.fna (cdu1 -> hypothetical_protein). Grabbed the results from command line blast: data/10_paloc/promoters_vs_paloc.csv. The first position in the blast is cdu1_start (so 1 == 785233).  

promoter_positions <- read.csv("../../data/10_paloc/promoters_vs_paloc.csv", 
                               sep = ",",
                               stringsAsFactors = FALSE, 
                               header = FALSE)
tcdR_promoter1_start <- 904 
tcdR_promoter1_stop <- 946  
tcdR_promoter2_start <- 1068 
tcdR_promoter2_stop <- 1108 

tcdA_promoter_start <- 10337 
tcdA_promoter_stop <- 10376

tcdB_promoter2_start <- 1916 
tcdB_promoter2_stop <- 1963 

tcdC_promoter_start <- 19935 
tcdC_promoter_stop <- 19897 

tcdE_promoter_start <- 9163 
tcdE_promoter_stop <- 9203 

# I had to do an additional command line blast to accomodate less stringent alignment to ge these positions
tcdB_promoter1_start <- 1961 
tcdB_promoter1_stop <- 1996 

# For the tcdR sigD promoter location I looked at two papers: Monot2015 (Fig S6) and El Meouche2013
tcdR_sigD_promoter_start <- tcdR_start - 111 - cdu1_start
tcdR_sigD_promoter_stop <- tcdR_start - 76 - cdu1_start

# Generate matrix of variants
# goal is to get a matrix with columns named start_position_number:end_position_number
# and 4 rows called significant_snps, sig_indel, tested_snps, tested_indels; 0 for absent, 1 for present
# tcdB_mapping_mat

mapping_mat <- matrix(0, nrow = 4, ncol = length(cdu1_start:hyp_protein_stop))
colnames(mapping_mat) <- cdu1_start:hyp_protein_stop
row.names(mapping_mat) <- c("Significant SNPs", "Tested SNPs", "Significant Indels", "Tested Indels")

mapping_mat[1, significant_snp - cdu1_start + 1] <- 2
mapping_mat[2, tested_snps - cdu1_start + 1] <- 1
mapping_mat[3, significant_indel - cdu1_start + 1] <- 2
mapping_mat[4, tested_indels - cdu1_start + 1] <- 1

# Generate annotation dataframe
# row.names(df) == colnames(mat)
# Columns: Gene / intergenic & color genes by name in the color key


my_colour = list(
  Domain = toxin_domain_palette, 
  Gene = just_gene_palette, 
  Promoter = promoter_palette
)

annot_df <- as.data.frame(matrix(NA, nrow = ncol(mapping_mat), ncol = 3))
row.names(annot_df) <- colnames(mapping_mat)
colnames(annot_df) <- c("Promoter", "Gene", "Domain")

annot_df$Promoter <- rep("Not promoter", nrow(annot_df))
annot_df$Promoter[tcdR_promoter1_start:tcdR_promoter1_stop] <- "Promoter"
annot_df$Promoter[tcdR_promoter2_start:tcdR_promoter2_stop] <- "Promoter"
annot_df$Promoter[tcdB_promoter2_start:tcdB_promoter2_stop] <- "Promoter"
annot_df$Promoter[tcdB_promoter1_start:tcdB_promoter1_stop] <- "Promoter"
annot_df$Promoter[tcdE_promoter_start:tcdE_promoter_stop] <- "Promoter"
annot_df$Promoter[tcdA_promoter_start:tcdA_promoter_stop] <- "Promoter"
annot_df$Promoter[tcdC_promoter_start:tcdC_promoter_stop] <- "Promoter"

# We know that there is some overlap in the tcdB promoters so lets handle that
# overlapping_tcdB_promoter_region <- intersect((tcdB_promoter2_start:tcdB_promoter2_stop), (tcdB_promoter1_start:tcdB_promoter1_stop))
# annot_df$Promoter[overlapping_tcdB_promoter_region] <- "Overlapping promoters"
# Actually, let's ignore the three nucleotide overlap b/c it's too hard to see on the plot

# And lets add the tcdR sigmaD promoter
annot_df$Promoter[tcdR_sigD_promoter_start:tcdR_sigD_promoter_stop] <- "sigD promoter"


annot_df$Domain <- rep("NA_domain", nrow(annot_df))
tcdA_index_1 <- tcdA_start - cdu1_start + 1
annot_df$Domain[tcdA_index_1:(tcdA_index_1 + 542 * 3)] <- "Glucosyltransferase"
annot_df$Domain[(tcdA_index_1 + 542 * 3 + 1):(tcdA_index_1 + 801 * 3)] <- "Autoprotease"
annot_df$Domain[(tcdA_index_1 + 801 * 3 + 1):(tcdA_index_1 + 1832 * 3)] <- "Pore-forming"
annot_df$Domain[(tcdA_index_1 + 1832 * 3 + 1):(tcdA_index_1 + 2710 * 3)] <- "CROPs"

tcdB_index_1 <- tcdB_start - cdu1_start + 1
annot_df$Domain[tcdB_index_1:(tcdB_index_1 + 543 * 3)] <- "Glucosyltransferase"
annot_df$Domain[(tcdB_index_1 + 543 * 3 + 1):(tcdB_index_1 + 799 * 3)] <- "Autoprotease"
annot_df$Domain[(tcdB_index_1 + 799 * 3 + 1):(tcdB_index_1 + 1834 * 3)] <- "Pore-forming"
annot_df$Domain[(tcdB_index_1 + 1834 * 3 + 1):(tcdB_index_1 + 2366 * 3)] <- "CROPs"

annot_df$Gene <- rep("intergenic", nrow(annot_df))
annot_df$Gene[(tcdR_start - cdu1_start + 1):(tcdR_stop - cdu1_start + 1)] <- "tcdR"
annot_df$Gene[(tcdB_start - cdu1_start + 1):(tcdB_stop - cdu1_start + 1)] <- "tcdB"
annot_df$Gene[(tcdE_start - cdu1_start + 1):(tcdE_stop - cdu1_start + 1)] <- "tcdE"
annot_df$Gene[(tcdA_start - cdu1_start + 1):(tcdA_stop - cdu1_start + 1)] <- "tcdA"
annot_df$Gene[(tcdC_comp_start - cdu1_start + 1):(tcdC_comp_stop - cdu1_start + 1)] <- "tcdC"
annot_df$Gene[(cdu1_start - cdu1_start + 1):(cdu1_stop - cdu1_start + 1)] <- "cdu1"
annot_df$Gene[(uviB_pseudo_start - cdu1_start + 1):(uviB_pseudo_stop - cdu1_start + 1)] <- "uviB pseudogene"
annot_df$Gene[(hyp_protein_start - cdu1_start + 1):(hyp_protein_stop - cdu1_start + 1)] <- "cdd1"

# I noticed when plotting that the individual locus lines are hard to read so we 
# need to enhance the figure for plotting -- we need to make the bar wider. 
# To do this I think the best approach is to take the mapping_mat and make each 
# variant also take up the 4 positions on either side. Meaning, a variant at pos
# 10 would now extend from 6 - 14. 

widen_bars <- function(map_mat) {
  min_col <- 1
  max_col <- ncol(map_mat)
  
  for (r in 1:nrow(map_mat)) {
    val <- 1
    if (r %in% c(1, 3)) {
      val <- 2
    } 
    col_index <- unname(which(map_mat[r, ] == val))
    
    new_indices <- NULL
    for (i in 1:length(col_index)) {
      index <- col_index[i]

      # catch cases when we're at the ends
      temp_min <- index - 4
      temp_max <- index + 4
      
      while (temp_min < min_col) {
        temp_min <- temp_min + 1
      }
      
      while (temp_max > max_col) {
        temp_max <- temp_max - 1
      }
      
      new_indices <- c(new_indices, temp_min:temp_max)
    }
    new_indices
    map_mat[r, new_indices] <- val
  }
  return(map_mat)
}

mapping_mat_with_wider_bars <- widen_bars(mapping_mat)

# Plot together
halfway <- ncol(mapping_mat_with_wider_bars) / 2
end <- ncol(mapping_mat_with_wider_bars)
pheatmap::pheatmap(mapping_mat_with_wider_bars[, 1:halfway],  
                   width = 11, 
                   height = 8.5,
                   file = "../../figures/paloc_name_tbd_1.pdf",
                   show_colnames = FALSE,  
                   show_rownames = TRUE, 
                   cluster_rows = FALSE,
                   cellheight = 40, 
                   cluster_cols = FALSE, 
                   color = c("white", "grey", "black"), 
                   annotation_col = annot_df, 
                   legend = FALSE,
                   annotation_legend = FALSE,
                   annotation_names_col = TRUE,
                   annotation_colors = my_colour)

pheatmap::pheatmap(mapping_mat_with_wider_bars[, (halfway + 1):end],  
                   width = 11, 
                   height = 8.5,
                   file = "../../figures/paloc_name_tbd_2.pdf",
                   show_colnames = FALSE,  
                   show_rownames = TRUE, 
                   cluster_rows = FALSE,
                   cellheight = 40, 
                   cluster_cols = FALSE, 
                   color = c("white", "grey"), 
                   annotation_col = annot_df, 
                   legend = FALSE,
                   annotation_legend = FALSE,
                   annotation_names_col = TRUE,
                   annotation_colors = my_colour)

pheatmap::pheatmap(mapping_mat_with_wider_bars,  
                   width = 11, 
                   height = 8.5,
                   file = "../../figures/paloc_name_tbd_3.pdf",
                   show_colnames = FALSE,  
                   show_rownames = TRUE, 
                   cluster_rows = FALSE,
                   cellheight = 40, 
                   cluster_cols = FALSE, 
                   color = c("white", "grey", "black"), 
                   annotation_col = annot_df, 
                   legend = FALSE,
                   annotation_names_col = TRUE,
                   annotation_colors = my_colour)


# repeat exact code from 226 -239
pheatmap_plot <- pheatmap::pheatmap(mapping_mat_with_wider_bars,  
                                    width = 11, 
                                    height = 8.5,
                                    show_colnames = FALSE,  
                                    show_rownames = TRUE, 
                                    cluster_rows = FALSE,
                                    cellheight = 40, 
                                    cluster_cols = FALSE, 
                                    color = c("white", "grey", "black"), 
                                    annotation_col = annot_df, 
                                    legend = FALSE,
                                    annotation_names_col = TRUE,
                                    annotation_colors = my_colour)
save(pheatmap_plot, file = "../../figures/paloc_name_tbd_3.RData")


# Does the tcdR-tcdB SNP (individual results) fall in the tcdB promoter? Yes!
  
# Do the cdu1-tcdR variants (grouped results) also fall in the promoter(s)? Yes! 3 fall in the first promoter and 2 fall in the second promoter.
mapping_df <- as.data.frame(t(mapping_mat))
mapping_df <- cbind(mapping_df, annot_df)

# mapping_df %>% filter(`Significant SNPs` == 2 & Promoter == "Promoter")
# mapping_df %>% filter(`Significant Indels` == 2 & Promoter == "Promoter")


tcdR_promoter_df <- mapping_df %>%
  mutate("genomic position" = row.names(mapping_df)) %>%       
  filter(`Tested SNPs` == 1 & Promoter %in%  c("Promoter", "sigD promoter") & `genomic position` < tcdR_start & `genomic position` > cdu1_stop) 
write_csv(tcdR_promoter_df, path = "../../data/13_summaries/tcdR_promoter_summary_for_individual_test.csv", col_names = TRUE)

all_promoter_df <- mapping_df %>%
  mutate("genomic position" = row.names(mapping_df)) %>%       
  filter(`Tested SNPs` == 1 & Promoter !=  "Not promoter")
write_csv(all_promoter_df, path = "../../data/13_summaries/all_promoter_summary_for_individual_test.csv", col_names = TRUE)

#117 deletion in TcdC
tcdC_pvals <- hogwash_continuous$hit_pvals[grepl("tcdC", row.names(hogwash_continuous$hit_pvals)),, drop = FALSE]
tcdC_pvals <- tcdC_pvals[grepl("117delA", row.names(tcdC_pvals)),, drop = FALSE]
tcdC_pvals <- tcdC_pvals$fdr_corrected_pvals

# Make a summary
pval_ep <- cbind(hogwash_continuous$hit_pvals$fdr_corrected_pvals, hogwash_continuous$convergence$epsilon)
colnames(pval_ep) <- c("pval", "epsilon")
pval_ep <- as.data.frame(pval_ep)
paloc_pval_ep <- cbind(pval_ep, row.names(hogwash_continuous$hit_pvals))
colnames(paloc_pval_ep)[3] <- "locus"
paloc_pval_ep <- paloc_pval_ep %>% filter(locus %in% colnames(paloc_genotype))

paloc_summary <- matrix(NA, nrow = 14, ncol = 2)
paloc_summary[, 1] <- c("Number of SNPs assigned p-value in paloc", 
                        "Number of Indels assigned p-value in paloc", 
                        "Total number of indiv loci assigned p-value in paloc", 
                        "Number of indiv loci above p-value threshold in paloc", 
                        "Number of SNPS above p-value threshold in paloc", 
                        "Number of indels above p-value threshold in paloc", 
                        "Number of indiv loci above the epsilon threshold in paloc ", 
                        "Number of indiv loci above both p-value and epsilon threshold in paloc", 
                        "Number of indiv loci in tcdR-tcdB intergenic region", 
                        "Number of indiv loci in tcdR-tcdB above p-value threshold & in promoter", 
                        "Number of indiv loci in cdu1-tcdR intergenic region", 
                        "Number of indiv loci in tcdR (not sigD) promoters", 
                        "Number of indiv loci in tcdR-tcdB above p-value threshold", 
                        "P-value of adenosine deletion at 117 in tcdC")
paloc_summary[1, 2] <- length(tested_snps) # unique SNPs, there are a coupel SNPs at the same loci
paloc_summary[2, 2] <- length(tested_indels)
paloc_summary[3, 2] <- sum(length(tested_snps), length(tested_indels))
paloc_summary[4, 2] <- sum(length(significant_snp), length(significant_indel))
paloc_summary[5, 2] <- length(significant_snp)
paloc_summary[6, 2] <- length(significant_indel)
paloc_summary[7, 2] <- paloc_pval_ep %>% filter(epsilon > indiv_epsilon_thresh) %>% nrow()
paloc_summary[8, 2] <- paloc_pval_ep %>% filter(pval > -log(indiv_pval_thresh) & epsilon > indiv_epsilon_thresh) %>% nrow()
paloc_summary[9, 2] <- length(intersect(tested_snps, tcdR_stop:tcdB_start))
paloc_summary[10, 2] <- length(intersect(significant_snp - cdu1_start, c(tcdB_promoter2_start:tcdB_promoter2_stop, tcdB_promoter1_start:tcdB_promoter1_stop)))
paloc_summary[11, 2] <- length(intersect(tested_snps, cdu1_stop:tcdR_start))
paloc_summary[12, 2] <- length(intersect(tested_snps - cdu1_start, c(tcdR_promoter2_start:tcdR_promoter2_stop, tcdR_promoter1_start:tcdR_promoter1_stop)))
paloc_summary[13, 2] <-  length(intersect(significant_snp, c(tcdR_stop:tcdB_start)))
paloc_summary[14, 2] <- round(tcdC_pvals, 3)

write.csv(paloc_summary, "../../data/13_summaries/paloc_indiv_summary.csv", row.names = FALSE, quote = FALSE)

                         