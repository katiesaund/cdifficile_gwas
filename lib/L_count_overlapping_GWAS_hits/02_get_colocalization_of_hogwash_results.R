# See if any GWAS results overlap between toxin and the other phenotypes. Of particular interest is severity and sporulation. 

source("../I_summarize_all_gwas_results/gwas_analysis_lib.R") 
phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")

# Read in genotype results for each phenotype
# find intersection of results & make a table
# report result intersection 

thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv",
                       stringsAsFactors = FALSE)
p_thresh <- thresholds$indiv_pval[thresholds$shorthand == "log_toxin"]
e_thresh <- thresholds$indiv_epsilon[thresholds$shorthand == "log_toxin"]

get_sig_hits_df <- function(shorthand, pv_thresh, ep_thresh, phenotypes){
  temp_shorthand <- shorthand
  bin_or_cont <- phenotypes %>% filter(shorthand == temp_shorthand) %>% pull(bin_cont)

  if (bin_or_cont == "continuous") {
    indiv_ranked_df <-
    read_tsv(
      paste0(
        "../../data/13_summaries/",
        temp_shorthand,
        "_indiv_genotypes_ranked_by_epsilon_then_p.tsv"
      ),
      col_names = TRUE
    )
    indiv_ranked_df <- indiv_ranked_df %>% 
      mutate(test = "cont", 
             pheno = temp_shorthand) %>% 
      filter(`-ln(FDR P-value)` > -log(pv_thresh), 
             Epsilon > ep_thresh)
    return(indiv_ranked_df)
  } else {
    indiv_ranked_sync_df <-
      read_tsv(
        paste0(
          "../../data/13_summaries/",
          temp_shorthand,
          "_synchronous_indiv_genotypes_ranked_by_epsilon_then_p.tsv"
        ),
        col_names = TRUE
      )
    indiv_ranked_sync_df <- indiv_ranked_sync_df %>% 
      mutate(test = "sync", 
             pheno = temp_shorthand)
    
    indiv_disc_df <- indiv_ranked_sync_df %>% 
      filter(`-ln(FDR P-value)` > -log(pv_thresh), 
             Epsilon > ep_thresh)
    
    return(indiv_disc_df)
  }
}

cfe_df <- get_sig_hits_df("log_cfe", p_thresh, e_thresh, phenotypes)
germ_tc_and_gly_df <- get_sig_hits_df("log_germ_tc_and_gly", p_thresh, e_thresh, phenotypes)
germ_tc_df <- get_sig_hits_df("log_germ_tc", p_thresh, e_thresh, phenotypes)
growth_df <- get_sig_hits_df("log_growth", p_thresh, e_thresh, phenotypes)
sporulation_df <- get_sig_hits_df("log_sporulation", p_thresh, e_thresh, phenotypes)
toxin_df <- get_sig_hits_df("log_toxin", p_thresh, e_thresh, phenotypes)
fqR_df <- get_sig_hits_df("fqR", p_thresh, e_thresh, phenotypes)
severity_df <- get_sig_hits_df("severity", p_thresh, e_thresh, phenotypes)



assign("log_cfe_var", cfe_df)
assign("log_germ_tc_and_gly_var", germ_tc_and_gly_df)
assign("log_germ_tc_var", germ_tc_df)
assign("log_growth_var", growth_df)
assign("log_sporulation_var", sporulation_df)
assign("log_toxin_var", toxin_df)
assign("fqR_var", fqR_df)
assign("severity_var", severity_df)

#shared data

intersect_mat <- matrix(NA, nrow = 8, ncol = 8)
colnames(intersect_mat) <- row.names(intersect_mat) <- phenotypes$shorthand
names_vec <- rep("", 28)
intersect_df_list <- list()
counter <- 0
for (i in 1:8) {
  for (j in 1:8) {
    # if (j < i) {
    counter <- counter + 1
    shared_geno <- intersect(
      eval(parse(text = paste0(phenotypes$shorthand[i] , "_var")))$Genotype, 
      eval(parse(text = paste0(phenotypes$shorthand[j] , "_var")))$Genotype)
    intersect_mat[i, j] <- length(shared_geno)
    
    shared_df_1 <- eval(parse(text = paste0(phenotypes$shorthand[i] , "_var"))) %>% 
      filter(Genotype %in% shared_geno) 
    
    shared_df_2 <- eval(parse(text = paste0(phenotypes$shorthand[j] , "_var"))) %>% 
      filter(Genotype %in% shared_geno)
    
    shared_df <- rbind(shared_df_1, shared_df_2)
    # names(shared_df) <- paste0(phenotypes$shorthand[i], "-", phenotypes$shorthand[j])
    intersect_df_list[[counter]] <- shared_df
    
    names_vec[counter] <- paste0(phenotypes$shorthand[i], "-", phenotypes$shorthand[j])
   #  }
  }
}

names(intersect_df_list) <- names_vec
save(intersect_df_list, file = "../../data/13_summaries/overlapping_gwas_hits.rda")
write.csv(intersect_mat, file = "../../data/13_summaries/number_overlapping_gwas_hits.csv")


# drop fqR
intersect_mat <- intersect_mat[2:8, 2:8]
# reorder matrix 
intersect_mat <- intersect_mat[, 7:1]
# only show half of the matrix
intersect_mat[1, 7] <- NA
intersect_mat[2, 6] <- NA
intersect_mat[3, 5] <- NA
intersect_mat[4, 4] <- NA
intersect_mat[5, 3] <- NA
intersect_mat[6, 2] <- NA
intersect_mat[7, 1] <- NA


for (i in 1:nrow(intersect_mat)) {
  for (j in 1:nrow(phenotypes)) {
    if (row.names(intersect_mat)[i] == phenotypes$shorthand[j]) {
      row.names(intersect_mat)[i] <- phenotypes$longhand[j]
    }
    if (colnames(intersect_mat)[i] == phenotypes$shorthand[j]) {
      colnames(intersect_mat)[i] <- phenotypes$longhand[j]
    }
  }
}

pheatmap::pheatmap(intersect_mat,
                   width = 8, 
                   height = 3, 
                   filename = "../../figures/overlapping/hogwash_colocalization_heatmap_row_scaled.jpg",
                   color = colorRampPalette(c("white", "red"))(100), 
                   display_numbers = FALSE,
                   na_col = "grey",
                   cluster_rows = FALSE,
                   scale = "row", 
                   angle_col = 45,
                   cluster_cols = FALSE)

pheatmap::pheatmap(intersect_mat,
                   filename = "../../figures/overlapping/hogwash_colocalization_heatmap_raw_values.jpg",
                   width = 8,
                   height = 3,
                   color = colorRampPalette(c("white", "red"))(100), 
                   display_numbers = FALSE,
                   na_col = "grey",
                   cluster_rows = FALSE,
                   angle_col = 45,
                   cluster_cols = FALSE)


# Add row and column annotations with numbers of significant hits for each phenotype
num_sig_hits_df <- as.data.frame(matrix(c(rep(0, 8), phenotypes$longhand), ncol = 2), stringsAsFactors = FALSE)
colnames(num_sig_hits_df) <- c("Loci", "Phenotype")
num_sig_hits_df$Loci[num_sig_hits_df$Phenotype == "Severity"] <- nrow(severity_df)
num_sig_hits_df$Loci[num_sig_hits_df$Phenotype == "ln(spore viability)"] <- nrow(cfe_df)
num_sig_hits_df$Loci[num_sig_hits_df$Phenotype == "ln(germination in Tc and Gly)"] <- nrow(germ_tc_and_gly_df)
num_sig_hits_df$Loci[num_sig_hits_df$Phenotype == "ln(germination in Tc)"] <- nrow(germ_tc_df)
num_sig_hits_df$Loci[num_sig_hits_df$Phenotype == "ln(growth rate)"] <- nrow(growth_df)
num_sig_hits_df$Loci[num_sig_hits_df$Phenotype == "ln(# spores)"] <- nrow(sporulation_df)
num_sig_hits_df$Loci[num_sig_hits_df$Phenotype == "ln(toxin activity)"] <- nrow(toxin_df)
row.names(num_sig_hits_df) <- num_sig_hits_df$Phenotype
num_sig_hits_df <- num_sig_hits_df[, 1, drop = FALSE]
num_sig_hits_df$Loci <- as.numeric(num_sig_hits_df$Loci)

pheatmap::pheatmap(intersect_mat,
                   annotation_col = num_sig_hits_df,
                   annotation_row = num_sig_hits_df,
                   filename = "../../figures/overlapping/DRAFT_hogwash_colocalization_heatmap_raw_values_and_annotations_with_num_loci_included.jpg",
                   width = 8,
                   height = 3, 
                   color = colorRampPalette(c("white", "red"))(100), 
                   display_numbers = FALSE,
                   na_col = "grey",
                   cluster_rows = FALSE,
                   angle_col = 45,
                   cluster_cols = FALSE)

write.table(num_sig_hits_df, file = "../../data/13_summaries/number_significant_and_convergent_loci_indiv_gwas_for_overlap_plots.tsv", sep = "\t")


pheatmap::pheatmap(intersect_mat,
                   annotation_col = num_sig_hits_df,
                   annotation_row = num_sig_hits_df,
                   filename = "../../figures/overlapping/DRAFT_hogwash_colocalization_heatmap_raw_values_and_annotations_with_num_loci_included.jpg",
                   width = 8,
                   height = 3, 
                   color = colorRampPalette(c("white", "red"))(100), 
                   display_numbers = FALSE,
                   na_col = "grey",
                   cluster_rows = FALSE,
                   angle_col = 45,
                   cluster_cols = FALSE)
# remove zero phenotypes

bad_phenotypes <- c("ln(germination in Tc and Gly)", "ln(spore viability)")
foo <- intersect_mat[!(row.names(intersect_mat) %in% bad_phenotypes), 
                     !(colnames(intersect_mat) %in% bad_phenotypes)]

pheatmap::pheatmap(foo,
                   annotation_col = num_sig_hits_df,
                   annotation_row = num_sig_hits_df,
                   filename = "../../figures/overlapping/DRAFT_hogwash_colocalization_heatmap_raw_values_and_annotations_with_num_loci_included_zero_phenotypes_removed.jpg",
                   width = 8,
                   height = 3, 
                   color = colorRampPalette(c("white", "red"))(100), 
                   display_numbers = FALSE,
                   na_col = "grey",
                   cluster_rows = FALSE,
                   angle_col = 45,
                   cluster_cols = FALSE)

# Permutation analysis in 06_permutation_of_overlap_counts.R


