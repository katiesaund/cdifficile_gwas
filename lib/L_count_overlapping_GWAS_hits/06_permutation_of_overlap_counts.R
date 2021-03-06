# See if any degree of GWAS results overlap between phenotypes is significant.

source("../I_summarize_all_gwas_results/gwas_analysis_lib.R") 
phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")

# Read in genotype results for each phenotype
# find intersection of all hits (not just significant)
# figure out if true values more or less than expected by chance.

thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv",
                       stringsAsFactors = FALSE)
p_thresh <- thresholds$indiv_pval[thresholds$shorthand == "log_toxin"]
e_thresh <- thresholds$indiv_epsilon[thresholds$shorthand == "log_toxin"]

get_hits_names_df <- function(shorthand, pv_thresh, ep_thresh, phenotypes){
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
             pheno = temp_shorthand)
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
    indiv_disc_df <- indiv_ranked_sync_df
    return(indiv_disc_df)
  }
}

cfe_df <- get_hits_names_df("log_cfe", p_thresh, e_thresh, phenotypes)
germ_tc_and_gly_df <- get_hits_names_df("log_germ_tc_and_gly", p_thresh, e_thresh, phenotypes)
germ_tc_df <- get_hits_names_df("log_germ_tc", p_thresh, e_thresh, phenotypes)
growth_df <- get_hits_names_df("log_growth", p_thresh, e_thresh, phenotypes)
sporulation_df <- get_hits_names_df("log_sporulation", p_thresh, e_thresh, phenotypes)
toxin_df <- get_hits_names_df("log_toxin", p_thresh, e_thresh, phenotypes)
fqR_df <- get_hits_names_df("fqR", p_thresh, e_thresh, phenotypes)
severity_df <- get_hits_names_df("severity", p_thresh, e_thresh, phenotypes)



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
save(intersect_df_list, file = "../../data/13_summaries/overlapping_permutation_tested_genotype_names.rda")
write.csv(intersect_mat, file = "../../data/13_summaries/overlapping_permutation_tested_genotype_names.csv")


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
                   filename = "../../figures/overlapping/hogwash_colocalization_heatmap_raw_values_all_overlapping_tested_genotypes.jpg",
                   width = 8,
                   height = 3,
                   color = colorRampPalette(c("white", "red"))(100), 
                   display_numbers = FALSE,
                   na_col = "grey",
                   cluster_rows = FALSE,
                   angle_col = 45,
                   cluster_cols = FALSE)


# Add row and column annotations with numbers of significant hits for each phenotype
num_sig_hits_df <- read.table(file = "../../data/13_summaries/number_significant_and_convergent_loci_indiv_gwas_for_overlap_plots.tsv",
                              sep = "\t")

num_loci_df <- as.data.frame(matrix(c(rep(0, 8), phenotypes$longhand), ncol = 2), stringsAsFactors = FALSE)
colnames(num_loci_df) <- c("Tested_Loci", "Phenotype")
num_loci_df$Tested_Loci[num_loci_df$Phenotype == "Severity"] <- nrow(severity_df)
num_loci_df$Tested_Loci[num_loci_df$Phenotype == "ln(spore viability)"] <- nrow(cfe_df)
num_loci_df$Tested_Loci[num_loci_df$Phenotype == "ln(germination in Tc and Gly)"] <- nrow(germ_tc_and_gly_df)
num_loci_df$Tested_Loci[num_loci_df$Phenotype == "ln(germination in Tc)"] <- nrow(germ_tc_df)
num_loci_df$Tested_Loci[num_loci_df$Phenotype == "ln(growth rate)"] <- nrow(growth_df)
num_loci_df$Tested_Loci[num_loci_df$Phenotype == "ln(# spores)"] <- nrow(sporulation_df)
num_loci_df$Tested_Loci[num_loci_df$Phenotype == "ln(toxin activity)"] <- nrow(toxin_df)
row.names(num_loci_df) <- num_loci_df$Phenotype
num_loci_df <- num_loci_df[, 1, drop = FALSE]
num_loci_df$Tested_Loci <- as.numeric(num_loci_df$Tested_Loci)
num_loci_df <- cbind(num_loci_df, num_sig_hits_df)
names(num_loci_df)[2] <- "Sig_Loci"
# drop fqR, spore viability, and germination in Tc and Gly because not interested in fqR and don't need to do permutation
# for phenotypes with no significant loci. 
num_loci_df <- num_loci_df[4:8, ]

intersect_mat <- intersect_mat[row.names(intersect_mat) %in% row.names(num_loci_df),
                               colnames(intersect_mat) %in% row.names( num_loci_df)]
intersect_mat
sig_overlapping_mat <- read.csv(file = "../../data/13_summaries/number_overlapping_gwas_hits.csv", row.names = 1)

# What is this permutation test? 
# Let's say for toxin + severity:
# We're asking: what are the chances that all 8 of the sig. toxin loci are also significant in
# the severity gwas? 
num_perm <- 1000

num_sev_all <- 9115
num_sev_sig <- 1052
num_tox_all <- 18362
num_tox_sig <- 8
num_overlap_sev_tox <- 9114
num_sig_overlap_sev_tox <- 7
# 9114/18362 toxin tested genes were also tested in severity 
# 9114/9115 sev tested genes were also tested in toxin. 
# 7/8 toxin hits are also a hit in severity. 
# Is that more or less than expected by chance given the number of shared loci tested?

# Let's simplify and say that 9114/9115 == 1

null_sev_tox_sig_overlap_vec <- rep(NA, num_perm)
set.seed(1)
for (i in 1:num_perm) {
  null_sig_sev_index <- sample(1:num_sev_all, size = num_sev_sig, replace = FALSE)
  null_sig_tox_index <- sample((num_sev_all - num_overlap_sev_tox + 1):num_tox_all, size = num_tox_sig, replace = FALSE)
  num_null_sig_overlap <- length(intersect(null_sig_sev_index, null_sig_tox_index))
  null_sev_tox_sig_overlap_vec[i] <- num_null_sig_overlap
}

hist(null_sig_tox_index, breaks = 100, col = rgb(1, 0, 0, 0.5), xlim = c(0, 20000), ylim = c(0, 200))
hist(null_sig_sev_index, breaks = 10, col = rgb(0, 1, 0, 0.5), add = TRUE)

hist(null_sev_tox_sig_overlap_vec)
abline(v = num_sig_overlap_sev_tox, col = "red")

tox_sev_pval <-
  (1 + sum(null_sev_tox_sig_overlap_vec >= num_sig_overlap_sev_tox)) / (1 + num_perm)

# Alright, now how do I implement this for all pairwise interactions? 

perm_pval_mat <- matrix(data = NA, nrow = nrow(intersect_mat), ncol = nrow(intersect_mat))
colnames(perm_pval_mat) <- colnames(intersect_mat)
row.names(perm_pval_mat) <- row.names(intersect_mat)

for (i in 1:nrow(sig_overlapping_mat)) {
  for (j in 1:nrow(phenotypes)) {
    if (row.names(sig_overlapping_mat)[i] == phenotypes$shorthand[j]) {
      row.names(sig_overlapping_mat)[i] <- phenotypes$longhand[j]
    }
  }
}

for (i in 1:ncol(sig_overlapping_mat)) {
  for (j in 1:nrow(phenotypes)) {
    if (colnames(sig_overlapping_mat)[i] == phenotypes$shorthand[j]) {
      colnames(sig_overlapping_mat)[i] <- phenotypes$longhand[j]
    }
  }
}

set.seed(1)
for (r in 1:nrow(perm_pval_mat)) {
  for (c in 1:ncol(perm_pval_mat)) {
    print("start")
    print(row.names(perm_pval_mat)[r])
    print(colnames(perm_pval_mat)[c])
    
    if (row.names(perm_pval_mat)[r] != colnames(perm_pval_mat)[c]) {
      num_r_all <- num_loci_df$Tested_Loci[row.names(num_loci_df) == row.names(perm_pval_mat)[r]]
      num_r_sig <- num_loci_df$Sig_Loci[row.names(num_loci_df) == row.names(perm_pval_mat)[r]]
      
      num_c_all <- num_loci_df$Tested_Loci[row.names(num_loci_df) == colnames(perm_pval_mat)[c]]
      num_c_sig <- num_loci_df$Sig_Loci[row.names(num_loci_df) == colnames(perm_pval_mat)[c]]
      
      num_overlap_r_c <- intersect_mat[r, c]
      
      num_sig_overlap_r_c <- sig_overlapping_mat[row.names(sig_overlapping_mat) == row.names(perm_pval_mat)[r], 
                                                 colnames(sig_overlapping_mat) == colnames(perm_pval_mat)[c]]
      
      null_r_c_sig_overlap_vec <- rep(NA, num_perm)
      for (i in 1:num_perm) {
        
        if (num_r_all < num_c_all) {
          null_sig_r_index <- sample(1:num_r_all, size = num_r_sig, replace = FALSE)
          null_sig_c_index <- sample((num_r_all - num_overlap_r_c + 1):num_c_all, size = num_c_sig, replace = FALSE)
        } else {
          null_sig_c_index <- sample(1:num_c_all, size = num_c_sig, replace = FALSE)
          null_sig_r_index <- sample((num_c_all - num_overlap_r_c + 1):num_r_all, size = num_r_sig, replace = FALSE)
        }
        
        num_null_sig_overlap <- length(intersect(null_sig_r_index, null_sig_c_index))
        null_r_c_sig_overlap_vec[i] <- num_null_sig_overlap
      }
      
      hist(null_r_c_sig_overlap_vec, xlim = c(0, max(c(null_r_c_sig_overlap_vec, num_sig_overlap_r_c))))
      abline(v = num_sig_overlap_r_c, col = "red")
      
      r_c_pval <-
        (1 + sum(null_r_c_sig_overlap_vec >= num_sig_overlap_r_c)) / (1 + num_perm)
      
      perm_pval_mat[r, c] <- r_c_pval
      
    }
  }
}

perm_pval_mat

# Next, get the non-duplicated matrix
# do bonferroni multiple test correction. 

perm_pval_mat[5, ] <- c(NA, NA, NA, NA, NA)
perm_pval_mat[4, 3:5] <- c(NA, NA, NA)
perm_pval_mat[3, 4:5] <- c(NA, NA)
perm_pval_mat[2, 5] <- NA


round(p.adjust(perm_pval_mat[!is.na(perm_pval_mat)], method = "bonf"), 7)
# Because there are 10 p-values, the bonf correct pvals are just these all x10

perm_pval_mat <- 10 * perm_pval_mat
perm_pval_mat[perm_pval_mat == 10] <- 1
perm_pval_mat <- round(perm_pval_mat, 5)

write.csv(perm_pval_mat, file = "../../data/13_summaries/overlapping_permutation_bonf_corrected_p_values.csv")


