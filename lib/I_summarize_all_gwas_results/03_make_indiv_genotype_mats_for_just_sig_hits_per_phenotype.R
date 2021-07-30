# Goals:
#   - Make a matrix that contains just those genotypes that are of interest: by significance and convergence. 
#   - Record some summary statistics (since I'm loading in this huge geno mat anyways): 
#       - # of variants input into hogwash (and type) 
#       - # of variants in final output (and type)
#       - # of variant that are significant (and type)

source("gwas_analysis_lib.R")
phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")

# individual loci
var_type_summary_df <- NULL
for (i in 1:nrow(phenotypes)) {
  temp_shorthand <- phenotypes$shorthand[i]
  genotype <- data.table::fread(file = paste0("../../data/1_combine_genotype_matrices/combinedmat/cleaned_", 
                                              temp_shorthand, 
                                              "_pan_and_rereferenced_snp_indel.tsv"),
                                sep = "\t", 
                                data.table = FALSE, 
                                stringsAsFactors = FALSE, 
                                check.names = FALSE)
  genotype <- as.matrix(genotype)
  row.names(genotype) <- genotype[, 1]
  genotype <- genotype[, 2:ncol(genotype)]
  
  data <- load_indiv_hogwash_rda_for_pheno(temp_shorthand, "../../data/5_hogwash/")

  thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv", stringsAsFactors = FALSE)
  thresholds <- thresholds %>% filter(shorthand == temp_shorthand)
  pval_thresh <- thresholds %>% pull(indiv_pval)
  epsilon_thresh <- thresholds %>% pull(indiv_epsilon)

  if (phenotypes$bin_cont[i] == "continuous") {
    hogwash_continuous <- data$hogwash_continuous
    data <- NULL
    
    # get size of genotype matrix and types of variants from output data ----
    indiv_pvalue_df <- data.frame(hogwash_continuous$hit_pvals)
    indiv_pvalue_df <- indiv_pvalue_df %>% 
      mutate("Variant type" = row.names(indiv_pvalue_df))
    indiv_pvalue_df$`Variant type`[grepl("Indel", indiv_pvalue_df$`Variant type`)] <- "Indel"
    indiv_pvalue_df$`Variant type`[grepl("SNP", indiv_pvalue_df$`Variant type`)] <- "SNP"
    indiv_pvalue_df$`Variant type`[!grepl("SNP", indiv_pvalue_df$`Variant type`) & 
                                     !grepl("Indel", indiv_pvalue_df$`Variant type`)] <- "Accessory Gene"
    variant_type_df <- indiv_pvalue_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Individual", 
             "bin_cont" = "continuous", 
             "Test" = "continuous", 
             "Input or output" = "output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_type_df)
    
    # Repeat for all input data ----
    indiv_input_df <- data.frame(colnames(genotype), row.names = colnames(genotype))
    indiv_input_df <- indiv_input_df %>% 
      mutate("Variant type" = row.names(indiv_input_df))
    
    indiv_input_df$`Variant type`[grepl("Indel", indiv_input_df$`Variant type`)] <- "Indel"
    indiv_input_df$`Variant type`[grepl("SNP", indiv_input_df$`Variant type`)] <- "SNP"
    indiv_input_df$`Variant type`[!grepl("SNP", indiv_input_df$`Variant type`) & 
                                     !grepl("Indel", indiv_input_df$`Variant type`)] <- "Accessory Gene"
    variant_input_df <- indiv_input_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Individual", 
             "bin_cont" = "continuous", 
             "Test" = "continuous", 
             "Input or output" = "input")
    var_type_summary_df <- rbind(var_type_summary_df, variant_input_df)
    # ----
    
    
    sig_index <- hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(pval_thresh) &
      hogwash_continuous$convergence$epsilon > epsilon_thresh
    
    if (!identical(length(hogwash_continuous$hit_pvals$fdr_corrected_pvals),
                   length(hogwash_continuous$convergence$epsilon))) {
      stop("length mismatch")
    }
    # Reorder genotype matrix to match order in hogwash_continuous$hit_pvals$fdr_corrected_pvals
    genotype <- genotype[, colnames(genotype) %in% row.names(hogwash_continuous$hit_pvals), drop = FALSE]
    genotype <- genotype[, match(row.names(hogwash_continuous$hit_pvals), colnames(genotype)), drop = FALSE]
    
    if (!identical(colnames(genotype), row.names(hogwash_continuous$hit_pvals))) {stop("mismatch A")}
    if (!identical(length(sig_index), ncol(genotype))) {stop("mismatch B")}
    
    sig_genotype <- genotype[, sig_index, drop = FALSE]
    
    write.table(x = sig_genotype, 
                quote = TRUE,
                file = paste0("../../data/1_combine_genotype_matrices/combinedmat/", 
                              temp_shorthand, 
                              "_sig_indiv_pan_and_rereferenced_snp_indel.tsv"),
                sep = "\t", 
                col.names = TRUE, 
                row.names = TRUE)
    
    # Repeat for all sig_genotype ----
    indiv_input_df <- data.frame(colnames(sig_genotype), row.names = colnames(sig_genotype))
    indiv_input_df <- indiv_input_df %>% 
      mutate("Variant type" = row.names(indiv_input_df))
    
    indiv_input_df$`Variant type`[grepl("Indel", indiv_input_df$`Variant type`)] <- "Indel"
    indiv_input_df$`Variant type`[grepl("SNP", indiv_input_df$`Variant type`)] <- "SNP"
    indiv_input_df$`Variant type`[!grepl("SNP", indiv_input_df$`Variant type`) & 
                                    !grepl("Indel", indiv_input_df$`Variant type`)] <- "Accessory Gene"
    variant_input_df <- indiv_input_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Individual", 
             "bin_cont" = "continuous", 
             "Test" = "continuous", 
             "Input or output" = "sig. output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_input_df)
    # ----
    
  } else {
    hogwash_phyc <- data$hogwash_phyc
    hogwash_synchronous <- data$hogwash_synchronous
    data <- NULL
    
    # get size of genotype matrix and types of variants ----
    indiv_pvalue_df <- data.frame(hogwash_synchronous$hit_pvals)
    indiv_pvalue_df <- indiv_pvalue_df %>% 
      mutate("Variant type" = row.names(indiv_pvalue_df))
    indiv_pvalue_df$`Variant type`[grepl("Indel", indiv_pvalue_df$`Variant type`)] <- "Indel"
    indiv_pvalue_df$`Variant type`[grepl("SNP", indiv_pvalue_df$`Variant type`)] <- "SNP"
    indiv_pvalue_df$`Variant type`[!grepl("SNP", indiv_pvalue_df$`Variant type`) & 
                                     !grepl("Indel", indiv_pvalue_df$`Variant type`)] <- "Accessory Gene"
    variant_type_df <- indiv_pvalue_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Individual", 
             "bin_cont" = "binary", 
             "Test" = "synchronous", 
             "Input or output" = "output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_type_df)
    # ----
    
    # Repeat for all input data ----
    indiv_input_df <- data.frame(colnames(genotype), row.names = colnames(genotype))
    indiv_input_df <- indiv_input_df %>% 
      mutate("Variant type" = row.names(indiv_input_df))
    indiv_input_df$`Variant type`[grepl("Indel", indiv_input_df$`Variant type`)] <- "Indel"
    indiv_input_df$`Variant type`[grepl("SNP", indiv_input_df$`Variant type`)] <- "SNP"
    indiv_input_df$`Variant type`[!grepl("SNP", indiv_input_df$`Variant type`) & 
                                    !grepl("Indel", indiv_input_df$`Variant type`)] <- "Accessory Gene"
    variant_input_df <- indiv_input_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Individual", 
             "bin_cont" = "binary", 
             "Test" = "synchronous", 
             "Input or output" = "input")
    var_type_summary_df <- rbind(var_type_summary_df, variant_input_df)
    variant_input_df2 <- variant_input_df
    variant_input_df2 <- variant_input_df2 %>% 
      mutate("Test" = "phyc")
    var_type_summary_df <- rbind(var_type_summary_df, variant_input_df2)
    
    # ------
    
    sig_index_sync <- hogwash_synchronous$hit_pvals$fdr_corrected_pvals > -log(pval_thresh) &
      hogwash_synchronous$convergence$epsilon > epsilon_thresh
    
    genotype <- genotype[, colnames(genotype) %in% row.names(hogwash_synchronous$hit_pvals), drop = FALSE]
    genotype <- genotype[, match(row.names(hogwash_synchronous$hit_pvals), colnames(genotype)), drop = FALSE]
    
    if (!identical(colnames(genotype), row.names(hogwash_synchronous$hit_pvals))) {stop("mismatch A")}
    if (!identical(length(sig_index_sync), ncol(genotype))) {stop("mismatch B")}
    
    sig_genotype_sync <- genotype[, sig_index_sync, drop = FALSE]
    
    write.table(x = sig_genotype_sync, 
                quote = TRUE,
                file = paste0("../../data/1_combine_genotype_matrices/combinedmat/", 
                              temp_shorthand, 
                              "_synchronous_sig_indiv_pan_and_rereferenced_snp_indel.tsv"),
                sep = "\t", 
                col.names = TRUE, 
                row.names = TRUE)
    
    
    # Repeat for significant output  ----
    indiv_input_df <- data.frame(colnames(sig_genotype_sync), row.names = colnames(sig_genotype_sync))
    indiv_input_df <- indiv_input_df %>% 
      mutate("Variant type" = row.names(indiv_input_df))
    indiv_input_df$`Variant type`[grepl("Indel", indiv_input_df$`Variant type`)] <- "Indel"
    indiv_input_df$`Variant type`[grepl("SNP", indiv_input_df$`Variant type`)] <- "SNP"
    indiv_input_df$`Variant type`[!grepl("SNP", indiv_input_df$`Variant type`) & 
                                    !grepl("Indel", indiv_input_df$`Variant type`)] <- "Accessory Gene"
    variant_input_df <- indiv_input_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Individual", 
             "bin_cont" = "binary", 
             "Test" = "synchronous", 
             "Input or output" = "sig. output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_input_df)
    # ---- 
    
    # Phyc
    # get size of genotype matrix and types of variants ----
    indiv_pvalue_df <- data.frame(hogwash_phyc$hit_pvals)
    indiv_pvalue_df <- indiv_pvalue_df %>% 
      mutate("Variant type" = row.names(indiv_pvalue_df))
    indiv_pvalue_df$`Variant type`[grepl("Indel", indiv_pvalue_df$`Variant type`)] <- "Indel"
    indiv_pvalue_df$`Variant type`[grepl("SNP", indiv_pvalue_df$`Variant type`)] <- "SNP"
    indiv_pvalue_df$`Variant type`[!grepl("SNP", indiv_pvalue_df$`Variant type`) & 
                                     !grepl("Indel", indiv_pvalue_df$`Variant type`)] <- "Accessory Gene"
    variant_type_df <- indiv_pvalue_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Individual", 
             "bin_cont" = "binary", 
             "Test" = "phyc", 
             "Input or output" = "output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_type_df)
    # ----
    
    sig_index_phyc <- hogwash_phyc$hit_pvals$fdr_corrected_pvals > -log(pval_thresh) &
      hogwash_phyc$convergence$epsilon > epsilon_thresh
    
    genotype <- genotype[, colnames(genotype) %in% row.names(hogwash_phyc$hit_pvals), drop = FALSE]
    genotype <- genotype[, match(row.names(hogwash_phyc$hit_pvals), colnames(genotype)), drop = FALSE]
    
    if (!identical(colnames(genotype), row.names(hogwash_phyc$hit_pvals))) {stop("mismatch A")}
    if (!identical(length(sig_index_phyc), ncol(genotype))) {stop("mismatch B")}
    
    sig_genotype_phyc <- genotype[, sig_index_phyc, drop = FALSE]

    write.table(x = sig_genotype_phyc, 
                quote = TRUE,
                file = paste0("../../data/1_combine_genotype_matrices/combinedmat/", 
                              temp_shorthand, 
                              "_phyc_sig_indiv_pan_and_rereferenced_snp_indel.tsv"),
                sep = "\t", 
                col.names = TRUE, 
                row.names = TRUE)
    
    # Repeat for significant output  ----
    indiv_input_df <- data.frame(colnames(sig_genotype_phyc), row.names = colnames(sig_genotype_phyc))
    indiv_input_df <- indiv_input_df %>% 
      mutate("Variant type" = row.names(indiv_input_df))
    indiv_input_df$`Variant type`[grepl("Indel", indiv_input_df$`Variant type`)] <- "Indel"
    indiv_input_df$`Variant type`[grepl("SNP", indiv_input_df$`Variant type`)] <- "SNP"
    indiv_input_df$`Variant type`[!grepl("SNP", indiv_input_df$`Variant type`) & 
                                    !grepl("Indel", indiv_input_df$`Variant type`)] <- "Accessory Gene"
    variant_input_df <- indiv_input_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Individual", 
             "bin_cont" = "binary", 
             "Test" = "phyc", 
             "Input or output" = "sig. output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_input_df)
    # ---- 
  }
}

write_tsv(var_type_summary_df, 
          "../../data/12_variant_types/variant_type_indiv.tsv", 
          col_names = TRUE)

