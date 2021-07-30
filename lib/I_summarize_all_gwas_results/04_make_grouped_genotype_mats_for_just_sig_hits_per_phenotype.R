# Goals:
#   - Make a matrix that contains just those genotypes that are of interest: by significance and convergence. 
#   - Record some summary statistics (since I'm loading in this huge geno mat anyways): 
#       - # of variants input into hogwash (and type) 
#       - # of variants in final output (and type)
#       - # of variant that are significant (and type)

source("gwas_analysis_lib.R")
phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")

# Grouped loci
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
  
  key <- read.table(paste0("../../data/2_gene_keys/cleaned_", 
                    temp_shorthand,
                    "_gene_and_intergenic_key.tsv"),
                    sep = "\t",
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
  key <- as.matrix(key)
  
  # Make geno mat where we collapse variants into groups ----
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
  # -----
  data <- load_post_ar_hogwash_rda_for_pheno(temp_shorthand, "../../data/5_hogwash/")

  thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv", stringsAsFactors = FALSE)
  thresholds <- thresholds %>% filter(shorthand == temp_shorthand)
  pval_thresh <- thresholds %>% pull(group_pval)
  epsilon_thresh <- thresholds %>% pull(group_epsilon)
  
  if (phenotypes$bin_cont[i] == "continuous") {
    hogwash_continuous <- data$hogwash_continuous
    data <- NULL
    
    # get size of genotype matrix and types of variants from output data ----
    group_pvalue_df <- data.frame(hogwash_continuous$hit_pvals)
    group_pvalue_df <- group_pvalue_df %>% 
      mutate("Variant type" = row.names(group_pvalue_df))
    group_pvalue_df$`Variant type`[!grepl("-", group_pvalue_df$`Variant type`)] <- "Gene"
    group_pvalue_df$`Variant type`[grepl("-", group_pvalue_df$`Variant type`)] <- "Intergenic region"
    
    variant_type_df <- group_pvalue_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Grouped", 
             "bin_cont" = "continuous", 
             "Test" = "continuous", 
             "Input or output" = "output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_type_df)
    
    # Repeat for all input data ----
    group_input_df <- data.frame(colnames(samples_by_genes), row.names = colnames(samples_by_genes))
    group_input_df <- group_input_df %>% 
      mutate("Variant type" = row.names(group_input_df))
    group_input_df$`Variant type`[!grepl("-", group_input_df$`Variant type`)] <- "Gene"
    group_input_df$`Variant type`[grepl("-", group_input_df$`Variant type`)] <- "Intergenic region"
    variant_input_df <- group_input_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Grouped", 
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
    # Reorder samples_by_genes matrix to match order in hogwash_continuous$hit_pvals$fdr_corrected_pvals
    samples_by_genes <- samples_by_genes[, colnames(samples_by_genes) %in% row.names(hogwash_continuous$hit_pvals), drop = FALSE]
    samples_by_genes <- samples_by_genes[, match(row.names(hogwash_continuous$hit_pvals), colnames(samples_by_genes)), drop = FALSE]
    
    if (!identical(colnames(samples_by_genes), row.names(hogwash_continuous$hit_pvals))) {stop("mismatch A")}
    if (!identical(length(sig_index), ncol(samples_by_genes))) {stop("mismatch B")}
    
    sig_genotype <- samples_by_genes[, sig_index, drop = FALSE]
    
    write.table(x = sig_genotype, 
                quote = TRUE,
                file = paste0("../../data/1_combine_genotype_matrices/combinedmat/", 
                              temp_shorthand, 
                              "_sig_group_pan_and_rereferenced_snp_indel.tsv"),
                sep = "\t", 
                col.names = TRUE, 
                row.names = TRUE)
    
    # Repeat for all sig_genotype ----
    group_input_df <- data.frame(colnames(sig_genotype), row.names = colnames(sig_genotype))
    group_input_df <- group_input_df %>% 
      mutate("Variant type" = row.names(group_input_df))
    group_input_df$`Variant type`[!grepl("-", group_input_df$`Variant type`)] <- "Gene"
    group_input_df$`Variant type`[grepl("-", group_input_df$`Variant type`)] <- "Intergenic region"
    variant_input_df <- group_input_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Grouped", 
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
    group_pvalue_df <- data.frame(hogwash_synchronous$hit_pvals)
    group_pvalue_df <- group_pvalue_df %>% 
      mutate("Variant type" = row.names(group_pvalue_df))
    group_pvalue_df$`Variant type`[!grepl("-", group_pvalue_df$`Variant type`)] <- "Gene"
    group_pvalue_df$`Variant type`[grepl("-", group_pvalue_df$`Variant type`)] <- "Intergenic region"
    variant_type_df <- group_pvalue_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Grouped", 
             "bin_cont" = "binary", 
             "Test" = "synchronous", 
             "Input or output" = "output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_type_df)
    # ----
    
    # Repeat for all input data ----
    group_input_df <- data.frame(colnames(samples_by_genes), row.names = colnames(samples_by_genes))
    group_input_df <- group_input_df %>% 
      mutate("Variant type" = row.names(group_input_df))
    group_input_df$`Variant type`[!grepl("-", group_input_df$`Variant type`)] <- "Gene"
    group_input_df$`Variant type`[grepl("-", group_input_df$`Variant type`)] <- "Intergenic region"
    variant_input_df <- group_input_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Grouped", 
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
    
    samples_by_genes <- samples_by_genes[, colnames(samples_by_genes) %in% row.names(hogwash_synchronous$hit_pvals), drop = FALSE]
    samples_by_genes <- samples_by_genes[, match(row.names(hogwash_synchronous$hit_pvals), colnames(samples_by_genes)), drop = FALSE]
    
    if (!identical(colnames(samples_by_genes), row.names(hogwash_synchronous$hit_pvals))) {stop("mismatch A")}
    if (!identical(length(sig_index_sync), ncol(samples_by_genes))) {stop("mismatch B")}
    
    sig_genotype_sync <- samples_by_genes[, sig_index_sync, drop = FALSE]
    
    write.table(x = sig_genotype_sync, 
                quote = TRUE,
                file = paste0("../../data/1_combine_genotype_matrices/combinedmat/", 
                              temp_shorthand, 
                              "_synchronous_sig_group_pan_and_rereferenced_snp_indel.tsv"),
                sep = "\t", 
                col.names = TRUE, 
                row.names = TRUE)
    
    
    # Repeat for significant output  ----
    group_input_df <- data.frame(colnames(sig_genotype_sync), row.names = colnames(sig_genotype_sync))
    group_input_df <- group_input_df %>% 
      mutate("Variant type" = row.names(group_input_df))
    group_input_df$`Variant type`[!grepl("-", group_input_df$`Variant type`)] <- "Gene"
    group_input_df$`Variant type`[grepl("-", group_input_df$`Variant type`)] <- "Intergenic region"
    
    variant_input_df <- group_input_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Grouped", 
             "bin_cont" = "binary", 
             "Test" = "synchronous", 
             "Input or output" = "sig. output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_input_df)
    # ---- 
    
    # Phyc
    # get size of genotype matrix and types of variants ----
    group_pvalue_df <- data.frame(hogwash_phyc$hit_pvals)
    group_pvalue_df <- group_pvalue_df %>% 
      mutate("Variant type" = row.names(group_pvalue_df))
    group_pvalue_df$`Variant type`[!grepl("-", group_pvalue_df$`Variant type`)] <- "Gene"
    group_pvalue_df$`Variant type`[grepl("-", group_pvalue_df$`Variant type`)] <- "Intergenic region"
    variant_type_df <- group_pvalue_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Grouped", 
             "bin_cont" = "binary", 
             "Test" = "phyc", 
             "Input or output" = "output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_type_df)
    # ----
    
    sig_index_phyc <- hogwash_phyc$hit_pvals$fdr_corrected_pvals > -log(pval_thresh) &
      hogwash_phyc$convergence$epsilon > epsilon_thresh
    
    samples_by_genes <- samples_by_genes[, colnames(samples_by_genes) %in% row.names(hogwash_phyc$hit_pvals), drop = FALSE]
    samples_by_genes <- samples_by_genes[, match(row.names(hogwash_phyc$hit_pvals), colnames(samples_by_genes)), drop = FALSE]
    
    if (!identical(colnames(samples_by_genes), row.names(hogwash_phyc$hit_pvals))) {stop("mismatch A")}
    if (!identical(length(sig_index_phyc), ncol(samples_by_genes))) {stop("mismatch B")}
    
    sig_genotype_phyc <- samples_by_genes[, sig_index_phyc, drop = FALSE]
    
    write.table(x = sig_genotype_phyc, 
                quote = TRUE,
                file = paste0("../../data/1_combine_genotype_matrices/combinedmat/", 
                              temp_shorthand, 
                              "_phyc_sig_group_pan_and_rereferenced_snp_indel.tsv"),
                sep = "\t", 
                col.names = TRUE, 
                row.names = TRUE)
    
    # Repeat for significant output  ----
    group_input_df <- data.frame(colnames(sig_genotype_phyc), row.names = colnames(sig_genotype_phyc))
    group_input_df <- group_input_df %>% 
      mutate("Variant type" = row.names(group_input_df))
    group_input_df$`Variant type`[!grepl("-", group_input_df$`Variant type`)] <- "Gene"
    group_input_df$`Variant type`[grepl("-", group_input_df$`Variant type`)] <- "Intergenic region"
    variant_input_df <- group_input_df %>%  count(`Variant type`) %>% 
      mutate("shorthand" = temp_shorthand,
             "Indiv. or Group." = "Grouped", 
             "bin_cont" = "binary", 
             "Test" = "phyc", 
             "Input or output" = "sig. output")
    var_type_summary_df <- rbind(var_type_summary_df, variant_input_df)
    # ---- 
  }
}

write_tsv(var_type_summary_df, 
          "../../data/12_variant_types/variant_type_group.tsv", 
          col_names = TRUE)

