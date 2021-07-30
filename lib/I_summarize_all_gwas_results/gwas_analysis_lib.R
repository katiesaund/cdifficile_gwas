library(tidyverse)

# color palettes
cbf_1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

indiv_color_palette <- cbf_1[c(3, 4, 8)] # blue, green, pink
names(indiv_color_palette) <- c("Accessory", "Indel", "SNP")

group_color_palette <- cbf_1[c(2, 6)] # orange and navy
names(group_color_palette) <- c("Gene", "Intergenic region")

# Functions --------------------------------------------------------------------	
format_indiv_for_manhattan <- function(indiv_obj) {	
  indiv_pvalue_df <- data.frame(indiv_obj$hit_pvals)	
  
  # Assign genotype number to order x-axis by	
  indiv_pvalue_df <- cbind(1:nrow(indiv_pvalue_df), indiv_pvalue_df)	
  colnames(indiv_pvalue_df)[1] <- "Locus Significance"	
  
  # Classify genotypes by type	
  indiv_pvalue_df <- indiv_pvalue_df %>% 	
    mutate("Variant type" = row.names(indiv_pvalue_df))	
  indiv_pvalue_df$`Variant type`[grepl("Indel", indiv_pvalue_df$`Variant type`)] <- "Indel"	
  indiv_pvalue_df$`Variant type`[grepl("SNP", indiv_pvalue_df$`Variant type`)] <- "SNP"	
  indiv_pvalue_df$`Variant type`[!grepl("SNP", indiv_pvalue_df$`Variant type`) & !grepl("Indel", indiv_pvalue_df$`Variant type`)] <- "Accessory"	
  
  return(indiv_pvalue_df)	
}	

format_group_for_manhattan <- function(group_obj) {	
  group_pvalue_df <- data.frame(group_obj$hit_pvals)	
  
  # Assign genotype number to order x-axis by	
  group_pvalue_df <- cbind(1:nrow(group_pvalue_df), group_pvalue_df)	
  colnames(group_pvalue_df)[1] <- "Locus Significance"	
  
  # Classify genotypes by type	
  group_pvalue_df <- group_pvalue_df %>% 	
    mutate("Group type" = as.character(row.names(group_pvalue_df)))	
  group_pvalue_df$`Group type`[!grepl("-", group_pvalue_df$`Group type`)] <- "Gene"	
  group_pvalue_df$`Group type`[grepl("-", group_pvalue_df$`Group type`)] <- "Intergenic region"	
  return(group_pvalue_df)	
}	

format_indiv_for_pval_vs_epsilon <- function(indiv_obj){	
  indiv_pval_convergence <- 	
    cbind(indiv_obj$hit_pvals, 	
          indiv_obj$convergence$epsilon)	
  indiv_pval_convergence <- 	
    cbind(indiv_pval_convergence, 	
          row.names(indiv_pval_convergence))	
  colnames(indiv_pval_convergence) <- c("pval", "epsilon", "genotype")	
  
  indiv_pval_convergence <- indiv_pval_convergence %>% 	
    mutate("Variant type" = row.names(indiv_pval_convergence))	
  indiv_pval_convergence$`Variant type`[grepl("Indel", indiv_pval_convergence$`Variant type`)] <- "Indel"	
  indiv_pval_convergence$`Variant type`[grepl("SNP", indiv_pval_convergence$`Variant type`)] <- "SNP"	
  indiv_pval_convergence$`Variant type`[!grepl("SNP", indiv_pval_convergence$`Variant type`) & 	
                                          !grepl("Indel", indiv_pval_convergence$`Variant type`)] <- "Accessory"	
  return(indiv_pval_convergence)	
}	

format_group_for_pval_vs_epsilon <- function(group_obj){	
  group_pval_convergence <- 	
    cbind(group_obj$hit_pvals, 	
          group_obj$convergence$epsilon)	
  group_pval_convergence <- 	
    cbind(group_pval_convergence, 	
          row.names(group_pval_convergence))	
  colnames(group_pval_convergence) <- c("pval", "epsilon", "genotype")	
  
  group_pval_convergence <- group_pval_convergence %>% 	
    mutate("Group type" = as.character(row.names(group_pval_convergence)))	
  group_pval_convergence$`Group type`[!grepl("-", group_pval_convergence$`Group type`)] <- "Gene"	
  group_pval_convergence$`Group type`[grepl("-", group_pval_convergence$`Group type`)] <- "Intergenic region"	
  return(group_pval_convergence)	
}	

plot_indiv_manhattan <- function(indiv_pval_df, name_char, neg_log_pval){	
  indiv_pval_df %>% 	
    ggplot(mapping = aes(x = `Locus Significance`,	
                         y = fdr_corrected_pvals)) + 	
    geom_jitter(alpha = 0.25, aes(color = `Variant type`), height = 0) + 	
    xlab("Genetic loci (unordered)") +	
    ylab("-ln(FDR Corrected P-value)") + 	
    scale_color_manual(values = c(indiv_color_palette)) + 	
    theme_bw() + 	
    ggtitle(paste0(name_char, " individual loci")) + 
    geom_hline(yintercept = neg_log_pval, col = "red") + 
    ylim(c(0, max(indiv_pval_df$fdr_corrected_pvals, neg_log_pval)))
}	

plot_indiv_pval_vs_epsilon <- function(indiv_p_ep_df, name_char, neg_log_pval, ep_num){	
  indiv_p_ep_df %>% 	
    ggplot(aes(x = epsilon, y = pval)) + 	
    geom_jitter(alpha = 0.15, aes(color = `Variant type`), height = 0) + 	
    theme_bw() + 	
    xlab("Convergence metric") + 	
    ylab("-ln(FDR corrected P-value)") +	
    ggtitle(paste0(name_char, " individual loci")) +	
    scale_color_manual(values = c(indiv_color_palette))	+ 
    geom_hline(yintercept = neg_log_pval, col = "red") +
    geom_vline(xintercept = ep_num, col = "red") + 
    ylim(c(0, max(indiv_p_ep_df$pval, neg_log_pval))) + 
    xlim(c(0, max(indiv_p_ep_df$epsilon, ep_num)))
}	

plot_group_manhattan <- function(group_pval_df, name_char, neg_log_pval){	
  group_pval_df %>% 	
    ggplot(mapping = aes(x = `Locus Significance`,	
                         y = fdr_corrected_pvals)) + 	
    geom_jitter(alpha = 0.25, aes(color = `Group type`), height = 0) + 	
    xlab("Genetic loci (unordered)") +	
    ylab("-ln(FDR Corrected P-value)") + 	
    scale_color_manual(values = c(group_color_palette)) + 	
    theme_bw() + 	
    ggtitle(paste0(name_char, " grouped loci"))	+ 
    geom_hline(yintercept = neg_log_pval, col = "red") + 
    ylim(c(0, max(group_pval_df$fdr_corrected_pvals, neg_log_pval)))
}	

plot_group_pval_vs_epsilon <- function(group_p_ep_df, name_char, neg_log_pval, ep_num){	
  group_p_ep_df %>% 	
    ggplot(aes(x = epsilon, y = pval)) + 	
    geom_jitter(alpha = 0.15, aes(col = `Group type`), height = 0) + 	
    theme_bw() + 	
    xlab("Convergence metric") + 	
    ylab("-ln(FDR corrected P-value)") +	
    ggtitle(paste0(name_char, " grouped loci")) +	
    scale_color_manual(values = c(group_color_palette))	+ 
    geom_hline(yintercept = neg_log_pval, col = "red") +
    geom_vline(xintercept = ep_num, col = "red") + 
    ylim(c(0, max(group_p_ep_df$pval, neg_log_pval))) + 
    xlim(c(0, max(group_p_ep_df$epsilon, ep_num)))
}

load_indiv_hogwash_rda_for_pheno <- function(pheno_shorthand, path_to_hog_data) {
  if (pheno_shorthand %in% c("fqR", "severity")) {
    load(paste0(path_to_hog_data, pheno_shorthand, "/hogwash_phyc_phyc_", pheno_shorthand, ".rda"))
    load(paste0(path_to_hog_data, pheno_shorthand, "/hogwash_synchronous_synchronous_", pheno_shorthand, ".rda"))
    results <- list("hogwash_phyc" = hogwash_phyc, "hogwash_synchronous" = hogwash_synchronous)
  } else if (pheno_shorthand %in% c("log_cfe", "log_germ_tc_and_gly",
                                    "log_germ_tc", "log_growth",
                                    "log_sporulation", "log_toxin")) {
    load(paste0(path_to_hog_data, pheno_shorthand, "/hogwash_continuous_continuous_", pheno_shorthand, ".rda"))
    results <- list("hogwash_continuous" = hogwash_continuous)
  } else {
    results <- NULL
    stop("Phenotype name not in set of Hanna in vitro phenotypes")
  }
  return(results)
}

load_post_ar_hogwash_rda_for_pheno <- function(pheno_shorthand, path_to_hog_data) {
  if (pheno_shorthand %in% c("fqR", "severity")) {
    load(paste0(path_to_hog_data, pheno_shorthand, "/hogwash_phyc_grouped_phyc_", pheno_shorthand, "_grouped_post-ar.rda"))
    load(paste0(path_to_hog_data, pheno_shorthand, "/hogwash_synchronous_grouped_synchronous_", pheno_shorthand, "_grouped_post-ar.rda"))
    results <- list("hogwash_phyc" = hogwash_phyc, "hogwash_synchronous" = hogwash_synchronous)
  } else if (pheno_shorthand %in% c("log_cfe", "log_germ_tc_and_gly",
                                    "log_germ_tc", "log_growth",
                                    "log_sporulation", "log_toxin")) {
    load(paste0(path_to_hog_data, pheno_shorthand, "/hogwash_continuous_grouped_continuous_", pheno_shorthand, "_grouped_post-ar.rda"))
    results <- list("hogwash_continuous" = hogwash_continuous)
  } else {
    results <- NULL
    stop("Phenotype name not in set of Hanna in vitro phenotypes")
  }
  return(results)
}

rank_by_epsilon_and_p <- function(hogwash_obj) {
  df <- data.frame(row.names(hogwash_obj$hit_pvals), 
                   hogwash_obj$hit_pvals$fdr_corrected_pvals, 
                   hogwash_obj$convergence$epsilon)
  colnames(df) <- c("Genotype", "-ln(FDR P-value)", "Epsilon")
  df <- df[order(-df$Epsilon, -df$`-ln(FDR P-value)`), ]
  row.names(df) <- NULL
  df <- cbind(df, 1:nrow(df))
  colnames(df) <- c("Genotype", "-ln(FDR P-value)", "Epsilon", "Rank")
  return(df)
}

plot_ribo_tree_and_heatmap <- function(pheno_short_name, test_name, group_or_indiv, phenotypes) {
  temp_short_hand <- pheno_short_name
  temp_long_hand <- phenotypes %>% filter(shorthand == pheno_short_name) %>% pull(longhand) 
  temp_type <- phenotypes %>% filter(shorthand == pheno_short_name) %>% pull(bin_cont)
  tree <- read.tree(paste0(data_dir, "3_trees/", temp_short_hand, ".tree"))
  ribotype <- read.table(paste0(data_dir, "6_ribotype/", temp_short_hand, "_ribotype.tsv"), 
                         sep = "\t", 
                         stringsAsFactors = FALSE)       
  phenotype <- read.table(paste0(data_dir, "4_phenotypes/cleaned_", temp_short_hand, ".tsv"),
                          sep = "\t",
                          row.names = 1,
                          header = TRUE,
                          stringsAsFactors = FALSE,
                          check.names = FALSE)
  phenotype <- as.matrix(phenotype)
  
  
  if (test_name == "continuous") {
    sig_genotype <- data.table::fread(file = paste0(data_dir, "1_combine_genotype_matrices/combinedmat/", temp_short_hand, "_sig_", group_or_indiv, "_pan_and_rereferenced_snp_indel.tsv"),
                                      sep = "\t", 
                                      data.table = FALSE, 
                                      stringsAsFactors = FALSE, 
                                      check.names = FALSE)
    
  } else if (test_name == "phyc") {
    sig_genotype <- data.table::fread(file = paste0(data_dir, "1_combine_genotype_matrices/combinedmat/", temp_short_hand, "_phyc_sig_", group_or_indiv, "_pan_and_rereferenced_snp_indel.tsv"),
                                      sep = "\t", 
                                      data.table = FALSE, 
                                      stringsAsFactors = FALSE, 
                                      check.names = FALSE)
  } else {
    sig_genotype <- data.table::fread(file = paste0(data_dir, "1_combine_genotype_matrices/combinedmat/", temp_short_hand, "_synchronous_sig_", group_or_indiv, "_pan_and_rereferenced_snp_indel.tsv"),
                                      sep = "\t",
                                      data.table = FALSE,
                                      stringsAsFactors = FALSE,
                                      check.names = FALSE)
  }
  
  sig_genotype <- as.matrix(sig_genotype)
  row.names(sig_genotype) <- sig_genotype[, 1]
  sig_genotype <- sig_genotype[, 2:ncol(sig_genotype)]
  storage.mode(sig_genotype) <- "numeric"
  
  if (!is.null(dim(sig_genotype))) {
    ribotype <- as.data.frame(cbind(row.names(ribotype), ribotype), 
                              stringsAsFactors = FALSE)
    colnames(ribotype) <- c("ID", "ribotype", "ribotype_color")
    
    phenotype <- as.data.frame(cbind(row.names(phenotype), phenotype))
    colnames(phenotype) <- c("ID", "pheno")
    phenotype$ID <- as.character(phenotype$ID)
    phenotype$pheno <- as.numeric(as.character(phenotype$pheno))
    
    ribotype_colors <- unique(ribotype$ribotype_color)
    names(ribotype_colors) <- unique(ribotype$ribotype)
    
    p <- ggtree(tree) %<+% ribotype + 
      geom_tippoint(aes(color = ribotype)) + 
      scale_color_manual(values = ribotype_colors)
    
    heatmap.colours <- c("white","black")
    names(heatmap.colours) <- 0:1
    
    bar_plot_p1 <- ggplot(phenotype, aes(ID, pheno)) + 
      geom_col() + 
      coord_flip() + 
      ggtitle(temp_long_hand) + 
      theme_tree2()
    
    d <- dist(t(sig_genotype), method = "euclidean") # distance matrix
    fit <- hclust(d, method = "ward.D")
    
    id_df <- as.data.frame(row.names(sig_genotype), stringsAsFactors = FALSE)
    colnames(id_df)[1] <- "ID"
    
    sig_genotype_df <- as.data.frame(sig_genotype, stringsAsFactors = FALSE)
    row.names(sig_genotype_df) <- NULL
    sig_genotype_df <- as.data.frame(sapply(sig_genotype_df, as.numeric)) #<- sapply is here
    sig_genotype_df <- cbind(id_df, sig_genotype_df)
    
    sig_genotype_df_long <- sig_genotype_df %>%  
      pivot_longer(-ID, names_to = "locus", values_to = "presence")
    
    sig_genotype_df_long$ID <- factor(sig_genotype_df_long$ID, levels = rownames(sig_genotype), labels =  rownames(sig_genotype))
    sig_genotype_df_long$locus <- factor( sig_genotype_df_long$locus, levels = colnames(sig_genotype)[fit$order],  labels = colnames(sig_genotype)[fit$order])
    
    heatmap_plot_p3 <- ggplot(sig_genotype_df_long, aes(x = locus, y = ID)) + 
      geom_tile(aes(fill = factor(presence))) + 
      scale_fill_manual(values = c("white", "black")) + 
     # ggtitle(paste0(group_or_indiv, ". geno. presence")) + 
      guides(fill = FALSE) + 
      theme_void()
    
    png(paste0(fig_dir, 
               temp_short_hand,
               "_", 
               test_name, 
               "_", 
               group_or_indiv, 
               "_results_tree_w_ribo_sig_pval_hi_epsilon_hits_heatmap_and_pheno_barplot.png"), 
        units = "in", 
        width = 9, height = 6, res = 250)
    print(heatmap_plot_p3 %>% insert_left(p) %>% insert_right(bar_plot_p1))
    dev.off()
    tree_htmp_pheno_plot <- heatmap_plot_p3 %>%
      insert_left(p) %>%
      insert_right(bar_plot_p1)
    save(tree_htmp_pheno_plot,
         file = paste0(fig_dir, 
                       temp_short_hand,
                       "_", 
                       test_name, 
                       "_", 
                       group_or_indiv, 
                       "_results_tree_w_ribo_sig_pval_hi_epsilon_hits_heatmap_and_pheno_barplot.RData"))
    }
}
