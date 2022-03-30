# Manhattan Plot and P-value vs. Epsilon (convergence) for each phenotype's 
# hogwash results (continuous, phyc, and synchronous & grouped and individual)
source("gwas_analysis_lib.R")
library(tidyverse)
# Make all plots ---------------------------------------------------------------
hogwash_data_dir <- "../../data/5_hogwash/"
fig_dir <- "../../figures/gwas_results/"

phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")
binary_phenotypes <- continuous_phenotypes <- phenotypes
colnames(binary_phenotypes)[1] <- colnames(continuous_phenotypes)[1] <- 'name'
binary_phenotypes <- binary_phenotypes %>% filter(bin_cont == "binary")
continuous_phenotypes <- continuous_phenotypes %>% filter(bin_cont == "continuous")
thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv", stringsAsFactors = FALSE)

# Binary phenotypes -----
for (i in 1:2) {
    id <- binary_phenotypes$name[i]
    temp_thresh <- thresholds %>% filter(shorthand == id)

    
    # Individual ---------------------------------------------------------------
    pval_thresh <- temp_thresh %>% pull(indiv_pval)
    epsilon_thresh <- temp_thresh %>% pull(indiv_epsilon)
    
    load(paste0(hogwash_data_dir,id, "/hogwash_phyc_phyc_", id, ".rda"))
    indiv_phyc <- hogwash_phyc
    hogwash_phyc <- NULL
    
    indiv_phyc_pval_df <- format_indiv_for_manhattan(indiv_phyc)
    plot_indiv_manhattan(indiv_phyc_pval_df, id, -log(pval_thresh))
    ggsave(paste0(fig_dir, id, "_indiv_phyc_manhattan.png"), 
           units = "in",
           height = 5,
           width = 5) 
    
    indiv_pval_ep_df <- format_indiv_for_pval_vs_epsilon(indiv_phyc)
    indiv_plot <- plot_indiv_pval_vs_epsilon(indiv_pval_ep_df, id, -log(pval_thresh), epsilon_thresh)
    ggsave(plot = indiv_plot, 
           filename = paste0(fig_dir, id, "_indiv_phyc_pval_vs_convergence.png"), 
           units = "in",
           height = 5,
           width = 5)
    save(indiv_plot,
         file = paste0(fig_dir, id, "_indiv_phyc_pval_vs_convergence.RData"))
    
    load(paste0(hogwash_data_dir, id, "/hogwash_synchronous_synchronous_", id, ".rda"))
    indiv_sync <- hogwash_synchronous
    hogwash_synchronous <- NULL
    
    indiv_sync_pval_df <- format_indiv_for_manhattan(indiv_sync)
    plot_indiv_manhattan(indiv_sync_pval_df, id, -log(pval_thresh))
    ggsave(paste0(fig_dir, id, "_indiv_sync_manhattan.png"), 
           units = "in",
           height = 5,
           width = 5) 
    
    indiv_sync_ep_df <- format_indiv_for_pval_vs_epsilon(indiv_sync)
    indiv_plot <- plot_indiv_pval_vs_epsilon(indiv_sync_ep_df, id, -log(pval_thresh), epsilon_thresh)
    ggsave(plot = indiv_plot, 
           filename = paste0(fig_dir, id, "_indiv_sync_pval_vs_convergence.png"), 
           units = "in",
           height = 5,
           width = 5) 
    save(indiv_plot,
         file = paste0(fig_dir, id, "_indiv_sync_pval_vs_convergence.RData"))
    
    
    # Grouped ------------------------------------------------------------------
    # Grouped PhyC ----
    pval_thresh <- temp_thresh %>% pull(group_pval)
    epsilon_thresh <- temp_thresh %>% pull(group_epsilon)
    
    load(paste0(hogwash_data_dir,id, "/hogwash_phyc_grouped_phyc_", id, "_grouped_post-ar.rda"))
    group_phyc <- hogwash_phyc
    hogwash_phyc <- NULL
    
    group_pval_df <- format_group_for_manhattan(group_phyc)
    plot_group_manhattan(group_pval_df, id, -log(pval_thresh))
    ggsave(paste0(fig_dir, id, "_group_phyc_manhattan.png"), 
           units = "in",
           height = 5,
           width = 5)  
    
    group_pval_epsilon_df <- format_group_for_pval_vs_epsilon(group_phyc)
    group_plot <- plot_group_pval_vs_epsilon(group_pval_epsilon_df, id, -log(pval_thresh), epsilon_thresh)
    ggsave(plot <- group_plot, 
           filename = paste0(fig_dir, id, "_group_phyc_pval_vs_convergence.png"), 
           units = "in",
           height = 5,
           width = 5)
    save(group_plot,
         file = paste0(fig_dir, id, "_group_phyc_pval_vs_convergence.RData"))
    
    # Grouped Synchronous ---- 
    load(paste0(hogwash_data_dir,id, "/hogwash_synchronous_grouped_synchronous_", id, "_grouped_post-ar.rda"))
    group_sync <- hogwash_synchronous
    hogwash_synchronous <- NULL
    
    group_pval_df <- format_group_for_manhattan(group_sync)
    plot_group_manhattan(group_pval_df, id, -log(pval_thresh))
    ggsave(paste0(fig_dir, id, "_group_sync_manhattan.png"), 
           units = "in",
           height = 5,
           width = 5)  
    
    group_pval_epsilon_df <- format_group_for_pval_vs_epsilon(group_sync)
    group_plot <- plot_group_pval_vs_epsilon(group_pval_epsilon_df, id, -log(pval_thresh), epsilon_thresh)
    ggsave(plot = group_plot, 
           filename = paste0(fig_dir, id, "_group_sync_pval_vs_convergence.png"), 
           units = "in",
           height = 5,
           width = 5) 
    save(group_plot,
         file = paste0(fig_dir, id, "_group_sync_pval_vs_convergence.RData"))
}

# Continuous phenotypes --------------------------------------------------------
for (i in 1:6) {
  id <- continuous_phenotypes$name[i]
  load(paste0(hogwash_data_dir, id, "/hogwash_continuous_continuous_", id, ".rda"))
  indiv_obj <- hogwash_continuous
  hogwash_continuous <- NULL
  temp_thresh <- thresholds %>% filter(shorthand == id)
  
  # Individual ---------------------------------------------------------------
  pval_thresh <- temp_thresh %>% pull(indiv_pval)
  epsilon_thresh <- temp_thresh %>% pull(indiv_epsilon)  
  indiv_cont_pval_df <- format_indiv_for_manhattan(indiv_obj)
  plot_indiv_manhattan(indiv_cont_pval_df, id, -log(pval_thresh))
  ggsave(paste0(fig_dir, id, "_indiv_manhattan.png"), 
         units = "in",
         height = 5,
         width = 5) 
  
  indiv_cont_ep_df <- format_indiv_for_pval_vs_epsilon(indiv_obj)
  indiv_plot <- plot_indiv_pval_vs_epsilon(indiv_cont_ep_df, id, -log(pval_thresh), epsilon_thresh)
  ggsave(plot = indiv_plot, 
         filename = paste0(fig_dir, id, "_indiv_pval_vs_convergence.png"), 
         units = "in",
         height = 5,
         width = 5) 
  save(indiv_plot,
       file = paste0(fig_dir, id, "_indiv_pval_vs_convergence.RData"))
  
 
  # Grouped results ------------------------------------------------------------
  pval_thresh <- temp_thresh %>% pull(group_pval)
  epsilon_thresh <- temp_thresh %>% pull(group_epsilon)
  
  load(paste0(hogwash_data_dir, id, "/hogwash_continuous_grouped_continuous_", id, "_grouped_post-ar.rda"))
  group_obj <- hogwash_continuous
  hogwash_continuous <- NULL
  
  group_pval_df <- format_group_for_manhattan(group_obj)
  plot_group_manhattan(group_pval_df, id, -log(pval_thresh))
  ggsave(paste0(fig_dir, id, "_group_manhattan.png"), 
         units = "in",
         height = 5,
         width = 5)  
  
  group_pval_epsilon_df <- format_group_for_pval_vs_epsilon(group_obj)
  group_plot <- plot_group_pval_vs_epsilon(group_pval_epsilon_df, id, -log(pval_thresh), epsilon_thresh)
  ggsave(plot = group_plot,
         filename = paste0(fig_dir, id, "_group_pval_vs_convergence.png"), 
         units = "in",
         height = 5,
         width = 5) 
  save(group_plot,
       file = paste0(fig_dir, id, "_group_pval_vs_convergence.RData"))
}
