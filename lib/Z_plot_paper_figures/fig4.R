# A) Toxin gwas p-ep plot
# B) tcdb hits / loci
# C) indiv. signifant hits with heatmap and toxin phenotype 
# D) group results 

library(aplot)
library(patchwork)
library(ggplotify)

#indiv_plot
# Made in I_summarize_all_gwas_results/02_plot_manhattan_and_pval_vs_epsilon
load("../../figures/gwas_results/log_toxin_indiv_pval_vs_convergence.RData")
new_indiv_plot <- indiv_plot +
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold")) + 
  ggtitle("") + 
  ylab("-ln(P-value)") +
  xlim(c(0, 0.20)) + 
  ylim(c(0, 2.25)) 

#tree_htmp_pheno_plot
load("../../figures/gwas_results/log_toxin_continuous_indiv_results_tree_w_ribo_sig_pval_hi_epsilon_hits_heatmap_and_pheno_barplot.RData")
gg_tree_htmp_pheno_plot <- as.ggplot(tree_htmp_pheno_plot)
new_tree_htmp_pheno_plot <- gg_tree_htmp_pheno_plot +
  labs(tag = "C") +
  theme(plot.tag = element_text(face = "bold")) + 
  ggtitle("")

#group_plot
load("../../figures/gwas_results/log_toxin_group_pval_vs_convergence.RData")
new_group_plot <- group_plot +
  labs(tag = "D") +
  theme(plot.tag = element_text(face = "bold")) + 
  ggtitle("") + 
  xlim(c(0, 0.20)) + 
  ylim(c(0, 2.25)) + 
  ylab("-ln(P-value)")

# Made in J_tcdb_enrichment/06_plot_paloc_figure.R
# pheatmap_plot
load("../../figures/paloc_name_tbd_3.RData")
gg_pheatmap_plot <- as.ggplot(pheatmap_plot)
new_gg_pheatmap_plot <- gg_pheatmap_plot + 
  labs(tag = "B") + 
  theme(plot.tag = element_text(face = "bold"))

# method 1
# final_plot <- 
#   new_indiv_plot +
#   new_tree_htmp_pheno_plot +
#   new_group_plot +
#   new_gg_pheatmap_plot +
#   plot_layout(ncol = 2) 

# method 2
# (new_indiv_plot + new_tree_htmp_pheno_plot + new_group_plot) / (new_gg_pheatmap_plot)

# method 3
# layout <- '
# ABC
# D
# '
# wrap_plots(A = new_indiv_plot,
#            B = new_tree_htmp_pheno_plot, 
#            C = new_group_plot, 
#            D = new_gg_pheatmap_plot, 
#            design = layout)

# method 4
layout <- c(
  area(1, 1), 
  area(1, 2, 1, 3), 
  area(1, 4), 
  area(2, 1, 2, 4)
)

# Check out layout design
# plot(layout)

new_indiv_plot + new_tree_htmp_pheno_plot + new_group_plot + new_gg_pheatmap_plot + plot_layout(design = layout)

ggsave(filename = "../../figures/Fig4.png",
       units = "in",
       height = 12,
       width = 14)
ggsave(filename = "../../figures/Fig4.pdf",
       units = "in",
       height = 12,
       width = 14) 
