# A) Toxin gwas p-ep plot
# B) tcdb hits / loci
# C) indiv. signifant hits with heatmap and toxin phenotype 
# D) group results 

library(aplot)
library(patchwork)
library(ggplotify)
library(ggplot2)
library(tidyverse)
library(ape)
library(ggtree)
source("../I_summarize_all_gwas_results/gwas_analysis_lib.R")
data_dir <- "../../data/"
# cytoxicity == log_toxin == ln(toxin activity)
hogwash_data_dir <- "../../data/5_hogwash/"
thresholds <- read.csv("../../data/7_thresholds/phenotype_thresholds.csv", stringsAsFactors = FALSE)
id <- "log_toxin"
load(paste0(hogwash_data_dir, id, "/hogwash_continuous_continuous_", id, ".rda"))
indiv_obj <- hogwash_continuous
hogwash_continuous <- NULL
temp_thresh <- thresholds %>% filter(shorthand == id)

# Individual ---------------------------------------------------------------
pval_thresh <- temp_thresh %>% pull(indiv_pval)
epsilon_thresh <- temp_thresh %>% pull(indiv_epsilon)  
indiv_cont_ep_df <- format_indiv_for_pval_vs_epsilon(indiv_obj)


verbose_names <- indiv_cont_ep_df$genotype
verbose_names[5278]
indiv_cont_ep_df$pretty_genotype <- verbose_names %>% 
  gsub("Coding ", "", .) %>% 
  gsub(" at .*locus_tag=.*Strand.*c[.]", " ", .) %>% 
  gsub("[|][0-9].*", "", .) %>% 
  gsub("[|]p[.]", " ", .) %>% 
  gsub("[|].*", "", .) %>% 
  gsub("fs", "frameshift", .)

indiv_cont_ep_df$pretty_genotype

gene_name <- verbose_names %>% gsub(".*[0-9][|]", "", .) %>% gsub(";", "", .) %>% gsub("[|].*", "", .)

for (i in 1:nrow(indiv_cont_ep_df)) {
  if (nchar(indiv_cont_ep_df$pretty_genotype[i]) > 10) {
    indiv_cont_ep_df$pretty_genotype[i] <- paste(gene_name[i], indiv_cont_ep_df$pretty_genotype[i])
  }
}

indiv_cont_ep_df  <- 
  indiv_cont_ep_df %>% 
  mutate(even_shorter_name = gsub("Indel.*", "Indel", pretty_genotype))

indiv_cont_ep_df %>% 
  filter(epsilon > 0.15 & pval > 2)

new_indiv_plot <- 
  indiv_cont_ep_df %>% 	
  ggplot(aes(x = epsilon, y = pval)) + 	
  geom_jitter(alpha = 0.50, aes(color = `Variant type`), height = 0) + 	
  theme_bw() + 	
  xlab(expression(epsilon)) + 	
  ylab("-ln(FDR corrected P-value)") +	
  scale_color_manual(values = c(indiv_color_palette))	+ 
  geom_hline(yintercept = -log(pval_thresh), col = "red") +
  geom_vline(xintercept = epsilon_thresh, col = "red") + 
  #ylim(c(0, max(indiv_cont_ep_df$pval, -log(pval_thresh)))) + 
 # xlim(c(0, max(indiv_cont_ep_df$epsilon, epsilon_thresh))) + 
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold"), 
        legend.position = "bottom") + 
  ggtitle("")  + 
  ylab("-ln(P-value)") +
  coord_cartesian(xlim = c(0, .25), 
                  ylim = c(0, 2.5), 
                           clip = "off") + 
  #geom_repel(
    ggrepel::geom_text_repel(direction = "y", 
      xlim = c(0.2, 0.3), 
      ylim = c(1, 2.5), 
    data = subset(
    indiv_cont_ep_df, epsilon > epsilon_thresh & pval > -log(pval_thresh)), 
    aes(x = epsilon, y = pval, label = even_shorter_name), 
    max.overlaps = 15) 

#tree_htmp_pheno_plot

y_width <- 10
y_fontsize <- 9
y_num_size <- 6
y_angle <- 0 
y_title_angle <- 0

# Prep data ----
tree <- read.tree("../../data/3_trees/log_toxin.tree")
ribotype <- read.table("../../data/6_ribotype/log_toxin_ribotype.tsv", 
                       sep = "\t", 
                       stringsAsFactors = FALSE)
ribotype <- as.data.frame(cbind(row.names(ribotype), ribotype), 
                          stringsAsFactors = FALSE)
colnames(ribotype) <- c("ID", "Ribotype", "ribotype_color")
ribotype_colors <- unique(ribotype$ribotype_color)
names(ribotype_colors) <- unique(ribotype$Ribotype)

tree_plot_wo_leg <- ggtree(tree) %<+% ribotype + 
  geom_tippoint(aes(color = Ribotype), size = .9) + 
  scale_color_manual(values = ribotype_colors) + 
  coord_flip() + 
  theme(legend.position = "none")

tree_plot_w_leg <- ggtree(tree) %<+% ribotype + 
  geom_tippoint(aes(color = Ribotype), size = .9) + 
  scale_color_manual(values = ribotype_colors) + 
  coord_flip() + 
  theme(legend.position = c(0.65, 1.0),
        legend.direction = "horizontal")
ribo_legend <- ggpubr::get_legend(tree_plot_w_leg)

# Toxin Bar Plot ----
toxin <- read.table("../../data/4_phenotypes/cleaned_log_toxin.tsv",
                    sep = "\t",
                    row.names = 1,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
toxin <- as.matrix(toxin)
toxin <- as.data.frame(cbind(row.names(toxin), toxin))
colnames(toxin) <- c("ID", "toxin")
toxin$ID <- as.character(toxin$ID)
toxin$toxin <- as.numeric(as.character(toxin$toxin))

toxin_plot <- 
  ggplot(toxin, aes(ID, toxin)) + 
  geom_col() + 
  ggtitle("Cytotoxicity") + 
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(angle = y_angle, size = y_num_size), 
        axis.title.y = element_text(size = y_fontsize, angle = y_title_angle)) + 
  ylab(lapply(strwrap("Equivalent Toxin B Activity ln(ng/ml)",
                      width = y_width, 
                      simplify = FALSE),
              paste, collapse = "\n")) + 
  annotate("rect",  # 078
           xmin = 1.5,
           xmax = 5.5,
           ymin = -Inf,
           ymax = Inf,
           alpha = .2,
           fill = "darkorchid") + 
  annotate("rect",  # 017
           xmin = 5.5,
           xmax = 6.5,
           ymin = -Inf,
           ymax = Inf,
           alpha = .2,
           fill = "blue") +
  annotate("rect", # 027 
           xmin = 14.5,
           xmax = 40.5,
           ymin = -Inf,
           ymax = Inf,
           alpha = .2,
           fill = "deepskyblue1") +
  annotate("rect", # 013
           xmin = 40.5,
           xmax = 41.5,
           ymin = -Inf,
           ymax = Inf,
           alpha = .2,
           fill = "darkgreen") +
  annotate("rect", # 053
           xmin = 41.5,
           xmax = 45.5,
           ymin = -Inf,
           ymax = Inf,
           alpha = .2,
           fill = "gold1") +
  annotate("rect", # 003
           xmin = 58.5,
           xmax = 59.5,
           ymin = -Inf,
           ymax = Inf,
           alpha = .2,
           fill = "orange") +
  annotate("rect", # 001
           xmin = 61.5,
           xmax = 65.5,
           ymin = -Inf,
           ymax = Inf,
           alpha = .2,
           fill = "red") +
  annotate("rect", # 014 1
           xmin = 82.5,
           xmax = 83.5,
           ymin = -Inf,
           ymax = Inf,
           alpha = .2,
           fill = "brown") + 
  annotate("rect", # 014 2
           xmin = 84.5,
           xmax = 102.5,
           ymin = -Inf,
           ymax = Inf,
           alpha = .2,
           fill = "brown") + 
  annotate("rect", # 014 3
           xmin = 103.5,
           xmax = 107.5,
           ymin = -Inf,
           ymax = Inf,
           alpha = .2,
           fill = "brown") +
  labs(tag = "C") +
  theme(plot.tag = element_text(face = "bold")) + 
  ggtitle("")

# sig toxin gene heatmap
sig_genotype <- data.table::fread(file = paste0(data_dir, "1_combine_genotype_matrices/combinedmat/log_toxin_sig_indiv_pan_and_rereferenced_snp_indel.tsv"),
                                  sep = "\t", 
                                  data.table = FALSE, 
                                  stringsAsFactors = FALSE, 
                                  check.names = FALSE)
sig_genotype <- as.matrix(sig_genotype)
row.names(sig_genotype) <- sig_genotype[, 1]
sig_genotype <- sig_genotype[, 2:ncol(sig_genotype)]
storage.mode(sig_genotype) <- "numeric"

heatmap.colours <- c("white","black")
names(heatmap.colours) <- 0:1

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

verbose_names <- sig_genotype_df_long$locus
sig_genotype_df_long$pretty_genotype <- verbose_names %>% 
  gsub("Coding ", "", .) %>% 
  gsub(" at .*locus_tag=.*Strand.*c[.]", " ", .) %>% 
  gsub("[|][0-9].*", "", .) %>% 
  gsub("[|]p[.]", " ", .) %>% 
  gsub("[|].*", "", .) %>% 
  gsub("fs", "frameshift", .)

sig_genotype_df_long$pretty_genotype

gene_name <- verbose_names %>% gsub(".*[0-9][|]", "", .) %>% gsub(";", "", .) %>% gsub("[|].*", "", .)

for (i in 1:nrow(sig_genotype_df_long)) {
  if (nchar(sig_genotype_df_long$pretty_genotype[i]) > 10) {
    sig_genotype_df_long$pretty_genotype[i] <- paste(gene_name[i], sig_genotype_df_long$pretty_genotype[i])
  }
}

sig_genotype_df_long  <- 
  sig_genotype_df_long %>% 
  mutate(even_shorter_name = gsub("Indel.*", "Indel", pretty_genotype))


heatmap_plot_p3 <- 
  ggplot(sig_genotype_df_long, aes(y = even_shorter_name, x = ID)) + 
  geom_tile(aes(fill = factor(presence))) + 
  scale_fill_manual(values = c("white", "black")) + 
  guides(fill = "none") + 
  ylab("") + 
    theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.line = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        rect = element_blank()) + 
  xlab("")


tree_htmp_pheno_plot <-
    heatmap_plot_p3 %>% 
    insert_bottom(tree_plot_wo_leg, height = 0.75) %>%
    insert_top(toxin_plot, height = 0.5) %>% 
    insert_bottom(ggpubr::as_ggplot(ribo_legend))
  

# Made in J_tcdb_enrichment/06_plot_paloc_figure.R
# pheatmap_plot
load("../../figures/paloc_name_tbd_3.RData")
gg_pheatmap_plot <- as.ggplot(pheatmap_plot)
new_gg_pheatmap_plot <- gg_pheatmap_plot + 
  labs(tag = "B") + 
  theme(plot.tag = element_text(face = "bold"),
        plot.margin = margin(t = 0, r = 0, b = 1, l = 0, unit =  "cm")) 
  

layout <- c(
  area(t = 1, l = 1, b = 1, r = 1), 
  area(t = 1, l = 2, b = 1, r = 3), 
  area(t = 2, l = 1, b = 2, r = 3)
)


new_indiv_plot + 
  as.grob(tree_htmp_pheno_plot) + 
  new_gg_pheatmap_plot + 
  plot_layout(design = layout)

ggsave(filename = "../../figures/Fig4.png",
       units = "in",
       height = 12,
       width = 14)
ggsave(filename = "../../figures/Fig4.pdf",
       units = "in",
       height = 12,
       width = 14) 
