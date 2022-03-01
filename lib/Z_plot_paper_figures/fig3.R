# Figure 5
# Overlapping GWAS results.
# (A) Heatmap indicates the number of shared individual locus GWAS results with
# significant P-values (as defined in threshold file) and high levels of
# convergence (as defined in threshold file) in Continuous Test (continuous
# phenotypes) or Synchronous Test (Severity).
# (B) Shared hits: toxin activity and severe infection GWASs. 

options(scipen = 999) # turn off scientific notation
#source("../summarize_all_phenotype_data/all_phenotype_analysis_lib.R") 
# source("../analysis_and_plots_for_thesis_chapter3/color_palettes.R")
library(tidyverse)
library(patchwork)
library(dplyr)
library(phytools)
library(ape)
library(ggtree)
library(aplot)
library(ggpubr)
library(ggplotify)
phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")


# A ----
# Data and original plot from summarize_all_phenotype_data/get_colocalization_of_howgwash_results.R
intersect_mat <- read.csv("../../data/13_summaries/number_overlapping_gwas_hits.csv", row.names = 1)
intersect_mat <- as.matrix(intersect_mat)


# make the pretty names
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

# Remove fqr because not relevant to paper, remove other two because they don't have any significant hits at all
bad_phenotypes <- c("ln(germination in Tc and Gly)", "ln(spore viability)", "fqR")
intersect_mat <- intersect_mat[!(row.names(intersect_mat) %in% bad_phenotypes), 
                     !(colnames(intersect_mat) %in% bad_phenotypes)]



colnames(intersect_mat) <- gsub("ln[(]", "", colnames(intersect_mat)) %>% gsub("[)]", "", .) %>% gsub(" and ", "&", .) %>% gsub("Severity", "severe infection", .)
row.names(intersect_mat) <- gsub("ln[(]", "", row.names(intersect_mat)) %>% gsub("[)]", "", .) %>% gsub(" and ", "&", .) %>% gsub("Severity", "severe infection", .)

# remove self-self diagonal
intersect_mat[1, 1] <- NA
intersect_mat[2, 2] <- NA
intersect_mat[3, 3] <- NA
intersect_mat[4, 4] <- NA
intersect_mat[5, 5] <- NA


# reorder martix to match other plots
# Want to change plot order to: 
# 1. germ tc
# 2. germ TG
# 3. growth
# 4. # sportes
# 5. spore viability
# 6. toxin 
name_order <- c("germination in Tc", #"germination in Tc&Gly", 
                "growth rate", "# spores", 
               # "spore viability", 
                "toxin activity", "severe infection")
intersect_mat <- intersect_mat[match(name_order, row.names(intersect_mat)),, drop = FALSE]
intersect_mat <- intersect_mat[, match(name_order, colnames(intersect_mat)), drop = FALSE]

# Remove redundant half
intersect_mat[lower.tri(intersect_mat)] <- NA

# Add asterisks for significantly enriched values
perm_pval_mat <- read.csv(row.names = 1, file = "../../data/13_summaries/overlapping_permutation_bonf_corrected_p_values.csv")
perm_pval_mat <- perm_pval_mat[c(1, 3, 4), 1:2] # keep just rows with significant results

# plot
melted_intersect_mat <- reshape2::melt(intersect_mat, na.rm = TRUE)
melted_intersect_mat <- cbind(melted_intersect_mat, rep(NA, nrow(melted_intersect_mat)))
colnames(melted_intersect_mat)[4] <- "Significant"

melted_intersect_mat$Significant[
  melted_intersect_mat$Var1 == "germination in Tc" &
    melted_intersect_mat$Var2 == "severe infection"] <- "*"

melted_intersect_mat$Significant[
  melted_intersect_mat$Var1 == "# spores" &
    melted_intersect_mat$Var2 == "severe infection"] <- "*"

melted_intersect_mat$Significant[
  melted_intersect_mat$Var1 == "toxin activity" &
    melted_intersect_mat$Var2 == "severe infection"] <- "*"

melted_intersect_mat$Significant[
  melted_intersect_mat$Var1 == "germination in Tc" &
    melted_intersect_mat$Var2 == "toxin activity"] <- "*"

melted_intersect_mat <- 
  melted_intersect_mat %>% 
  mutate(Var1 = case_when(Var1 == "germination in Tc" ~ "Germ. in TCA", 
                          Var1 == "germination in Tc&Gly" ~ "Germ. in TCA&Gly", 
                          Var1 == "fqR" ~ "Fluoroquinolone Resistance", 
                          Var1 == "spore viability" ~ "Spore viability", 
                          Var1 == "growth rate" ~ "Max growth rate", 
                          Var1 == "# spores" ~ "Total spores",
                          Var1 == "toxin activity" ~ "Cytotoxicity", 
                          Var1 == "severe infection" ~ "Severe infection")) %>% 
  
  mutate(Var2 = case_when(Var2 == "germination in Tc" ~ "Germ. in TCA", 
                          Var2 == "germination in Tc&Gly" ~ "Germ. in TCA&Gly", 
                          Var2 == "fqR" ~ "Fluoroquinolone Resistance", 
                          Var2 == "spore viability" ~ "Spore viability", 
                          Var2 == "growth rate" ~ "Max growth rate", 
                          Var2 == "# spores" ~ "Total spores",
                          Var2 == "toxin activity" ~ "Cytotoxicity", 
                          Var2 == "severe infection" ~ "Severe infection")) %>% 
  mutate(plot_order_x = 5, 
         plot_order_y = 5) %>% 
  mutate(plot_order_x = case_when(
                          Var2 ==  "Max growth rate" ~ 1, 
                          Var2 == "Total spores" ~ 2,
                          Var2 == "Cytotoxicity" ~ 3, 
                          Var2 == "Severe infection" ~ 4)) %>% 
  mutate(plot_order_y = case_when(
    Var1 == "Total spores" ~ 2,
    Var1 == "Cytotoxicity" ~ 3, 
    Var1 == "Germ. in TCA" ~ 1))

htmp_plot <-
  ggplot(data = melted_intersect_mat,
         aes(x = fct_reorder(Var2, plot_order_x), 
             y = fct_reorder(Var1, plot_order_y),
             fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    high = "red", mid = "grey",
    name = "Overlapping\nresults"
  ) +
  theme_minimal() +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_fixed() +
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold")) + 
  geom_text(aes(label = Significant), size = 5)#, nudge_y = -0.25)

# B ----

y_width <- 15
y_fontsize <- 9
y_num_size <- 6
y_angle <- 0 
y_title_angle <- 0

# Prep data ----
tree <- read.tree("../../data/3_trees/log_toxin.tree")
ribotype <- read.table("../../data/6_ribotype/log_toxin_ribotype.tsv", 
                       sep = "\t", 
                       stringsAsFactors = FALSE)

phenotype <- read.table("../../data/4_phenotypes/cleaned_log_toxin.tsv",
                        sep = "\t",
                        row.names = 1,
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        check.names = FALSE)
phenotype <- as.matrix(phenotype)

# sig_genotype
load("../../data/1_combine_genotype_matrices/combinedmat/cleaned_toxin_sig_indiv_pan_and_rereferenced_snp_indel.RData")
dim(sig_genotype) # 107 220

# intersect_df_list
load("../../data/13_summaries/overlapping_gwas_hits.rda")
tox_sev_df <- intersect_df_list$`log_toxin-severity`  %>% 
  filter(pheno == "log_toxin") # filter so genos not duplicated for ploting 
sig_genotype <- 
  sig_genotype[, colnames(sig_genotype) %in% tox_sev_df$Genotype, drop = FALSE]
ribotype <- as.data.frame(cbind(row.names(ribotype), ribotype), 
                          stringsAsFactors = FALSE)
colnames(ribotype) <- c("ID", "Ribotype", "ribotype_color")

phenotype <- as.data.frame(cbind(row.names(phenotype), phenotype))
colnames(phenotype) <- c("ID", "toxin")
phenotype$ID <- as.character(phenotype$ID)
phenotype$toxin <- as.numeric(as.character(phenotype$toxin))

ribotype_colors <- unique(ribotype$ribotype_color)
names(ribotype_colors) <- unique(ribotype$Ribotype)

# severity on / off  ---- 
sev_df <- read.table("../../data/4_phenotypes/cleaned_severity.tsv")
sev_df <- as.data.frame(sev_df)
sev_df <- cbind(row.names(sev_df), sev_df)
colnames(sev_df) <- c("ID", "Severe infection")

sev_plot <-
  ggplot(sev_df, aes(ID, `Severe infection`)) +
  geom_col() +
  ggtitle("Severe infection") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = y_fontsize, angle = y_title_angle)
  ) +
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold")) + 
  ylab("Severe infection") + 
  geom_text(data = NULL, x = 5,  y = nd_text_y, label = "ND", size = nd_text_size) + # Cd054 078-126
  geom_text(data = NULL, x = 25, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd093 027
  geom_text(data = NULL, x = 59, y = nd_text_y, label = "ND", size = nd_text_size) + # VPI 003
  geom_text(data = NULL, x = 86, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd182 014-020
  geom_text(data = NULL, x = 97, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd072 014-020
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
           fill = "brown")


# Tree plot ----
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
bar_plot_p1 <- 
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
           fill = "brown")


# Presence /absence matrix ----
heatmap.colours <- c("white", "black")
names(heatmap.colours) <- 0:1

if (ncol(sig_genotype) < 2) {
  sig_genotype <- cbind(sig_genotype, sig_genotype)
}
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
sig_genotype_df_long$ID <- factor(sig_genotype_df_long$ID, levels = rownames(sig_genotype), labels = rownames(sig_genotype))
sig_genotype_df_long$locus <- factor(sig_genotype_df_long$locus, levels = colnames(sig_genotype)[fit$order], labels = colnames(sig_genotype)[fit$order])

# Clean up locus names
verbose_names <- sig_genotype_df_long$locus
sig_genotype_df_long$locus <- verbose_names %>% 
  gsub("Coding ", "", .) %>% 
  gsub("at .*locus_tag=", "",. ) %>% 
  gsub("Strand.*c[.]", " ", .) %>% 
  gsub("[|][0-9].*", "", .) %>% 
  gsub("[|]p[.]", " ", .) %>% 
  gsub("[|].*", "", .) %>% 
  gsub("fs", "frameshift", .)

gene_name <- verbose_names %>% gsub(".*[0-9][|]", "", .) %>% gsub(";", "", .) %>% gsub("[|].*", "", .)
sig_genotype_df_long$locus <- paste(gene_name, sig_genotype_df_long$locus)

# get rid of duplicate groups
sig_genotype_df_long$locus <- gsub(" group_.*", "", sig_genotype_df_long$locus)

# remove locus tag for indel b/c we have gene name
sig_genotype_df_long$locus <- gsub("[ ]CD630_[0-9]+[ ]", "", sig_genotype_df_long$locus )

heatmap_plot_p3 <- ggplot(sig_genotype_df_long, aes(y = locus, x = ID)) + 
  geom_tile(aes(fill = factor(presence))) + 
  scale_fill_manual(values = c("white", "black")) + 
  scale_y_discrete(position = "right") + 
  guides(fill = "none") + 
  theme_bw() + 
  theme(axis.line = element_blank(),
        axis.text.y.right = element_text(size = 5, color = "black"),
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

# Combine tox / sev results 
combined_plot <- heatmap_plot_p3 %>% 
  insert_bottom(tree_plot_wo_leg, height = 0.75) %>%
  insert_top(bar_plot_p1, height = 0.5) %>% 
  insert_top(sev_plot, height = 0.05) %>% 
  insert_bottom(ggpubr::as_ggplot(ribo_legend))

# Plot together is too hard b/c of 
# htmp_plot + tree_plot + 
#    plot_layout(ncol = 2, widths = c(1, 2)) & 
#    theme(plot.tag = element_text(face = "bold"))

# Plot all together ---- 
combined_plot_grob <- ggplotify::as.grob(combined_plot)
htmp_plot + combined_plot_grob  + 
  plot_layout(guides = 'keep') + 
  plot_layout(ncol = 2, widths = c(1, 6))
ggsave("../../figures/Fig3.png", width = 12, height = 6, units = "in")
ggsave("../../figures/Fig3.pdf", width = 12, height = 6, units = "in")

