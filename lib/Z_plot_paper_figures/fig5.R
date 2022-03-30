library(aplot)
library(patchwork)
library(ggplotify)
library(ggplot2)
library(ape)
library(ggtree)
library(forcats)


y_width <- 15
y_fontsize <- 9
y_num_size <- 6
y_angle <- 0 
y_title_angle <- 0

# ../H_pyseer_all/compare_wgs_paloc_and_tcdb_pyseer_results.R
# pyseer_plot
load("../../figures/pyseer_results/all_three_EN_models_jitter_best_alpha_log_toxin.RData")
new_pyseer_plot <- pyseer_plot + 
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold")) + 
  theme(legend.position = "none")


# ribotype
ribotype <- read.table("../../data/6_ribotype/log_toxin_ribotype.tsv", 
                       sep = "\t", 
                       stringsAsFactors = FALSE)
ribotype <- as.data.frame(cbind(row.names(ribotype), ribotype), 
                          stringsAsFactors = FALSE)
colnames(ribotype) <- c("ID", "Ribotype", "ribotype_color")
ribotype_colors <- unique(ribotype$ribotype_color)
names(ribotype_colors) <- unique(ribotype$Ribotype)


# toxin
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

toxin_plot_B <- 
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
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold"))
  

toxin_plot_C <- 
  toxin_plot_B + 
  labs(tag = "C") +
  theme(plot.tag = element_text(face = "bold"))


# B - WG tree plot
# Tree
wg_tree <- read.tree("../../data/3_trees/log_toxin.tree")
wg_tree_plot_w_leg <- ggtree(wg_tree) %<+% ribotype + 
  geom_tippoint(aes(color = Ribotype), size = 0.5) + 
  scale_color_manual(values = ribotype_colors) + 
  coord_flip() + 
  theme(legend.title = element_text(size = 10), # legend.position = c(0.65, 0.0),   
    legend.position = c(0.40, 0.0), 
        legend.text = element_text(size = 8), 
        legend.direction = "horizontal") + 
  ggtitle("WG")

wg_tree_plot_wo_leg <- ggtree(wg_tree) %<+% ribotype + 
  geom_tippoint(aes(color = Ribotype), size = 0.5) + 
  scale_color_manual(values = ribotype_colors) + 
  coord_flip() + 
  theme(legend.position = "none") + 
  ggtitle("WG")

# # legend
ribo_legend <- ggpubr::get_legend(wg_tree_plot_w_leg)

wg_tree_w_tox_plot <- 
  toxin_plot_B %>% 
  insert_bottom(wg_tree_plot_wo_leg, height = 4) %>% 
  insert_bottom(ggpubr::as_ggplot(ribo_legend))

# C - tcdB tree plot 
# Tree
gene_tree <- read.tree("../../data/3_trees/gene_trees/cleaned_log_toxin_midpoint_root_tcdB_gene_prank.tree")
gene_tree_plot_wo_leg <- ggtree(gene_tree) %<+% ribotype + 
  geom_tippoint(aes(color = Ribotype), size = 0.5) + 
  scale_color_manual(values = ribotype_colors) + 
  coord_flip() + 
  theme(legend.position = "none") + 
  ggtitle(expression(italic("tcdB")))

# gene_tree_w_tox_plot <- 
#   toxin_plot_C %>% 
#   insert_bottom(gene_tree, height = 4)

gene_tree_w_tox_plot <- 
  toxin_plot_C %>% 
  insert_bottom(gene_tree_plot_wo_leg, height = 4) %>% 
  insert_bottom(ggpubr::as_ggplot(ribo_legend))

layout <- c(
  area(1, 1), 
  area(1, 2, 1, 3), 
  area(1, 4, 1, 5)
)

# Check out layout design
# plot(layout)

new_pyseer_plot + as.ggplot(wg_tree_w_tox_plot) + 
  as.ggplot(gene_tree_w_tox_plot) + plot_layout(design = layout) + 
  plot_layout(guides = 'keep')

ggsave("../../figures/Fig5.png", width = 14, height = 5, units = "in")
ggsave("../../figures/Fig5.pdf", width = 14, height = 5, units = "in")


