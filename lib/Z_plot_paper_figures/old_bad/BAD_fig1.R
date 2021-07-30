# this is all bad. use fig1.R instead

# # TODO:  phenotypes ARE  plotting wrong -- toxin is too low for 027 and severity is wrong. VPI (003) isn't correct - it is ND except for toxin and growth and that's not what the plot looks like.
# # TODO: add ND symbols once the plotting is correct.
# # TODO: double check -- but I think the non-functionalized version plots correct. WHAT?
# # For a check - there should be 10 severe 027 and 7 sever 014, 1 sever 053, and 1 sever 078
# source("tree_ribo_heatmap_pheno_helper_functions.R")
# library(ggpubr)
# 
# # example used with agr variants and the PSM dataset
# # tree <- read_and_midpoint_tree("../../data/3_trees/log_toxin.tree")
# # ribotype <- read_format_ribotype("../../data/6_ribotype/log_toxin_ribotype.tsv")
# # ribotype_colors <- create_ribotype_colors_obj(ribotype)
# # sev_df <- read_format_pheno("../../data/4_phenotypes/cleaned_severity.tsv", "severity")
# # tox_df <- read_format_pheno("../../data/4_phenotypes/cleaned_log_toxin.tsv", "log_toxin")
# # 
# # # genotype <- read_prep_geno_for_heatmap("../../data/1_combine_genotype_matrices/combinedmat/cleaned_log_toxin_pan_and_rereferenced_snp_indel.tsv")
# # 
# # 
# # tree_plot <- vert_tree_w_ribo(tree, ribotype, ribotype_colors)
# # sev_bar_plot <- vert_pheno_bar_plot_no_ribo(sev_df, "severity", plot_title="")
# # tox_bar_plot <- vert_pheno_bar_plot_no_ribo(tox_df, "log_toxin", plot_title="")
# # # heatmap_plot <- vert_heatmap(genotype)
# # # 
# # # 
# # # heatmap_plot %>% insert_bottom(tree_plot) %>% insert_top(bar_plot) 
# # tree_plot %>% insert_top(sev_bar_plot) %>% insert_top(tox_bar_plot) 
# 
# 
# phenotype <- read_tsv("../../data/4_phenotypes/all_phenos.tsv")
# tree <- read.tree("../../data/3_trees/log_toxin.tree")
# ribotype <- read.table("../../data/6_ribotype/log_toxin_ribotype.tsv", sep = "\t", stringsAsFactors = FALSE)       
# ribotype <- as.data.frame(cbind(row.names(ribotype), ribotype), stringsAsFactors = FALSE)
# colnames(ribotype) <- c("ID", "Ribotype", "ribotype_color")
# ribotype_colors <- unique(ribotype$ribotype_color)
# names(ribotype_colors) <- unique(ribotype$Ribotype)
# 
# fs <- 5
# fs2 <- 5
# nd_text_y <- 0.15
# nd_text_size <- 1.5
# 
# tree_plot_wo_leg <- ggtree(tree) %<+% ribotype + 
#   geom_tippoint(aes(color = Ribotype)) + 
#   scale_color_manual(values = ribotype_colors) + 
#   coord_flip() + 
#   theme(legend.position = "none")
# 
# tree_plot_w_leg <- ggtree(tree) %<+% ribotype + 
#   geom_tippoint(aes(color = Ribotype)) + 
#   scale_color_manual(values = ribotype_colors) + 
#   coord_flip() + 
#   theme(legend.position = c(0.55, 1.75),
#         legend.direction = "horizontal")
# 
# ribo_legend <- ggpubr::get_legend(tree_plot_w_leg)
# 
# make_bar_plot <- function(phenotype, pheno_name, yaxis_name){
#   bar_plot <- ggplot(phenotype) + 
#     geom_col(aes_string("ID", pheno_name)) + 
#     theme_classic() + 
#     ylab(yaxis_name) + 
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(), 
#           axis.title.y = element_text(size = fs), 
#           axis.text.y = element_text(size = fs2)) + 
#     annotate("rect",  # 078
#              xmin = 1.5,
#              xmax = 5.5,
#              ymin = -Inf,
#              ymax = Inf,
#              alpha = .2,
#              fill = "darkorchid") + 
#     annotate("rect",  # 017
#              xmin = 5.5,
#              xmax = 6.5,
#              ymin = -Inf,
#              ymax = Inf,
#              alpha = .2,
#              fill = "blue") +
#     annotate("rect", # 027 
#              xmin = 14.5,
#              xmax = 40.5,
#              ymin = -Inf,
#              ymax = Inf,
#              alpha = .2,
#              fill = "deepskyblue1") +
#     annotate("rect", # 013
#              xmin = 40.5,
#              xmax = 41.5,
#              ymin = -Inf,
#              ymax = Inf,
#              alpha = .2,
#              fill = "darkgreen") +
#     annotate("rect", # 053
#              xmin = 41.5,
#              xmax = 45.5,
#              ymin = -Inf,
#              ymax = Inf,
#              alpha = .2,
#              fill = "gold1") +
#     annotate("rect", # 003
#              xmin = 58.5,
#              xmax = 59.5,
#              ymin = -Inf,
#              ymax = Inf,
#              alpha = .2,
#              fill = "orange") +
#     annotate("rect", # 001
#              xmin = 61.5,
#              xmax = 65.5,
#              ymin = -Inf,
#              ymax = Inf,
#              alpha = .2,
#              fill = "red") +
#     annotate("rect", # 014 1
#              xmin = 82.5,
#              xmax = 83.5,
#              ymin = -Inf,
#              ymax = Inf,
#              alpha = .2,
#              fill = "brown") + 
#     annotate("rect", # 014 2
#              xmin = 84.5,
#              xmax = 102.5,
#              ymin = -Inf,
#              ymax = Inf,
#              alpha = .2,
#              fill = "brown") + 
#     annotate("rect", # 014 3
#              xmin = 103.5,
#              xmax = 107.5,
#              ymin = -Inf,
#              ymax = Inf,
#              alpha = .2,
#              fill = "brown")
#   
#   return(bar_plot)
# }
# 
# # add in ND points: 
# #   geom_text(data = NULL, x = 42, y = nd_text_y, label = "ND", size = nd_text_size) + 
# spor_nd <- phenotype %>% filter(is.na(log_sporulation)) %>% pull(ID) %>% as.character() # "Cd033" "Cd086" "VPI"  
# cfe_nd <- phenotype %>% filter(is.na(log_cfe)) %>% pull(ID) %>% as.character()# "DA00484" "VPI"     "DA00315"
# growth_nd <- phenotype %>% filter(is.na(log_growth)) %>% pull(ID) %>% as.character() # "Cd054" "Cd182" "Cd072" "Cd093"
# 
# germ_tc_and_gly_nd <- phenotype %>% filter(is.na(log_germ_tc_and_gly)) %>% pull(ID) %>% as.character() #  "Cd054" "VPI"   "Cd182" "Cd072" "Cd093"  
# germ_tc_only_nd <- phenotype %>% filter(is.na(log_germ_tc)) %>% pull(ID) %>% as.character() # "Cd054"   "VPI"     "Cd182"   "Cd072"   "Cd093"   "DA00330"
# 
# 
# 
# severity_bar_plot <- make_bar_plot(phenotype, "severity", "severe infection")
# toxin_bar_plot <- make_bar_plot(phenotype, "log_toxin", "toxin activity")
# germ_TG_bar_plot <- make_bar_plot(phenotype, "log_germ_tc_and_gly", "germination in Tc and Gly)")
# germ_bar_plot <- make_bar_plot(phenotype, "log_germ_tc", "germination in Tc")
# cfe_bar_plot <- make_bar_plot(phenotype, "log_cfe", "spore viability")
# spor_bar_plot <- make_bar_plot(phenotype, "log_sporulation", "sporulation")
# grow_bar_plot <- make_bar_plot(phenotype, "log_growth", "growth rate")
# 
# png("../../figures/Fig1.png"),
#      units = "in",
#      width = 9, height = 10, res = 250)
# print(cfe_bar_plot %>% 
#         insert_bottom(tree_plot_wo_leg, height = 4) %>% 
#         insert_top(germ_bar_plot) %>% 
#         insert_top(germ_TG_bar_plot) %>% 
#         insert_top(grow_bar_plot) %>% 
#         insert_top(spor_bar_plot) %>% 
#         insert_top(toxin_bar_plot) %>%
#         insert_top(severity_bar_plot) %>% 
#         insert_bottom(ggpubr::as_ggplot(ribo_legend)))
# dev.off()
# 
# print("TODO: add ND to plot")
# stop()
# 
# # Weirdly, I get an error if I try to put the plots left to right, but not if I plot them top to bottom. Oh well. 
# 
# 
# ## Let's do it the bad way
# bar_plot_toxin <- ggplot(phenotype) + 
#   geom_col(aes(ID, as.numeric(as.character(log_toxin)))) + 
#   # ggtitle("Toxin") + 
#   theme_classic() + 
#   ylab("ln(toxin activity)") + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), 
#         axis.title.y = element_text(size = fs), 
#         axis.text.y = element_text(size = fs2)) + 
#   annotate("rect",  # 017 1
#            xmin = 4.5,
#            xmax = 5.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkorchid") + 
#   annotate("rect",  # 027 1
#            xmin = 9.5,
#            xmax = 33.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 027 2
#            xmin = 34.5,
#            xmax = 36.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 053 1
#            xmin = 53.5,
#            xmax = 54.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 013 1
#            xmin = 54.5,
#            xmax = 55.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "orange") +
#   annotate("rect", # 053 2
#            xmin = 55.5,
#            xmax = 58.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 014 1
#            xmin = 77.5,
#            xmax = 86.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 001 1
#            xmin = 86.5,
#            xmax = 87.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 003 1
#            xmin = 87.5,
#            xmax = 88.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "red") +
#   annotate("rect", # 001 1
#            xmin = 88.5,
#            xmax = 91.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 014 2
#            xmin = 91.5,
#            xmax = 97.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 3
#            xmin = 98.5,
#            xmax = 104.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 4
#            xmin = 105.5,
#            xmax = 107.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 078 1
#            xmin = 40.5,
#            xmax = 41.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen") + 
#   annotate("rect", # 078 2
#            xmin = 42.5,
#            xmax = 45.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen")
# 
# 
# # geom_rect(fill = ribotype$ribotype_color, 
# #           xmin = tox_x_min, 
# #           xmax = tox_x_max,
# #           ymin = min(as.numeric(as.character(phenotype$log_toxin))), 
# #           ymax = max(as.numeric(as.character(phenotype$log_toxin))), 
# #           alpha = 0.2)
# 
# bar_plot_cfe <- ggplot(phenotype, aes(ID,  as.numeric(as.character(log_cfe)))) +
#   geom_col(na.rm = FALSE) +
#   #ggtitle("CFE") + 
#   theme_classic() +  
#   ylab("ln(spore viability)") + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), 
#         axis.title.y = element_text(size = fs), 
#         axis.text.y = element_text(size = fs2)) + 
#   # which(tree$tip.label %in% cfe_nd)
#   geom_text(data = NULL, x = 42, y = nd_text_y, label = "ND", size = nd_text_size) + 
#   geom_text(data = NULL, x = 88, y = nd_text_y, label = "ND", size = nd_text_size) + # VPI / RT003
#   geom_text(data = NULL, x = 67, y = nd_text_y, label = "ND", size = nd_text_size) + 
#   annotate("rect",  # 017 1
#            xmin = 4.5,
#            xmax = 5.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkorchid") + 
#   annotate("rect",  # 027 1
#            xmin = 9.5,
#            xmax = 33.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 027 2
#            xmin = 34.5,
#            xmax = 36.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 053 1
#            xmin = 53.5,
#            xmax = 54.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 013 1
#            xmin = 54.5,
#            xmax = 55.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "orange") +
#   annotate("rect", # 053 2
#            xmin = 55.5,
#            xmax = 58.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 014 1
#            xmin = 77.5,
#            xmax = 86.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 001 1
#            xmin = 86.5,
#            xmax = 87.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 003 1
#            xmin = 87.5,
#            xmax = 88.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "red") +
#   annotate("rect", # 001 1
#            xmin = 88.5,
#            xmax = 91.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 014 2
#            xmin = 91.5,
#            xmax = 97.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 3
#            xmin = 98.5,
#            xmax = 104.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 4
#            xmin = 105.5,
#            xmax = 107.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 078 1
#            xmin = 40.5,
#            xmax = 41.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen") + 
#   annotate("rect", # 078 2
#            xmin = 42.5,
#            xmax = 45.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen")
# 
# bar_plot_growth <- ggplot(phenotype, aes(ID,  as.numeric(as.character(log_growth)))) +
#   geom_col(na.rm = FALSE) +
#   #ggtitle("Growth") +  
#   ylab("ln(growth rate)") + 
#   theme_classic() + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), 
#         axis.title.y = element_text(size = fs), 
#         axis.text.y = element_text(size = fs2)) + 
#   geom_text(data = NULL, x = 14, y = -nd_text_y, label = "ND", size = nd_text_size) + # Cd093
#   geom_text(data = NULL, x = 43, y = -nd_text_y, label = "ND", size = nd_text_size) + # Cd054
#   geom_text(data = NULL, x = 83, y = -nd_text_y, label = "ND", size = nd_text_size) + # Cd072
#   geom_text(data = NULL, x = 103, y = -nd_text_y, label = "ND", size = nd_text_size) + # Cd182
#   annotate("rect",  # 017 1
#            xmin = 4.5,
#            xmax = 5.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkorchid") + 
#   annotate("rect",  # 027 1
#            xmin = 9.5,
#            xmax = 33.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 027 2
#            xmin = 34.5,
#            xmax = 36.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 053 1
#            xmin = 53.5,
#            xmax = 54.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 013 1
#            xmin = 54.5,
#            xmax = 55.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "orange") +
#   annotate("rect", # 053 2
#            xmin = 55.5,
#            xmax = 58.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 014 1
#            xmin = 77.5,
#            xmax = 86.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 001 1
#            xmin = 86.5,
#            xmax = 87.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 003 1
#            xmin = 87.5,
#            xmax = 88.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "red") +
#   annotate("rect", # 001 1
#            xmin = 88.5,
#            xmax = 91.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 014 2
#            xmin = 91.5,
#            xmax = 97.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 3
#            xmin = 98.5,
#            xmax = 104.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 4
#            xmin = 105.5,
#            xmax = 107.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 078 1
#            xmin = 40.5,
#            xmax = 41.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen") + 
#   annotate("rect", # 078 2
#            xmin = 42.5,
#            xmax = 45.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen")
# 
# 
# bar_plot_sporulation <- ggplot(phenotype, aes(ID,  as.numeric(as.character(log_sporulation)))) +
#   geom_col(na.rm = FALSE) +
#   #ggtitle("Sporulation") +  
#   ylab("ln(# spores)") + 
#   theme_classic() + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), 
#         axis.title.y = element_text(size = fs), 
#         axis.text.y = element_text(size = fs2))+ 
#   geom_text(data = NULL, x = 88, y = nd_text_y, label = "ND", size = nd_text_size) + # VPI / RT003
#   geom_text(data = NULL, x = 45, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd033
#   geom_text(data = NULL, x = 53, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd086
#   annotate("rect",  # 017 1
#            xmin = 4.5,
#            xmax = 5.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkorchid") + 
#   annotate("rect",  # 027 1
#            xmin = 9.5,
#            xmax = 33.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 027 2
#            xmin = 34.5,
#            xmax = 36.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 053 1
#            xmin = 53.5,
#            xmax = 54.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 013 1
#            xmin = 54.5,
#            xmax = 55.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "orange") +
#   annotate("rect", # 053 2
#            xmin = 55.5,
#            xmax = 58.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 014 1
#            xmin = 77.5,
#            xmax = 86.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 001 1
#            xmin = 86.5,
#            xmax = 87.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 003 1
#            xmin = 87.5,
#            xmax = 88.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "red") +
#   annotate("rect", # 001 1
#            xmin = 88.5,
#            xmax = 91.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 014 2
#            xmin = 91.5,
#            xmax = 97.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 3
#            xmin = 98.5,
#            xmax = 104.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 4
#            xmin = 105.5,
#            xmax = 107.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 078 1
#            xmin = 40.5,
#            xmax = 41.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen") + 
#   annotate("rect", # 078 2
#            xmin = 42.5,
#            xmax = 45.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen")
# 
# 
# bar_plot_germ_tc <- ggplot(phenotype, aes(ID,  as.numeric(as.character(log_germ_tc)))) +
#   geom_col(na.rm = FALSE) +
#   # ggtitle("Germ Tc") +  
#   ylab("ln(germination in Tc)") + 
#   theme_classic() + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), 
#         axis.title.y = element_text(size = fs), 
#         axis.text.y = element_text(size = fs2))+ 
#   geom_text(data = NULL, x = 14, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd093
#   geom_text(data = NULL, x = 43, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd054
#   geom_text(data = NULL, x = 83, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd072
#   geom_text(data = NULL, x = 103, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd182
#   geom_text(data = NULL, x = 88, y = nd_text_y, label = "ND", size = nd_text_size) + # VPI / RT003
#   geom_text(data = NULL, x = 5, y = nd_text_y, label = "ND", size = nd_text_size) + # DA00330 / RT 017
#   annotate("rect",  # 017 1
#            xmin = 4.5,
#            xmax = 5.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkorchid") + 
#   annotate("rect",  # 027 1
#            xmin = 9.5,
#            xmax = 33.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 027 2
#            xmin = 34.5,
#            xmax = 36.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 053 1
#            xmin = 53.5,
#            xmax = 54.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 013 1
#            xmin = 54.5,
#            xmax = 55.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "orange") +
#   annotate("rect", # 053 2
#            xmin = 55.5,
#            xmax = 58.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 014 1
#            xmin = 77.5,
#            xmax = 86.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 001 1
#            xmin = 86.5,
#            xmax = 87.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 003 1
#            xmin = 87.5,
#            xmax = 88.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "red") +
#   annotate("rect", # 001 1
#            xmin = 88.5,
#            xmax = 91.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 014 2
#            xmin = 91.5,
#            xmax = 97.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 3
#            xmin = 98.5,
#            xmax = 104.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 4
#            xmin = 105.5,
#            xmax = 107.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 078 1
#            xmin = 40.5,
#            xmax = 41.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen") + 
#   annotate("rect", # 078 2
#            xmin = 42.5,
#            xmax = 45.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen")
# 
# 
# bar_plot_germ_tc_and_gly <- ggplot(phenotype, aes(ID,  as.numeric(as.character(log_germ_tc_and_gly)))) +
#   geom_col(na.rm = FALSE) +
#   #ggtitle("Germ Tc & Gly") +  
#   ylab("ln(germination in Tc and Gly)") + 
#   theme_classic() + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), 
#         axis.title.y = element_text(size = fs), 
#         axis.text.y = element_text(size = fs2))+ 
#   geom_text(data = NULL, x = 14, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd093
#   geom_text(data = NULL, x = 43, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd054
#   geom_text(data = NULL, x = 83, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd072
#   geom_text(data = NULL, x = 103, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd182
#   geom_text(data = NULL, x = 88, y = nd_text_y, label = "ND", size = nd_text_size) + # VPI / RT003
#   annotate("rect",  # 017 1
#            xmin = 4.5,
#            xmax = 5.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkorchid") + 
#   annotate("rect",  # 027 1
#            xmin = 9.5,
#            xmax = 33.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 027 2
#            xmin = 34.5,
#            xmax = 36.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "blue") +
#   annotate("rect", # 053 1
#            xmin = 53.5,
#            xmax = 54.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 013 1
#            xmin = 54.5,
#            xmax = 55.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "orange") +
#   annotate("rect", # 053 2
#            xmin = 55.5,
#            xmax = 58.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "aquamarine") +
#   annotate("rect", # 014 1
#            xmin = 77.5,
#            xmax = 86.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 001 1
#            xmin = 86.5,
#            xmax = 87.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 003 1
#            xmin = 87.5,
#            xmax = 88.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "red") +
#   annotate("rect", # 001 1
#            xmin = 88.5,
#            xmax = 91.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deepskyblue1") +
#   annotate("rect", # 014 2
#            xmin = 91.5,
#            xmax = 97.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 3
#            xmin = 98.5,
#            xmax = 104.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 014 4
#            xmin = 105.5,
#            xmax = 107.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "deeppink") +
#   annotate("rect", # 078 1
#            xmin = 40.5,
#            xmax = 41.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen") + 
#   annotate("rect", # 078 2
#            xmin = 42.5,
#            xmax = 45.5,
#            ymin = -Inf,
#            ymax = Inf,
#            alpha = .2,
#            fill = "darkgreen")
# bar_plot_sev <- ggplot(phenotype, aes(ID,  as.numeric(as.character(severity)))) +
#   geom_col(na.rm = FALSE) +
#   #ggtitle("Germ Tc & Gly") +  
#   ylab("ln(germination in Tc and Gly)") + 
#   theme_classic() + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), 
#         axis.title.y = element_text(size = fs), 
#         axis.text.y = element_text(size = fs2))  
# print(bar_plot_cfe %>% 
#         insert_bottom(tree_plot_wo_leg, height = 4) %>% 
#         insert_top(bar_plot_germ_tc) %>% 
#         insert_top(bar_plot_germ_tc_and_gly) %>% 
#         insert_top(bar_plot_growth) %>% 
#         insert_top(bar_plot_sporulation) %>% 
#         insert_top(bar_plot_toxin) %>% 
#         insert_top(bar_plot_sev) %>% 
#         insert_bottom(ggpubr::as_ggplot(ribo_legend)))
