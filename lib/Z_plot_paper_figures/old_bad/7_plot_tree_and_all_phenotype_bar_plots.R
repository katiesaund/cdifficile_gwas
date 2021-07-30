# # TODO: 
# # Use Fig1.R instead 
# options(scipen = 999) # turn off scientific notation
# library(tidyverse)
# library(dplyr)
# library(phytools)
# library(ape)
# library(ggtree)
# library(aplot)
# library(ggpubr)
# source("all_phenotype_analysis_lib.R")
# 
# data_dir <- "../../data/"
# fig_dir <- "../../figures/"
# 
# cfe <-  read.table(file = "../../data/4_phenotypes/cleaned_log_cfe.tsv",
#                    sep = "\t", 
#                    header = TRUE,
#                    row.names = 1, 
#                    stringsAsFactors = FALSE)
# 
# germ_tc_and_gly <-  read.table(file = "../../data/4_phenotypes/cleaned_log_germ_tc_and_gly.tsv",
#                                sep = "\t", 
#                                header = TRUE,
#                                row.names = 1, 
#                                stringsAsFactors = FALSE)
# 
# germ_tc <-  read.table(file = "../../data/4_phenotypes/cleaned_log_germ_tc.tsv",
#                        sep = "\t", 
#                        header = TRUE,
#                        row.names = 1, 
#                        stringsAsFactors = FALSE)
# 
# growth <-  read.table(file = "../../data/4_phenotypes/cleaned_log_growth.tsv",
#                       sep = "\t", 
#                       header = TRUE,
#                       row.names = 1, 
#                       stringsAsFactors = FALSE)
# 
# sporulation <-  read.table(file = "../../data/4_phenotypes/cleaned_log_sporulation.tsv",
#                            sep = "\t", 
#                            header = TRUE,
#                            row.names = 1, 
#                            stringsAsFactors = FALSE)
# 
# toxin <-  read.table(file = "../../data/4_phenotypes/cleaned_log_toxin.tsv",
#                      sep = "\t", 
#                      header = TRUE,
#                      row.names = 1, 
#                      stringsAsFactors = FALSE)
# 
# cfe <- cbind(row.names(cfe), cfe)
# germ_tc_and_gly <- cbind(row.names(germ_tc_and_gly), germ_tc_and_gly)
# germ_tc <- cbind(row.names(germ_tc), germ_tc)
# growth <- cbind(row.names(growth), growth)
# sporulation <- cbind(row.names(sporulation), sporulation)
# toxin <- cbind(row.names(toxin), toxin)
# 
# colnames(cfe)[1] <- "ID"
# colnames(germ_tc_and_gly)[1] <- "ID"
# colnames(germ_tc)[1] <- "ID"
# colnames(growth)[1] <- "ID"
# colnames(sporulation)[1] <- "ID"
# colnames(toxin)[1] <- "ID"
# 
# pheno_df <- full_join(toxin, germ_tc_and_gly, by = "ID")
# pheno_df <- full_join(pheno_df, germ_tc, by = "ID")
# pheno_df <- full_join(pheno_df, growth, by = "ID")
# pheno_df <- full_join(pheno_df, sporulation, by = "ID")
# pheno_df <- full_join(pheno_df, cfe, by = "ID")
# 
# 
# 
# # Lets double check that the data looks like its plotted correctly
# # Is the 017 sample ~3 ln(tox) and either 0 or no data for germ in Tc? 
# key <- read.csv("../../data/6_ribotype/Hanna_in_vitro_data_plus_germination.txt") 
# key <- cbind(row.names(key), key)
# key <- as.data.frame(key)
# colnames(key)[1] <- "ID"
# key %>% 
#   filter(ribotype == "017")
# #  ID ribotype severe tcdA tcdB cdtA cdtB tcdC sporulation  toxin   CFE growth germ_tc_only
# # 1 DA00330      017      N    +    +    -    -  sc7        9633 18.410 44.47 0.0704            0
# # germ_tc_and_gly
# # 1           94.82
# # yup! log(18.4) = 2.9 
# # germ_tc_only = 0
# 
# # Also spot check the 013 sample
# key %>% 
#   filter(ribotype == "013")
# # ID ribotype severe tcdA tcdB cdtA cdtB tcdC sporulation toxin   CFE growth germ_tc_only germ_tc_and_gly
# # 1 Cd013      013      N    +    +    -    -  sc9        7550 0.240 51.32 0.0792          0.3           46.13
# 
# # tox log(0.24) == -1.4 looks good
# # germ tc log(0.3) = -1.2
# 
# 
# sum(is.na(key))
# 
# # get no data ids for each phenotype
# 
# sum(is.na(pheno_df$log_toxin)) # 0
# 
# sum(is.na(pheno_df$log_sporulation)) # 3
# spor_nd <- pheno_df %>% filter(is.na(log_sporulation)) %>% pull(ID) %>% as.character() # "Cd033" "Cd086" "VPI"  
# 
# sum(is.na(pheno_df$log_cfe)) # 3
# cfe_nd <- pheno_df %>% filter(is.na(log_cfe)) %>% pull(ID) %>% as.character()# "DA00484" "VPI"     "DA00315"
# 
# sum(is.na(pheno_df$log_growth)) # 4
# growth_nd <- pheno_df %>% filter(is.na(log_growth)) %>% pull(ID) %>% as.character() # "Cd054" "Cd182" "Cd072" "Cd093"
# 
# sum(is.na(pheno_df$log_germ_tc_and_gly)) #5
# germ_tc_and_gly_nd <- pheno_df %>% filter(is.na(log_germ_tc_and_gly)) %>% pull(ID) %>% as.character() #  "Cd054" "VPI"   "Cd182" "Cd072" "Cd093"  
# 
# sum(is.na(pheno_df$log_germ_tc)) # 6
# germ_tc_only_nd <- pheno_df %>% filter(is.na(log_germ_tc)) %>% pull(ID) %>% as.character() # "Cd054"   "VPI"     "Cd182"   "Cd072"   "Cd093"   "DA00330"
# 
# 
# tree <- read.tree(paste0(data_dir, "3_trees/log_toxin.tree"))
# ribotype <- read.table(paste0(data_dir, "6_ribotype/toxin_ribotype.tsv"), 
#                        sep = "\t", 
#                        stringsAsFactors = FALSE)       
# ribotype <- as.data.frame(cbind(row.names(ribotype), ribotype), 
#                           stringsAsFactors = FALSE)
# colnames(ribotype) <- c("ID", "Ribotype", "ribotype_color")
# 
# # Update ribotype colors
# ribotype$ribotype_color[ribotype$Ribotype == "017"] <- "darkorchid"
# ribotype$ribotype_color[ribotype$Ribotype == "027"] <- "blue"
# ribotype$ribotype_color[ribotype$Ribotype == "078-126"] <- "deepskyblue1"
# ribotype$ribotype_color[ribotype$Ribotype == "053-163"] <- "darkgreen"
# ribotype$ribotype_color[ribotype$Ribotype == "013"] <- "gold1"
# ribotype$ribotype_color[ribotype$Ribotype == "014-020"] <- "orange"
# ribotype$ribotype_color[ribotype$Ribotype == "001"] <- "red"
# ribotype$ribotype_color[ribotype$Ribotype == "003"] <- "brown"
# 
# 
# phenotype <- pheno_df
# phenotype$ID <- as.character(phenotype$ID)
# phenotype$log_toxin <- as.numeric(as.character(phenotype$log_toxin))
# phenotype$log_germ_tc_and_gly <- as.numeric(as.character(phenotype$log_germ_tc_and_gly))
# phenotype$log_germ_tc <- as.numeric(as.character(phenotype$log_germ_tc))
# phenotype$log_growth <- as.numeric(as.character(phenotype$log_growth))
# phenotype$log_sporulation <- as.numeric(as.character(phenotype$log_sporulation))
# phenotype$log_cfe <- as.numeric(as.character(phenotype$log_cfe))
# 
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
# sum(ribotype$ID != tree$tip.label) # check that these match
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
#   # geom_rect(fill = ribotype$ribotype_color, 
#   #           xmin = tox_x_min, 
#   #           xmax = tox_x_max,
#   #           ymin = min(as.numeric(as.character(phenotype$log_toxin))), 
#   #           ymax = max(as.numeric(as.character(phenotype$log_toxin))), 
#   #           alpha = 0.2)
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
# 
# 
# png(paste0(fig_dir,
#            "tree_w_ribo_and_pheno_barplots.png"),
#     units = "in",
#     width = 9, height = 10, res = 250)
# print(bar_plot_cfe %>% 
#         insert_bottom(tree_plot_wo_leg, height = 4) %>% 
#         insert_top(bar_plot_germ_tc) %>% 
#         insert_top(bar_plot_germ_tc_and_gly) %>% 
#         insert_top(bar_plot_growth) %>% 
#         insert_top(bar_plot_sporulation) %>% 
#         insert_top(bar_plot_toxin) %>% 
#         insert_bottom(ggpubr::as_ggplot(ribo_legend)))
# dev.off()
# 
# # Weirdly, I get an error if I try to put the plots left to right, but not if I plot them top to bottom. Oh well. 
# 
# 
# 
# # Which tree to use? The tree with the most samples; we'll have NA values for some samples 
# 
# 
