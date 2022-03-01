library(tidyverse)
library(patchwork)
library(dplyr)
library(phytools)
library(ape)
library(ggtree)
library(aplot)
library(ggpubr)
library(ggplotify)

# 
y_width <- 15
y_fontsize <- 9
y_num_size <- 6
y_angle <- 0 
y_title_angle <- 0

# Tree
tree <- read.tree("../../data/3_trees/log_toxin.tree")

# ribotype
ribotype <- read.table("../../data/6_ribotype/log_toxin_ribotype.tsv", 
                       sep = "\t", 
                       stringsAsFactors = FALSE)
ribotype <- as.data.frame(cbind(row.names(ribotype), ribotype), 
                          stringsAsFactors = FALSE)
colnames(ribotype) <- c("ID", "Ribotype", "ribotype_color")
ribotype_colors <- unique(ribotype$ribotype_color)
names(ribotype_colors) <- unique(ribotype$Ribotype)

# Tree plot ----
tree_plot_wo_leg <- ggtree(tree) %<+% ribotype + 
  geom_tippoint(aes(color = Ribotype)) + 
  scale_color_manual(values = ribotype_colors) + 
  coord_flip() + 
  theme(legend.position = "none")

tree_plot_w_leg <- ggtree(tree) %<+% ribotype + 
  geom_tippoint(aes(color = Ribotype)) + 
  scale_color_manual(values = ribotype_colors) + 
  coord_flip() + 
  theme(legend.position = c(0.65, 2.0),
        legend.direction = "horizontal")

# legend
ribo_legend <- ggpubr::get_legend(tree_plot_w_leg)

# Find ND points
key <- read.csv("../../data/6_ribotype/Hanna_in_vitro_data_plus_germination.txt") 
key <- cbind(row.names(key), key)
key <- as.data.frame(key)
colnames(key)[1] <- "ID"
key <- key %>% filter(!ID %in% c("630", "R20291", "ATCC_43255", "VPI10463_asm"))

sev_nd <- key %>% filter(!severe %in% c("N", "Y")) %>% pull(ID) %>% as.character(.)
toxin_nd <- key %>% filter(toxin == "ND") %>% pull(ID) %>% as.character(.)
germ_tc_nd <- key %>% filter(is.na(germ_tc_only)) %>% pull(ID) %>% as.character(.)
germ_tc_and_gly_nd <- key %>% filter(is.na(germ_tc_and_gly)) %>% pull(ID) %>% as.character(.)
growth_nd <- key %>% filter(growth == "ND") %>% pull(ID) %>% as.character(.)
sporulation_nd <- key %>% filter(sporulation == "ND" | sporulation == 0) %>% pull(ID) %>% as.character(.)
# We lost two sporulation data points by doing the log(0) so we need to add those as well

cfe_nd <- key %>% filter(CFE == "ND") %>% pull(ID) %>% as.character(.)

# No missing toxin data
identical(sev_nd, germ_tc_nd)
identical(sev_nd, germ_tc_and_gly_nd)
growth_nd
cfe_nd
key %>% filter(ID %in% sporulation_nd)

nd_text_y <- 0.3
nd_text_size <- 1.0
                
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
           fill = "brown")

# Germ. in Tc Bar Plot ----
germ_tc <- read.table("../../data/4_phenotypes/cleaned_log_germ_tc.tsv",
                    sep = "\t",
                    row.names = 1,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
germ_tc <- as.matrix(germ_tc)
germ_tc <- as.data.frame(cbind(row.names(germ_tc), germ_tc))
colnames(germ_tc) <- c("ID", "germ_tc")
germ_tc$ID <- as.character(germ_tc$ID)
germ_tc$germ_tc <- as.numeric(as.character(germ_tc$germ_tc))

germ_tc_plot <- 
  ggplot(germ_tc, aes(ID, germ_tc)) + 
  geom_col() + 
  ggtitle("Germination in TCA") + 
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(angle = y_angle, size = y_num_size), 
        axis.title.y = element_text(size = y_fontsize, angle = y_title_angle)) + 
  ylab(lapply(strwrap("Germination ln(%)",
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
  geom_text(data = NULL, x = 5,  y = nd_text_y, label = "ND", size = nd_text_size) + # Cd054 078-126
  geom_text(data = NULL, x = 25, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd093 027
  geom_text(data = NULL, x = 59, y = nd_text_y, label = "ND", size = nd_text_size) + # VPI 003
  geom_text(data = NULL, x = 86, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd182 014-020
  geom_text(data = NULL, x = 97, y = nd_text_y, label = "ND", size = nd_text_size)   # Cd072 014-020


# log_germ_tc_and_gly Bar Plot ----
germ_tc_and_gly <- read.table("../../data/4_phenotypes/cleaned_log_germ_tc_and_gly.tsv",
                    sep = "\t",
                    row.names = 1,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
germ_tc_and_gly <- as.matrix(germ_tc_and_gly)
germ_tc_and_gly <- as.data.frame(cbind(row.names(germ_tc_and_gly), germ_tc_and_gly))
colnames(germ_tc_and_gly) <- c("ID", "germ_tc_and_gly")
germ_tc_and_gly$ID <- as.character(germ_tc_and_gly$ID)
germ_tc_and_gly$germ_tc_and_gly <- as.numeric(as.character(germ_tc_and_gly$germ_tc_and_gly))

germ_tc_and_gly_plot <- 
  ggplot(germ_tc_and_gly, aes(ID, germ_tc_and_gly)) + 
  geom_col() + 
  ggtitle("Germination in TCA&Gly") + 
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(angle = y_angle, size = y_num_size), 
        axis.title.y = element_text(size = y_fontsize, angle = y_title_angle)) + 
  ylab(lapply(strwrap("Germination ln(%)",
                      width = y_width, 
                      simplify = FALSE),
              paste, collapse = "\n")) + 
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

# cfe Bar Plot ----
cfe <- read.table("../../data/4_phenotypes/cleaned_log_cfe.tsv",
                  sep = "\t",
                  row.names = 1,
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  check.names = FALSE)
cfe <- as.matrix(cfe)
cfe <- as.data.frame(cbind(row.names(cfe), cfe))
colnames(cfe) <- c("ID", "cfe")
cfe$ID <- as.character(cfe$ID)
cfe$cfe <- as.numeric(as.character(cfe$cfe))

cfe_plot <- 
  ggplot(cfe, aes(ID, cfe)) + 
  geom_col() + 
  ggtitle("Spore viability") + 
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(angle = y_angle, size = y_num_size), 
        axis.title.y = element_text(size = y_fontsize, angle = y_title_angle)) + 
  ylab(lapply(strwrap("Spore viability ln(%)",
                      width = y_width, 
                      simplify = FALSE),
              paste, collapse = "\n")) + 
  geom_text(data = NULL, x = 59, y = nd_text_y, label = "ND", size = nd_text_size) + # VPI 003
  geom_text(data = NULL, x = 74, y = nd_text_y, label = "ND", size = nd_text_size) + # DA00315 unique
  geom_text(data = NULL, x = 1, y = nd_text_y, label = "ND", size = nd_text_size) + # DA00484 unique
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

# sporulation Bar Plot ----
sporulation <- read.table("../../data/4_phenotypes/cleaned_log_sporulation.tsv",
                          sep = "\t",
                          row.names = 1,
                          header = TRUE,
                          stringsAsFactors = FALSE,
                          check.names = FALSE)
sporulation <- as.matrix(sporulation)
sporulation <- as.data.frame(cbind(row.names(sporulation), sporulation))
colnames(sporulation) <- c("ID", "sporulation")
sporulation$ID <- as.character(sporulation$ID)
sporulation$sporulation <- as.numeric(as.character(sporulation$sporulation))

sporulation_plot <- 
  ggplot(sporulation, aes(ID, sporulation)) + 
  geom_col() + 
  ggtitle("Total spores") + 
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(angle = y_angle, size = y_num_size), 
        axis.title.y = element_text(size = y_fontsize, angle = y_title_angle)) + 
  ylab(lapply(strwrap("Spore count ln(CFU/ml)",
                      width = y_width, 
                      simplify = FALSE),
              paste, collapse = "\n")) + 
  geom_text(data = NULL, x = 59, y = nd_text_y, label = "ND", size = nd_text_size) + # VPI 003
  geom_text(data = NULL, x = 55, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd086 078-126
  geom_text(data = NULL, x = 4, y = nd_text_y, label = "ND", size = nd_text_size) + # Cd033 Other
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

# growth Bar Plot ----
growth <- read.table("../../data/4_phenotypes/cleaned_log_growth.tsv",
                     sep = "\t",
                     row.names = 1,
                     header = TRUE,
                     stringsAsFactors = FALSE,
                     check.names = FALSE)
growth <- as.matrix(growth)
growth <- as.data.frame(cbind(row.names(growth), growth))
colnames(growth) <- c("ID", "growth")
growth$ID <- as.character(growth$ID)
growth$growth <- as.numeric(as.character(growth$growth))

growth_plot <- 
  ggplot(growth, aes(ID, growth)) + 
  geom_col() + 
  ggtitle("Maximum growth rate") + 
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(angle = y_angle, size = y_num_size), 
        axis.title.y = element_text(size = y_fontsize, angle = y_title_angle)) + 
  ylab(lapply(strwrap("Max growth rate ln(OD/hour)",
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
  geom_text(data = NULL, x = 5,  y = -nd_text_y, label = "ND", size = nd_text_size) + # Cd054 078-126
  geom_text(data = NULL, x = 25, y = -nd_text_y, label = "ND", size = nd_text_size) + # Cd093 027
  geom_text(data = NULL, x = 86, y = -nd_text_y, label = "ND", size = nd_text_size) + # Cd182 014-020
  geom_text(data = NULL, x = 97, y = -nd_text_y, label = "ND", size = nd_text_size)   # Cd072 014-020

# Put it all together
combined_plot <- toxin_plot %>% 
  insert_bottom(tree_plot_wo_leg, height = 4) %>% 
  insert_top(cfe_plot) %>% 
  insert_top(sporulation_plot) %>% 
  insert_top(growth_plot) %>% 
  insert_top(germ_tc_and_gly_plot) %>% 
  insert_top(germ_tc_plot) %>% 
  insert_top(sev_plot, height = 0.3) %>% 
  insert_bottom(ggpubr::as_ggplot(ribo_legend))
ggsave(plot = combined_plot,
       filename =  "../../figures/Fig1.png",
       width = 8,
       height = 8,
       units = "in")
ggsave(plot = combined_plot,
       filename =  "../../figures/Fig1.pdf",
       width = 8,
       height = 8,
       units = "in")
