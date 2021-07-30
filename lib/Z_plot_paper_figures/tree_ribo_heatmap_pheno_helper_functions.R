# subfunctions to help make the tree + heatmap + phenotype figure

# Based on code from severity_GWAS/2020-12-09_agr_variants/lib/plot_low_agr_variants_only.R
# copied from thesis code

library(tidyverse)
library(phytools)
library(ape)
library(ggtree)
library(aplot)

read_and_midpoint_tree <- function(tree_path){
  # Read tree
  tree <- read.tree(tree_path)
  
  # Midpoint root tree
  tree <- midpoint.root(tree)
  
  # Return midpoint rooted tree
  return(tree)
}

read_format_ribotype <- function(ribo_path){
  ribotype <- read.table(ribo_path, 
                         sep = "\t", 
                         stringsAsFactors = FALSE)       
  ribotype <- as.data.frame(cbind(row.names(ribotype), ribotype), 
                            stringsAsFactors = FALSE)
  colnames(ribotype) <- c("ID", "Ribotype", "ribotype_color")
  return(ribotype)
}

create_ribotype_colors_obj <- function(ribo_df){
  # Ribo_df should be: 
  
  # ID      Ribotype  ribotype_color
  # PSM041	053-163	  aquamarine
  # PSM023	Other	    grey
  #	PSM059	Other	    grey
  
  if (!identical(c("ID", "Ribotype", "ribotype_color"), 
                   colnames(ribo_df))) {
    stop("ribotype dataframe incorrectly formatted")
  }
  
  ribotype_col <- unique(ribo_df$ribotype_color)
  names(ribotype_col) <- unique(ribo_df$Ribotype)
  return(ribotype_col)
  
  # TODO update colors to be ROYGBIV
}

read_format_pheno <- function(pheno_path, pheno_name){
  mat <- read.table(pheno_path,
                    sep = "\t",
                    row.names = 1,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
  mat <- as.matrix(mat)
  mat <- as.data.frame(cbind(row.names(mat), mat))
  colnames(mat) <- c("ID", pheno_name)
  mat$ID <- as.character(mat$ID)
  mat[, 2] <- as.numeric(as.character(mat[, 2]))
  return(mat)
}

read_prep_geno_for_heatmap <- function(geno_mat_path) {
  geno_mat <- data.table::fread(file = geno_mat_path,
                                    sep = "\t", 
                                    data.table = FALSE, 
                                    stringsAsFactors = FALSE, 
                                    check.names = FALSE)
  geno_mat <- as.matrix(geno_mat)
  row.names(geno_mat) <- geno_mat[, 1]
  geno_mat <- geno_mat[, 2:ncol(geno_mat)]
  storage.mode(geno_mat) <- "numeric"
  geno_mat <- t(geno_mat) # now samples in rows, variants in columns
  
  d <- dist(t(geno_mat), method = "euclidean") # distance matrix
  fit <- hclust(d, method = "ward.D")
  
  id_df <- as.data.frame(row.names(geno_mat), stringsAsFactors = FALSE)
  colnames(id_df)[1] <- "ID"
  
  geno_df <- as.data.frame(geno_mat, stringsAsFactors = FALSE)
  row.names(geno_df) <- NULL
  geno_df <- as.data.frame(sapply(geno_df, as.numeric))
  geno_df <- cbind(id_df, geno_df)
  
  geno_df_long <- geno_df %>%  
    pivot_longer(-ID, names_to = "locus", values_to = "presence")
  
  geno_df_long$ID <- factor(geno_df_long$ID, levels = rownames(geno_mat), labels =  rownames(geno_mat))
  geno_df_long$locus <- factor( geno_df_long$locus, levels = colnames(geno_mat)[fit$order],  labels = colnames(geno_mat)[fit$order])
  return(geno_df_long)
}

vert_tree_w_ribo <- function(tr, ribo_df, ribo_col){
  tr_plot <- ggtree(tr) %<+% ribo_df + 
    geom_tippoint(aes(color = Ribotype), size = .9) + 
    scale_color_manual(values = ribo_col) + 
    coord_flip()
  return(tr_plot)
}


vert_pheno_bar_plot_no_ribo <- function(pheno_df, pheno_name, plot_title=""){
  pheno_df <- pheno_df %>% 
    select(ID, all_of(pheno_name))
  bar_plot <- ggplot(pheno_df, aes(x = ID, 
                                   y = as.numeric(as.character(get(pheno_name))))) +
    geom_col() +
    ggtitle(plot_title) +
    theme_tree2() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  return(bar_plot)
}

vert_heatmap <- function(geno_long_df){
  heatmap.colours <- c("white","black")
  names(heatmap.colours) <- 0:1
  htmp <- ggplot(geno_long_df, aes(y = locus, x = ID)) +
    geom_tile(aes(fill = factor(presence))) + 
    scale_fill_manual(values = c("white", "black")) + 
    scale_y_discrete(position = "right") + 
    guides(fill = FALSE) + 
    theme_void()
  return(htmp)
}

