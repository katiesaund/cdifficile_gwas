source("data_prep_lib.R")

phenotypes <- as.data.frame(matrix(NA, nrow = 8, ncol = 1))
colnames(phenotypes) <- "shorthand"
phenotypes[, 1] <- c("fqR", "log_cfe", "log_germ_tc_and_gly", "log_germ_tc", 
                     "log_growth", "log_sporulation", "log_toxin", "severity")

in_vitro_df <- read.csv("../../data/6_ribotype/Hanna_in_vitro_data_plus_germination.txt", stringsAsFactors = FALSE)
in_vitro_df <- in_vitro_df[, 1, drop = FALSE]
in_vitro_df$ribotype[grepl(pattern = "UM", in_vitro_df$ribotype)] <- "Other"
in_vitro_df$ribotype[grepl(pattern = "Unique", in_vitro_df$ribotype)] <- "Other"


for (i in 1:nrow(phenotypes)) {
  tree <- read.tree(paste0("../../data/3_trees/", phenotypes$shorthand[i], ".tree"))
  ribo_df <- in_vitro_df[row.names(in_vitro_df) %in% tree$tip.label, , drop = FALSE]
  
  ribo_df <- ribo_df[match( tree$tip.label, row.names(ribo_df)), , drop = FALSE]
  
  identical(row.names(ribo_df), tree$tip.label)
  
  ribo_df <- cbind(ribo_df, ribo_df)
  colnames(ribo_df) <- c("ribotype", "color")
  ribo_df$color[grepl(pattern = "078-126", ribo_df$color)] <- "darkorchid"
  ribo_df$color[grepl(pattern = "017", ribo_df$color)] <- "blue"
  ribo_df$color[grepl(pattern = "027", ribo_df$color)] <- "deepskyblue1"
  ribo_df$color[grepl(pattern = "013", ribo_df$color)] <- "darkgreen"
  ribo_df$color[grepl(pattern = "053-163", ribo_df$color)] <- "gold1"
  ribo_df$color[grepl(pattern = "003", ribo_df$color)] <- "orange"
  ribo_df$color[grepl(pattern = "001", ribo_df$color)] <- "red"
  ribo_df$color[grepl(pattern = "014-020", ribo_df$color)] <- "brown"
  ribo_df$color[grepl(pattern = "Other", ribo_df$color)] <- "grey"
  
  write.table(ribo_df,
              paste0("../../data/6_ribotype/", phenotypes$shorthand[i], "_ribotype.tsv"), 
              sep = "\t",
              row.names = TRUE, 
              col.names = TRUE)

  # add these files for plotting ribotype within hogwash
  write.table(ribo_df[, 1, drop=FALSE], 
              paste0("../../data/6_ribotype/", phenotypes$shorthand[i], "_ribotype_no_colors.tsv"),
              sep = "\t",
              row.names = TRUE,
              col.names = TRUE)
}
