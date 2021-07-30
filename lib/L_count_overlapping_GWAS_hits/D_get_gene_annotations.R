# Take in normally named genes and the best blast results. Look up their gene annotations if possible. 

#I might be able to skip blasting entirely 
# Check out this directory: 
# Sequence_data/Project_Hanna_collections/assembly_annotation/prokka_genbank_files
library(ape)
library(tidyverse)

data_dir <- "../../data/8_roary_hits_of_interest/"

roary_full_names <- read.csv(file = paste0(data_dir, "/log_toxin/toxin_indiv_hits_full_name.csv"),
                             stringsAsFactors = FALSE)

normal_names <- read.csv(file = paste0(data_dir, "/log_toxin/indiv_normal_names.csv"),
                             stringsAsFactors = FALSE)

# grap annotation from annotation file

#annot_dir <- "../../../../../Sequence_data/Project_Hanna_collections/assembly_annotation/prokka_genbank_files/"

gene_annotation_df <- as.data.frame(matrix(NA, ncol = 4, nrow = sum(nrow(normal_names), nrow(roary_full_names))))
colnames(gene_annotation_df) <- c("original_name", "roary_id", "prokka_annotation", "locus_tag")
gene_annotation_df$original_name <- c(roary_full_names$x, normal_names$genotype)

# for (i in 1:nrow(roary_full_names)) {
#   if (!grepl("Cdif_630", roary_full_names[i, 1])) {
#     group_name <- strsplit(roary_full_names[i, 1], " ")[[1]][2]
#     temp <- strsplit(roary_full_names[i, 1], " ")[[1]][1]
#     sample_name <- strsplit(temp, "_")[[1]][1]
#     prokka_name <- paste0(strsplit(temp, "_")[[1]][2], "_", strsplit(temp, "_")[[1]][3])
#     
#     print(sample_name)
#     print(prokka_name)
#     print(group_name)
#     # annot_df <- read.table(file = paste0(annot_dir, "/", sample_name, "/", sample_name, ".annot"), 
#     #                         sep = "\t", stringsAsFactors = FALSE, header = FALSE)
#     annot_df <- read_tsv(file = paste0(annot_dir, "/", sample_name, "/", sample_name, ".annot"), col_names = FALSE)
#     print(head(annot_df))
#     gene_annotation_df$prokka_annotation[i] <- annot_df %>% filter(X2 == prokka_name) %>% pull(X4)
#   }
# }


for (i in 1:nrow(normal_names)) {
  df_index <- nrow(roary_full_names) + i 
  gene_id <- gsub("_.*", "", normal_names$genotype[i])
  gff_df <- read.gff(file = "../../data/8_roary_hits_of_interest/trimmed_gff/Cdif_630_trimmed.gff")
  search_term <- paste0("gene=", gene_id)
  attribute <- gff_df$attributes[grepl(search_term, gff_df$attributes)]
  product <- gsub(".*product=", "", attribute)
  locus_tag <- gsub(".*locus_tag=", "", attribute) %>% gsub("product=.*", "", .)  %>% gsub(";", "", .)
  gene_annotation_df$prokka_annotation[df_index] <- product
  gene_annotation_df$locus_tag[df_index] <- locus_tag
  gene_annotation_df$roary_id[df_index] <- gene_id
}


for (i in 1:nrow(roary_full_names)) {
  if (grepl("Cdif_630", roary_full_names[i, 1])) {
    group_name <- strsplit(roary_full_names[i, 1], " ")[[1]][2]
    temp <- strsplit(roary_full_names[i, 1], " ")[[1]][1]
    sample_name <- paste0(strsplit(temp, "_")[[1]][1], "_", strsplit(temp, "_")[[1]][2])
    prokka_name <- paste0(strsplit(temp, "_")[[1]][3], "_", strsplit(temp, "_")[[1]][4])
  } else {
    group_name <- strsplit(roary_full_names[i, 1], " ")[[1]][2]
    temp <- strsplit(roary_full_names[i, 1], " ")[[1]][1]
    sample_name <- strsplit(temp, "_")[[1]][1]
    prokka_name <- paste0(strsplit(temp, "_")[[1]][2], "_", strsplit(temp, "_")[[1]][3])
  }
  gff_df <- read.gff(file = paste0("../../data/8_roary_hits_of_interest/trimmed_gff/", sample_name, "_trimmed.gff"))
  search_term <- paste0("ID=", sample_name, "_", prokka_name, ";")
  attribute <- gff_df$attributes[grepl(search_term, gff_df$attributes)]
  product <- gsub(".*product=", "", attribute)
  locus_tag <- gsub(".*locus_tag=", "", attribute) %>% gsub("product=.*", "", .)  %>% gsub(";", "", .)
  gene_annotation_df$prokka_annotation[i] <- product
  gene_annotation_df$locus_tag[i] <- locus_tag
  gene_annotation_df$roary_id[i] <- group_name
}

write.csv(gene_annotation_df, file = "../../data/8_roary_hits_of_interest/log_toxin/indiv_roary_hits_with_annotations.csv",
          quote = FALSE, row.names = FALSE)