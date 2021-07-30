# Goal: Save a csv with genotype name, P-value, epsilon, and rank for each phenotype

library(seqinr)
library(tidyverse)
source("../I_summarize_all_gwas_results/gwas_analysis_lib.R") 
phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")

load("../../data/1_combine_genotype_matrices/snpmat/snp_annots.RData") # snp_annotations
load("../../data/1_combine_genotype_matrices/indelmat/indel_annots.RData") # indel_annotations

# Combine indel and snp annotation matrices (kinda messy) so we can easily grab annotations below ---- 
# indel_annotations <- indel_annotations %>% rename("raw_rownames" = "full_annots")

snp_only_colnames <- colnames(snp_annotations)[!colnames(snp_annotations) %in% colnames(indel_annotations)]
indel_only_colnames <- colnames(indel_annotations)[!colnames(indel_annotations) %in% colnames(snp_annotations)]
indel_orig_colnames <- colnames(indel_annotations)
snp_orig_colnames <- colnames(snp_annotations)

indel_annotations <- cbind(indel_annotations, rep(NA, nrow(indel_annotations)))
colnames(indel_annotations) <- c(indel_orig_colnames, snp_only_colnames)

snp_annotations <- cbind(snp_annotations, rep(NA, nrow(snp_annotations)), rep(NA, nrow(snp_annotations)), rep(NA, nrow(snp_annotations)))
colnames(snp_annotations) <- c(snp_orig_colnames, indel_only_colnames)

# Reorder indel_annotations to same order and snp_annotations
indel_annotations <- indel_annotations[, match(colnames(snp_annotations), colnames(indel_annotations))]
identical(colnames(indel_annotations), colnames(snp_annotations))

annotation_df <- rbind(snp_annotations, indel_annotations)
# ----- 


# Make a "cleanish" version of annotation for locus tags. It won't be perfect, though. ---- 
group_annotation_df <- annotation_df %>% 
  select(locus_tag, annotation_1, annotation_2) %>% 
  unique()

row_to_drop <- NULL
for (i in 1:nrow(group_annotation_df)) {
  if (group_annotation_df$locus_tag[i] == group_annotation_df$annotation_1[i] & 
      group_annotation_df$locus_tag[i] == group_annotation_df$annotation_2[i]) {
    row_to_drop <- c(row_to_drop, i)
  }
}
group_annotation_df <- group_annotation_df[-row_to_drop, ]

first_appearance_index <- NULL
unique_locus_tags <- unique(group_annotation_df$locus_tag)
for (i in 1:length(unique_locus_tags)) {
  locus_tag_index <- NULL
  for (j in 1:nrow(group_annotation_df)) {
    if (unique_locus_tags[i] == group_annotation_df$locus_tag[j]) {
      locus_tag_index <- c(locus_tag_index, j)
    }
  }
  first_appearance_index <- c(first_appearance_index, locus_tag_index[1])
}
group_annotation_df <- group_annotation_df[first_appearance_index, ]
# ----- 

# Add roary information for individual locus stuff ----
pangenome_dna <- seqinr::read.fasta(file = "../../data/1_combine_genotype_matrices/pangenomemat/pan_genome_reference.fa", seqtype = "DNA")
full_names <- unlist(getAnnot(pangenome_dna))
split_names <- unlist(unlist(full_names) %>% str_split(., " "))
annot <- split_names[!grepl(">", split_names)]
prokka <- split_names[grepl(">", split_names)]
annot_prokka_key <- as.data.frame(cbind(annot, prokka))
annot_prokka_key$prokka <- gsub(">", "", annot_prokka_key$prokka)
colnames(annot_prokka_key)[2] <- "roary_prokka_ID"
# -----


for (i in 1:nrow(phenotypes)) {
  temp_shorthand <- phenotypes$shorthand[i]
  
  indiv_data <- load_indiv_hogwash_rda_for_pheno(temp_shorthand, "../../data/5_hogwash/")
  group_data <- load_post_ar_hogwash_rda_for_pheno(temp_shorthand, "../../data/5_hogwash/")

  
  if (phenotypes$bin_cont[i] == "continuous") {
    indiv_ranked_df <- rank_by_epsilon_and_p(indiv_data$hogwash_continuous)
    group_ranked_df <- rank_by_epsilon_and_p(group_data$hogwash_continuous)
    
    indiv_ranked_df <- left_join(indiv_ranked_df, annotation_df, 
                                 by = c("Genotype" = "raw_rownames"))
    indiv_ranked_df <- indiv_ranked_df %>% 
      select(Genotype, `-ln(FDR P-value)`, Epsilon, Rank, label, variant_type, 
             ref, var, aa_change, annotation_1, annotation_2, intergenic) %>% 
      rename("reference_allele"= "ref", 
             "variant_allele" = "var", 
             "amino_acid_change" = "aa_change", 
             "intergenic_variant_logical" = "intergenic", 
             "variant_type_1" = "label",  
             "variant_type_2" = "variant_type")
    indiv_ranked_df <- left_join(indiv_ranked_df, 
                                 annot_prokka_key, by = c("Genotype" = "annot"))
    
    group_ranked_df <- left_join(group_ranked_df, group_annotation_df, 
                                      by = c("Genotype" = "locus_tag"))    
    write_tsv(indiv_ranked_df, 
              paste0("../../data/13_summaries/", temp_shorthand, "_indiv_genotypes_ranked_by_epsilon_then_p.tsv"), 
              col_names = TRUE)
    write_tsv(group_ranked_df, 
              paste0("../../data/13_summaries/", temp_shorthand, "_group_genotypes_ranked_by_epsilon_then_p.tsv"), 
              col_names = TRUE)
    
  } else {
    indiv_ranked_phyc_df <- rank_by_epsilon_and_p(indiv_data$hogwash_phyc)
    indiv_ranked_sync_df <- rank_by_epsilon_and_p(indiv_data$hogwash_synchronous)
    
    # combine with annots: gene name, annotations, variant type and variant
    
    indiv_ranked_phyc_df <- left_join(indiv_ranked_phyc_df, annotation_df, 
                                      by = c("Genotype" = "raw_rownames"))
    indiv_ranked_phyc_df <- indiv_ranked_phyc_df %>% 
      select(Genotype, `-ln(FDR P-value)`, Epsilon, Rank, label, variant_type, 
             ref, var, aa_change, annotation_1, annotation_2, intergenic) %>% 
      rename("reference_allele"= "ref", 
             "variant_allele" = "var", 
             "amino_acid_change" = "aa_change", 
             "intergenic_variant_logical" = "intergenic", 
             "variant_type_1" = "label",  
             "variant_type_2" = "variant_type")
    indiv_ranked_phyc_df <- left_join(indiv_ranked_phyc_df, 
                                 annot_prokka_key, by = c("Genotype" = "annot"))
    
    indiv_ranked_sync_df
    indiv_ranked_sync_df <- left_join(indiv_ranked_sync_df, annotation_df,
                                      by = c("Genotype" = "raw_rownames"))
    indiv_ranked_sync_df <- indiv_ranked_sync_df %>% 
      select(Genotype, `-ln(FDR P-value)`, Epsilon, Rank, label, variant_type, 
             ref, var, aa_change, annotation_1, annotation_2, intergenic) %>% 
      rename("reference_allele"= "ref", 
             "variant_allele" = "var", 
             "amino_acid_change" = "aa_change", 
             "intergenic_variant_logical" = "intergenic", 
             "variant_type_1" = "label",  
             "variant_type_2" = "variant_type")
    indiv_ranked_sync_df <- left_join(indiv_ranked_sync_df, 
                                      annot_prokka_key, 
                                      by = c("Genotype" = "annot"))
    
    
    group_ranked_phyc_df <- rank_by_epsilon_and_p(group_data$hogwash_phyc)
    group_ranked_sync_df <- rank_by_epsilon_and_p(group_data$hogwash_synchronous)
    
    group_ranked_phyc_df <- left_join(group_ranked_phyc_df, group_annotation_df, 
                                      by = c("Genotype" = "locus_tag"))
    group_ranked_sync_df <- left_join(group_ranked_sync_df, group_annotation_df, 
                                      by = c("Genotype" = "locus_tag"))
    
    
    write_tsv(indiv_ranked_phyc_df, 
              paste0("../../data/13_summaries/", temp_shorthand, "_phyc_indiv_genotypes_ranked_by_epsilon_then_p.tsv"), 
              col_names = TRUE)
    write_tsv(indiv_ranked_sync_df, 
              paste0("../../data/13_summaries/", temp_shorthand, "_synchronous_indiv_genotypes_ranked_by_epsilon_then_p.tsv"), 
              col_names = TRUE)
    write_tsv(group_ranked_phyc_df, 
              paste0("../../data/13_summaries/", temp_shorthand, "_phyc_group_genotypes_ranked_by_epsilon_then_p.tsv"), 
              col_names = TRUE)
    write_tsv(group_ranked_sync_df, 
              paste0("../../data/13_summaries/", temp_shorthand, "_synchronous_group_genotypes_ranked_by_epsilon_then_p.tsv"), 
              col_names = TRUE)
  }
}