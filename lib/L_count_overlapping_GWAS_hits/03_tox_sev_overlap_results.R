load("../../data/13_summaries/overlapping_gwas_hits.rda")

accessory_genes <- unique(intersect_df_list$`severity-log_toxin`$Genotype)[1:6]

# Let's extract the sequences of these 6 genes, save them, and then blast them in data/17_blast_overlap_genes

write.csv(accessory_genes,
          "../../data/17_blast_overlap_genes/toxin_severity_sig_overlap_roary_names.csv", 
          quote = FALSE, 
          row.names = FALSE)
