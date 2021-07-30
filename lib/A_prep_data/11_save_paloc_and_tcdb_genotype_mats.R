genotype <- data.table::fread(file = "../../data/1_combine_genotype_matrices/combinedmat/cleaned_log_toxin_pan_and_rereferenced_snp_indel.tsv",
                              sep = "\t", 
                              data.table = FALSE, 
                              stringsAsFactors = FALSE, 
                              check.names = FALSE)
genotype <- as.matrix(genotype)
row.names(genotype) <- genotype[, 1]
genotype <- genotype[, 2:ncol(genotype)]
print(genotype[1:2, 1:2])

locus_tags_for_tcd_genes <- c("CD630_06590", 
                              "CD630_06600", 
                              "CD630_06610",
                              "CD630_06630", 
                              "CD630_06640")

toxin_index <- as.logical(grepl(locus_tags_for_tcd_genes[1], colnames(genotype)) +
                            grepl(locus_tags_for_tcd_genes[2], colnames(genotype)) + 
                            grepl(locus_tags_for_tcd_genes[3], colnames(genotype)) + 
                            grepl(locus_tags_for_tcd_genes[4], colnames(genotype)) + 
                            grepl(locus_tags_for_tcd_genes[5], colnames(genotype)))

tcdb_index <- as.logical(grepl("CD630_06600", colnames(genotype)))

paloc_genotype <- genotype[, toxin_index]
tcdb_genotype <- genotype[, tcdb_index]

print(paloc_genotype[1:2, 1:2])
print(sum(toxin_index))
print(dim(genotype))
print(dim(paloc_genotype))
write.table(x = paloc_genotype, 
            quote = TRUE,
            file = "../../data/1_combine_genotype_matrices/combinedmat/just_paloc_pan_and_rereferenced_snp_indel.tsv",
            sep = "\t", 
            col.names = TRUE, 
            row.names = TRUE)


print(tcdb_genotype[1:2, 1:2])
print(sum(tcdb_index))
print(dim(genotype))
print(dim(tcdb_genotype))
write.table(x = tcdb_genotype, 
            quote = TRUE,
            file = "../../data/1_combine_genotype_matrices/combinedmat/just_tcdb_pan_and_rereferenced_snp_indel.tsv",
            sep = "\t", 
            col.names = TRUE, 
            row.names = TRUE)
