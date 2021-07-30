tcdB_geno_mat <- data.table::fread(file = "../../../data/1_combine_genotype_matrices/combinedmat/just_tcdb_pan_and_rereferenced_snp_indel.tsv",
                                    sep = "\t", 
                                    data.table = FALSE, 
                                    stringsAsFactors = FALSE, 
                                    check.names = FALSE)

tcdB_geno_mat <- t(tcdB_geno_mat)
write.table(x = tcdB_geno_mat, 
            quote = FALSE,
            file = paste0("../../../data/15_pyseer/genotypes_for_pyseer/pyseer_cleaned_log_toxin_just_tcdB_pan_and_rereferenced_snp_indel.tsv"), 
            sep = "\t", 
            col.names = FALSE, 
            row.names = TRUE)
