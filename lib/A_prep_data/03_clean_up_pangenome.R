source("data_prep_lib.R")

pangenome_mat <- 
  read.table("../../data/1_combine_genotype_matrices/pangenomemat/gene_presence_absence.Rtab",
             header = TRUE, 
             row.names = 1, 
             sep = "\t", 
             quote = "", 
             stringsAsFactors = FALSE)

key <- read.table("../../data/0_sample_list/genomes_to_use_in_analysis.tsv",
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE)

pangenome_mat <- clean_up_mat_names(pangenome_mat, key)  
pangenome_mat <- pangenome_mat[, colnames(pangenome_mat) %in% key$Hanna_ID, drop = FALSE]

write.table(pangenome_mat, 
            sep = "\t", 
            quote = FALSE,
            row.names = TRUE, 
            col.names = TRUE, 
            file = "../../data/1_combine_genotype_matrices/pangenomemat/cleaned_names_pangenomemat.tsv")
