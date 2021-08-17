# Subset roary multifasta to just the genes of interest from GWAS runs

library(ape)

data_dir <- "../../data/17_blast_overlap_genes//"

roary_fasta <- read.FASTA(file = paste0(data_dir, "/pan_genome_reference.fa"), 
                          type = "DNA")


# Toxin specific stuff 
toxin_roary_names <- read.csv(paste0(data_dir, "toxin_severity_sig_overlap_roary_names.csv"), stringsAsFactors = FALSE)
indices <- rep(NA, nrow(toxin_roary_names))

for (i in 1:nrow(toxin_roary_names)) {
  indices[i] <- which(grepl(paste0(toxin_roary_names[i, 1], "$"), names(roary_fasta)))
}
toxin_fasta <- roary_fasta[indices]

write.FASTA(toxin_fasta, file = paste0(data_dir, "toxin_severity_sig_overlap_indiv_hits.fna"))
write.csv(names(toxin_fasta),
          file = paste0(data_dir, "toxin_severity_sig_overlap_indiv_hits_full_name.csv"), row.names = FALSE, quote = FALSE)

