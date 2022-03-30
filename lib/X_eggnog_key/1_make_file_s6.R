library(tidyverse)
library(magrittr)
library(ape)

eggnog_path <- "../../data/18_eggnog/eggnog_run_2_vs_all/"
files <- list.files(path = eggnog_path, 
                    pattern = "*emapper.annotations.tsv",
                    recursive = TRUE)

expected_colnames <- c("query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", 
                       "max_annot_lvl", "COG_category", "Description", "Preferred_name", 
                       "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", 
                       "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
)

temp_tbl <- NULL
full_tbl <- NULL
for (i in 1:length(files)) {
  temp_tbl <- read_tsv(col_names = expected_colnames, 
                  paste0(eggnog_path, files[i]), 
                  comment = "#") 
  full_tbl <- rbind(full_tbl, temp_tbl)
}

stopifnot(identical(length(unique(full_tbl$query)), nrow(full_tbl)))

full_tbl %<>% 
  select(query, seed_ortholog, Preferred_name, Description, PFAMs) %>% 
  mutate(name1 = sub("_1$", "", query)) %>% 
  select(-query)

full_tbl %<>%
  mutate(seed_ortholog = sub(".*[.]", "", seed_ortholog))

pangenome <- read.FASTA("../../data/1_combine_genotype_matrices/pangenomemat/pan_genome_reference.fa")
pangenome_name_df <- as.data.frame(names(pangenome))
pangenome_name_df %<>% separate(col = `names(pangenome)`, sep = " ", into = c("name1", "name2"))

joined_tbl <-  left_join(pangenome_name_df, full_tbl, by = "name1")
colnames(joined_tbl) <- c("Genome_ID_plus_Prokka_ID", "Roary_ID", "Eggnog_seed_ortholog", "Eggnog_preferred_name", "Eggnog_description", "Eggnog_PFAMs")

write_tsv(joined_tbl,
          file = paste0(eggnog_path, "file_s6_pangenome_eggnog_key.tsv"))

