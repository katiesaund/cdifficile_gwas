library(ape)
library(tidyverse)
library(phytools)
tcdb_prank_tree <- read.tree("../../data/3_trees/gene_trees/prank_tcdB.fasttree")
key <- read_tsv("../../data/0_sample_list/genomes_to_use_in_analysis.tsv")
tcdb_prank_tree$tip.label <- gsub("Cdif_", "", tcdb_prank_tree$tip.label) %>% 
  gsub("_genome", "", .) %>% 
  gsub("_", "", .)


toxin <-  read.table(file = "../../data/4_phenotypes/cleaned_log_toxin.tsv",
                     sep = "\t", 
                     header = TRUE,
                     row.names = 1, 
                     stringsAsFactors = FALSE)

for (i in 1:Ntip(tcdb_prank_tree)) {
  for (j in 1:nrow(key)) {
    if (tcdb_prank_tree$tip.label[i] == key$Unique_ID[j]) {
      tcdb_prank_tree$tip.label[i] <- key$Hanna_ID[j]
    }
  }
}

tips_to_drop <- tcdb_prank_tree$tip.label[!tcdb_prank_tree$tip.label %in% row.names(toxin)]

tcdb_prank_tree <- midpoint.root(tcdb_prank_tree)
tcdb_prank_tree <- drop.tip(tcdb_prank_tree, tips_to_drop)

write.tree(tcdb_prank_tree, "../../data/3_trees/gene_trees/cleaned_log_toxin_midpoint_root_tcdB_gene_prank.tree")

