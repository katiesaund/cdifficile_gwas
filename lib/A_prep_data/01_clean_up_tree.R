source("data_prep_lib.R")

# unrooted tree
tree <- read.tree("../../data/3_trees/unrooted.tree")
# outgroup <- "ERR232398"

key <- read.table("../../data/0_sample_list/genomes_to_use_in_analysis.tsv",
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE)

tree <- clean_up_tree_names(tree, key)

# Mid point here so that we get the really long branch to PH53 included in the midpoint rooting step
tree <- phytools::midpoint.root(tree)

# Drop bad samples
tips_to_drop <- tree$tip.label[!tree$tip.label %in% key$Hanna_ID]
tree <- ape::drop.tip(tree, tips_to_drop)

# Save updated tree
write.tree(tree, file = "../../data/3_trees/midpoint_rooted.tree")

