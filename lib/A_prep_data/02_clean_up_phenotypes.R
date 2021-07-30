# To be run after 1_clean_up_tree.R
source("data_prep_lib.R")

# Read in phenotypes
fqR <- read.table("../../data/4_phenotypes/fqR_pheno.tsv", sep = "\t", stringsAsFactors = FALSE)
cfe <- read.table("../../data/4_phenotypes/log_cfe_pheno.tsv", sep = "\t", stringsAsFactors = FALSE)
germ_tg <- read.table("../../data/4_phenotypes/log_germ_tc_and_gly_pheno.tsv", sep = "\t", stringsAsFactors = FALSE)
germ_t <- read.table("../../data/4_phenotypes/log_germ_tc_pheno.tsv", sep = "\t", stringsAsFactors = FALSE)
growth <- read.table("../../data/4_phenotypes/log_growth_pheno.tsv", sep = "\t", stringsAsFactors = FALSE)
spore <- read.table("../../data/4_phenotypes/log_sporulation_pheno.tsv", sep = "\t", stringsAsFactors = FALSE)
toxin <- read.table("../../data/4_phenotypes/log_toxin_pheno.tsv", sep = "\t", stringsAsFactors = FALSE)
severity <- read.table("../../data/4_phenotypes/severity_pheno.tsv", sep = "\t", stringsAsFactors = FALSE)

tree <- read.tree("../../data/3_trees/midpoint_rooted.tree")
key <- read.table("../../data/0_sample_list/genomes_to_use_in_analysis.tsv",
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE)

prep_phenotype(tree, fqR, key, "fqR") 
prep_phenotype(tree, cfe, key, "log_cfe") 
prep_phenotype(tree, germ_tg, key, "log_germ_tc_and_gly") 
prep_phenotype(tree, germ_t, key, "log_germ_tc") 
prep_phenotype(tree, growth, key, "log_growth") 
prep_phenotype(tree, spore, key, "log_sporulation") 
prep_phenotype(tree, toxin, key, "log_toxin") 
prep_phenotype(tree, severity, key, "severity") 
