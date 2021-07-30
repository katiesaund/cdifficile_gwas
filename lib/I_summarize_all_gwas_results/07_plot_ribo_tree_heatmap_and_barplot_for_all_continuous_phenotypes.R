# currently only works interactively, not from command line due to pivot_longer 

options(scipen = 999) # turn off scientific notation
library(tidyverse)
library(dplyr)
library(phytools)
library(ape)
library(ggtree)
library(aplot)
source("gwas_analysis_lib.R")

data_dir <- "../../data/"
fig_dir <- "../../figures/gwas_results/"
phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")

plot_ribo_tree_and_heatmap("fqR", "phyc", "indiv", phenotypes)
plot_ribo_tree_and_heatmap("fqR", "synchronous", "indiv", phenotypes)
plot_ribo_tree_and_heatmap("severity", "phyc", "indiv", phenotypes)
plot_ribo_tree_and_heatmap("severity", "synchronous", "indiv", phenotypes)
plot_ribo_tree_and_heatmap("log_cfe", "continuous", "indiv", phenotypes)
plot_ribo_tree_and_heatmap("log_germ_tc_and_gly", "continuous", "indiv", phenotypes)
plot_ribo_tree_and_heatmap("log_germ_tc", "continuous", "indiv", phenotypes)
plot_ribo_tree_and_heatmap("log_growth", "continuous", "indiv", phenotypes)
plot_ribo_tree_and_heatmap("log_sporulation", "continuous", "indiv", phenotypes)
plot_ribo_tree_and_heatmap("log_toxin", "continuous", "indiv", phenotypes)

plot_ribo_tree_and_heatmap("fqR", "phyc", "group", phenotypes)
plot_ribo_tree_and_heatmap("fqR", "synchronous", "group", phenotypes)
plot_ribo_tree_and_heatmap("severity", "phyc", "group", phenotypes)
plot_ribo_tree_and_heatmap("severity", "synchronous", "group", phenotypes)
plot_ribo_tree_and_heatmap("log_cfe", "continuous", "group", phenotypes)
plot_ribo_tree_and_heatmap("log_germ_tc_and_gly", "continuous", "group", phenotypes)
plot_ribo_tree_and_heatmap("log_germ_tc", "continuous", "group", phenotypes)
plot_ribo_tree_and_heatmap("log_growth", "continuous", "group", phenotypes)
plot_ribo_tree_and_heatmap("log_sporulation", "continuous", "group", phenotypes)
plot_ribo_tree_and_heatmap("log_toxin", "continuous", "group", phenotypes)
