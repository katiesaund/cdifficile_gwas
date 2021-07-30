# Goal: write a whole bunch of sbat files to submit and run simultaneously
#       These elastic net models will have many different options and then 
#       we'll compare the performance of models. 
# We're going to run for all of the 8 phenotypes
# We'll also include a negative controls for each of those phenotypes. 
# 3 negative controls: phenotype completely scrambeled on the tree. 
# Options for pyseer
# --wg enet fits an elastic net to all variants 
# --n-folds cross-validation (default 10-fold).
# --save-model saves the fitted model so it can be used for prediction.
# --cpu uses four cores efficiently during cross-valdation.
# --alpha controls the mixing between ridge regression and lasso regression. Above we have used a value of 1, which is lasso regression, selecting just a few variants. We can use a value closer to ridge regression if desired, which will select more variants with smaller effect sizes. "Sets the mixing between ridge regression (0) and lasso regression (1) in the above formula. Default is 0.0069 (closer to ridge regression)"
# --phenotypes  phenotype file
# --lineage-clusters has two effects. Cross-validation will be performed by leaving one strain out. This will usually take longer as there are more strains than folds, but may help reduce the number of lineage effects included. Also, training predition accuracy for each lineage will be reported, making it easier to see whether there are some parts of the data where the model is performing better. For binary phenotypes 𝑅2 can be difficult to interpret, so true/false positives/negatives are also reported.
# --sequence-reweighting has one further effect. Within each lineage, the weight 𝑤𝑖 given to each sample in the loss function
# --pres genotype presence/absence file

# --n-folds
n_folds_num <- 10

# --alpha
alpha_value_vec <- seq(0.01, 0.99, by = 0.245) 

# --lineage-clusters
# MLST_file_path <- "../../../data/MLST/cleaned_MLST_report.tsv"
# MLST_lineage_text <- paste0(" --lineage-clusters ",  MLST_file_path, " --sequence-reweighting ")
exclude_MLST_lineage_text <- " "
# In the future we could add clade or some other cluster information here instead of MLST
# lineage_vec <- c(MLST_lineage_text, exclude_MLST_lineage_text)
# names(lineage_vec) <- c("MLST", "no_lineage")
lineage_vec <- exclude_MLST_lineage_text
names(lineage_vec) <- "no_lineage"
# --pres
# # SNP only
# snp_path <- "../../data/snp_matrix/cleaned_parsed_bin_snp_matrix.tsv"
# # Indel only
# indel_path <- "" # once it gets made, it'll be here: "../data/indel_matrix/cleaned_parsed_bin_indel_matrix.tsv"
# # accessory genome only
# pangenome_path <- "../../data/pangenome_matrix/panaroo_gene_presence_absence.Rtab"
# SNP + Indel + accessory 
# paloc_mat_path <- "../../../data/15_pyseer/genotypes_for_pyseer/pyseer_cleaned_log_toxin_just_paloc_pan_and_rereferenced_snp_indel.tsv"
tcdb_mat_path <- "../../../data/15_pyseer/genotypes_for_pyseer/pyseer_cleaned_log_toxin_just_tcdB_pan_and_rereferenced_snp_indel.tsv"

# unique for each phenotype: cleaned_fqR_pan_and_rereferenced_snp_indel.tsv

# genotype_path_vec <- c(snp_path, indel_path, pangenome_path, snp_indel_pan_path)
genotype_path_vec <- c(tcdb_mat_path) #, tcdb_mat_path)
names(genotype_path_vec) <- "tcdb" # c("paloc", "tcdb")

# Phenotypes
ctrl_dir_path <- "../../../data/4_phenotypes/negative_controls/"
data_dir_path <-  "../../../data/15_pyseer/phenotypes_for_pyseer/"

true_pheno_paths <- paste0(data_dir_path, "cleaned_log_toxin.tsv")
names(true_pheno_paths) <- c("log_toxin")
ctrl_pheno_paths <- NULL
for (i in 1:length(true_pheno_paths)) {
  current_name <- names(true_pheno_paths)[i]
  ctrl_pheno_paths <- c(ctrl_pheno_paths, 
                        paste0(ctrl_dir_path, "neg_ctrl_", current_name, "_1.tsv"), 
                        paste0(ctrl_dir_path, "neg_ctrl_", current_name, "_2.tsv"), 
                        paste0(ctrl_dir_path, "neg_ctrl_", current_name, "_3.tsv"))
  names(ctrl_pheno_paths)[(length(ctrl_pheno_paths) - 2):length(ctrl_pheno_paths)] <- 
                             c(paste0("neg_ctrl_", current_name, "_1"), 
                              paste0("neg_ctrl_", current_name, "_2"),
                              paste0("neg_ctrl_", current_name, "_3"))
}

pheno_paths <- c(true_pheno_paths, ctrl_pheno_paths)

names(pheno_paths) <- c(names(true_pheno_paths), names(ctrl_pheno_paths))

# Covariates 
# Not including covariates for this analysis

counter <- 0

num_n_fold <- length(n_folds_num)
num_alpha <- length(alpha_value_vec)
num_lineage <- length(lineage_vec)
num_genotype <- length(genotype_path_vec)
num_phenotype <- length(pheno_paths)
# num_covariate <- length(covariate_vec)

for (i in 1:num_n_fold) {
  for (j in 1:num_alpha) {
    for (k in 1:num_lineage) {
      for (l in 1:num_genotype) {
        for (m in 1:num_phenotype) {
            counter <- counter + 1
            job_name <- paste0("EN", counter)
            
            model_name <- paste0("EN_risk_score_", n_folds_num[i], "_fold_", 
                                 alpha_value_vec[j], "_alpha_", 
                                 names(lineage_vec)[k], "_", 
                                 names(genotype_path_vec)[l], "_", 
                                 names(pheno_paths)[m]) #, "_", 
            out_name <- paste0(model_name, ".out")
            
            
            geno_mat_path <- genotype_path_vec[l]
            
            if (grepl("neg", names(pheno_paths)[m])) {
              geno_mat_path <- gsub("neg_ctrl_", "", geno_mat_path)
              geno_mat_path <- gsub("_[123]", "", geno_mat_path)
            }

            sbat_text <- paste0(
              "#!/bin/sh
#SBATCH --job-name=", 
              job_name, 
              "\n#SBATCH --output=", 
              out_name, 
              "\n#SBATCH --mail-user=katiephd@umich.edu
#SBATCH --mail-type=NONE
#SBATCH --export=ALL
#SBATCH --partition=standard
#SBATCH --account=esnitkin1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=500M --time=00:10:00
cd $SLURM_SUBMIT_DIR
echo $SLURM_SUBMIT_DIR
echo $SLURM_JOB_ID
pyseer --phenotypes ",  
              pheno_paths[m], 
              " --pres ", 
              geno_mat_path,
              " --wg enet --save-model  ../../../data/15_pyseer/tcdb_model_results/",  #tcdb_model_results
              model_name, 
              " --alpha ", 
              alpha_value_vec[j], 
              " --n-fold ", 
              n_folds_num[i],
              lineage_vec[k], 
              # covariate_vec[n],
              " > ../../../data/15_pyseer/tcdb_model_results/", #tcdb_model_results
              model_name, 
              ".txt ")
            file_name <- paste0(getwd(), "/", model_name, ".sbat")
            writeLines(sbat_text,
                       file_name,
                       sep = "\n") 
        }
      }
    }
  }
}


