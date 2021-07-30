# Write sbat -- this code is WAY too messy. time to refactor.

make_sbatch <- function(temp_name, temp_partition, temp_memory, temp_time, command){
  fname <- paste0(getwd(), "/", "03_hogwash_", temp_name, ".sbat")
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=", temp_name),
               paste0("#SBATCH --output=", temp_name, ".out"),
               "#SBATCH --mail-user=katiephd@umich.edu",
               "#SBATCH --mail-type=NONE",
               "#SBATCH --export=ALL",
               paste0("#SBATCH --partition=", temp_partition),
               "#SBATCH --account=esnitkin1",
               paste0("#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=", temp_memory, " --time=", temp_time),
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               command),
             fname,
             sep = "\n")
}

write_command <- function(temp_pheno,
                           temp_geno,
                           temp_tree,
                           temp_name,
                           temp_outdir,
                           temp_perm,
                           temp_fdr,
                           temp_bootstrap,
                           temp_test,
                           temp_group_method,
                           temp_key, 
			   temp_strain){
  command <- paste(paste0("Rscript run_hogwash_sbatch.R "),
                   temp_pheno,
                   temp_geno,
                   temp_tree,
                   temp_name,
                   temp_outdir,
                   temp_perm,
                   temp_fdr,
                   temp_bootstrap,
		   temp_test,
                   temp_group_method,
                   temp_strain,
                   temp_key,
                   sep = " ")
  return(command)
}

write_hogwash_sbat <- function(perm = 10000, 
                               group_fdr = 0.001,
                               indiv_fdr = 0.3,
                               bootstrap = 0.95,
                               binary_memory = "47G", 
                               continuous_memory = "60G", 
                               group_memory = "20G",
                               continuous_time = "240:00:00",
                               binary_time = "120:00:00", 
                               gene_key_path = "../../data/2_gene_keys/", 
                               data_path = "../../data/") {
  bin_phenos <- c("fqR", "severity")
  cont_phenos <- c("log_germ_tc",       
                   "log_toxin",    
                   "log_cfe", 
                   "log_growth",
                   "log_germ_tc_and_gly",
                   "log_sporulation")
  num_bin_pheno <- length(bin_phenos)
  num_cont_pheno <- length(cont_phenos)
  
  for (j in 1:num_cont_pheno) {
     phenos <- cont_phenos
     # continuous  - NO KEY
      temp_test <- "both"
      temp_outdir <- paste0(data_path, "/5_hogwash/", phenos[j], "/")
      temp_tree <- paste0(data_path, "/3_trees/", phenos[j], ".tree")
      temp_geno <- paste0(data_path, "/1_combine_genotype_matrices/combinedmat/cleaned_", phenos[j], "_pan_and_rereferenced_snp_indel.tsv")
      temp_pheno <- paste0(data_path, "/4_phenotypes/cleaned_", phenos[j], ".tsv")
      temp_name <- paste0("continuous_", phenos[j])
      temp_group_method <- "post-ar"
      temp_fdr <- indiv_fdr
      temp_key <- NULL
      temp_partition <- "standard" # "largemem"
      temp_memory <-  binary_memory # continuous_memory
      temp_time <- continuous_time
      temp_strain <- paste0(data_path, "/6_ribotype/", phenos[j], "_ribotype_no_colors.tsv")
      command <- write_command(temp_pheno, temp_geno, temp_tree, temp_name,
                               temp_outdir, perm, temp_fdr, bootstrap, temp_test,
                               temp_group_method, temp_key, temp_strain)
      make_sbatch(temp_name, temp_partition, temp_memory, temp_time, command)

      # continuous  - GROUPED: WITH KEY & post-AR
      temp_key <- paste0(gene_key_path, "/cleaned_", phenos[j], "_gene_and_intergenic_key.tsv")
      temp_name <- paste0("continuous_", phenos[j], "_grouped_", temp_group_method)
      temp_fdr <- group_fdr
      temp_memory <- group_memory
      temp_partition <- "standard"
      command <- write_command(temp_pheno, temp_geno, temp_tree, temp_name,
                               temp_outdir, perm, temp_fdr, bootstrap, temp_test,
                               temp_group_method, temp_key, temp_strain)
      make_sbatch(temp_name, temp_partition, temp_memory, temp_time, command)
  }
  
  for (j in 1:num_bin_pheno) {
      phenos <- bin_phenos
      # PHYC - Not grouped
      temp_test <- "phyc"
      temp_outdir <- paste0(data_path, "/5_hogwash/", phenos[j], "/")
      temp_tree <- paste0(data_path, "/3_trees/", phenos[j], ".tree")
      temp_geno <- paste0(data_path, "/1_combine_genotype_matrices/combinedmat/cleaned_", phenos[j], "_pan_and_rereferenced_snp_indel.tsv")
      temp_pheno <- paste0(data_path, "/4_phenotypes/cleaned_", phenos[j], ".tsv")
      temp_name <- paste0("phyc_", phenos[j])
      temp_group_method <- "post-ar"
      temp_memory <- binary_memory
      temp_time <- binary_time
      temp_partition <- "standard"
      temp_key <- NULL
      temp_fdr <- indiv_fdr
      temp_strain <- paste0(data_path, "/6_ribotype/", phenos[j], "_ribotype_no_colors.tsv")
      command <- write_command(temp_pheno, temp_geno, temp_tree, temp_name,
                               temp_outdir, perm, temp_fdr, bootstrap, temp_test,
                               temp_group_method, temp_key, temp_strain)
      make_sbatch(temp_name, temp_partition, temp_memory, temp_time, command)
      
      # PHYC - Grouped WITH KEY & post-AR
      temp_key <- paste0(gene_key_path, "/cleaned_", phenos[j], "_gene_and_intergenic_key.tsv")
      temp_name <- paste0("phyc_", phenos[j], "_grouped_", temp_group_method)
      temp_memory <- group_memory
      temp_fdr <- group_fdr
      command <- write_command(temp_pheno, temp_geno, temp_tree, temp_name,
                               temp_outdir, perm, temp_fdr, bootstrap, temp_test,
                               temp_group_method, temp_key, temp_strain)
      make_sbatch(temp_name, temp_partition, temp_memory, temp_time, command)

      # synchronous - NO KEY
      temp_key <- NULL
      temp_test <- "synchronous"
      temp_name <- paste0("synchronous_", phenos[j])
      temp_memory <- binary_memory
      temp_fdr <- indiv_fdr
      command <- write_command(temp_pheno, temp_geno, temp_tree, temp_name,
                               temp_outdir, perm, temp_fdr, bootstrap, temp_test,
                               temp_group_method, temp_key, temp_strain)
      make_sbatch(temp_name, temp_partition, temp_memory, temp_time, command)
      
      # synchronous - WITH KEY & post-ar
      temp_key <- paste0(gene_key_path, "/cleaned_", phenos[j], "_gene_and_intergenic_key.tsv")
      temp_name <- paste0("synchronous_", phenos[j], "_grouped_", temp_group_method)
      temp_memory <- group_memory
      temp_fdr <- group_fdr
      command <- write_command(temp_pheno, temp_geno, temp_tree, temp_name,
                               temp_outdir, perm, temp_fdr, bootstrap, temp_test,
                               temp_group_method, temp_key, temp_strain)
      make_sbatch(temp_name, temp_partition, temp_memory, temp_time, command)
  }
}

write_hogwash_sbat()

