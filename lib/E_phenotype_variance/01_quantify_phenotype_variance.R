library(tidyverse)
options(scipen = 999) # turn off scientific notation

# For log normal data you can calculate the geometic coefficient of variance
get_gcv <- function(ln_sd){
  # https://en.wikipedia.org/wiki/Coefficient_of_variation
  sqrt(exp(ln_sd^2) - 1) # square root (e^((standard deviation of the data after a natural log transformation)^2) - 1)
}

phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")
phenotypes <- phenotypes %>% filter(bin_cont == "continuous")
stats_df <- as.data.frame(matrix(NA, nrow = 6, ncol = 9))
colnames(stats_df) <- c("phenotype", "SD", "mean-SD", "mean", "mean+SD", "data_outside_of_1SD", "variance", "spread", "gcv")
for (i in 1:nrow(phenotypes)) {
  current_pheno <- read.table(paste0("../../data/4_phenotypes/cleaned_", 
                                     phenotypes$shorthand[i], 
                                     ".tsv"))
  current_pheno <- as.matrix(current_pheno)
  
  stats_df$phenotype[i] <- phenotypes$shorthand[i]
  sd_pheno <- sd(current_pheno)
  stats_df$SD[i] <- sd_pheno
  mean_pheno <- summary(current_pheno)[4] %>% gsub(".*:", "", .) %>% gsub(" ", "", .) %>% as.numeric()
  stats_df$mean[i] <- mean_pheno
  stats_df$`mean-SD`[i] <- mean_pheno - sd_pheno
  stats_df$`mean+SD`[i] <- mean_pheno + sd_pheno
  num_outside <- sum(current_pheno < mean_pheno - sd_pheno | 
                       current_pheno > mean_pheno + sd_pheno) 
  stats_df$data_outside_of_1SD[i] <- 100 * num_outside / nrow(current_pheno)
  stats_df$variance[i] <- var(current_pheno)
  stats_df$spread[i] <- max(current_pheno) - min(current_pheno)
  stats_df$gcv[i] <- get_gcv(sd_pheno)
}

write_tsv(stats_df,  "../../data/4_phenotypes/phenotypic_analysis_results/geometric_coefficient_of_variance.tsv")


