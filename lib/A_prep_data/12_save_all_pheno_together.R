cfe <-  read.table(file = "../../data/4_phenotypes/cleaned_log_cfe.tsv",
                   sep = "\t", 
                   header = TRUE,
                   row.names = 1, 
                   stringsAsFactors = FALSE)

germ_tc_and_gly <-  read.table(file = "../../data/4_phenotypes/cleaned_log_germ_tc_and_gly.tsv",
                               sep = "\t", 
                               header = TRUE,
                               row.names = 1, 
                               stringsAsFactors = FALSE)

germ_tc <-  read.table(file = "../../data/4_phenotypes/cleaned_log_germ_tc.tsv",
                       sep = "\t", 
                       header = TRUE,
                       row.names = 1, 
                       stringsAsFactors = FALSE)

growth <-  read.table(file = "../../data/4_phenotypes/cleaned_log_growth.tsv",
                      sep = "\t", 
                      header = TRUE,
                      row.names = 1, 
                      stringsAsFactors = FALSE)

sporulation <-  read.table(file = "../../data/4_phenotypes/cleaned_log_sporulation.tsv",
                           sep = "\t", 
                           header = TRUE,
                           row.names = 1, 
                           stringsAsFactors = FALSE)

toxin <-  read.table(file = "../../data/4_phenotypes/cleaned_log_toxin.tsv",
                     sep = "\t", 
                     header = TRUE,
                     row.names = 1, 
                     stringsAsFactors = FALSE)
sev <- read.table("../../data/4_phenotypes/cleaned_severity.tsv", 
                     sep = "\t", 
                     header = TRUE,
                     row.names = 1, 
                     stringsAsFactors = FALSE)

cfe <- cbind(row.names(cfe), cfe)
germ_tc_and_gly <- cbind(row.names(germ_tc_and_gly), germ_tc_and_gly)
germ_tc <- cbind(row.names(germ_tc), germ_tc)
growth <- cbind(row.names(growth), growth)
sporulation <- cbind(row.names(sporulation), sporulation)
toxin <- cbind(row.names(toxin), toxin)
sev <- cbind(row.names(sev), sev)


colnames(cfe)[1] <- "ID"
colnames(germ_tc_and_gly)[1] <- "ID"
colnames(germ_tc)[1] <- "ID"
colnames(growth)[1] <- "ID"
colnames(sporulation)[1] <- "ID"
colnames(toxin)[1] <- "ID"
colnames(sev)[1] <- "ID"


pheno_df <- full_join(toxin, germ_tc_and_gly, by = "ID")
pheno_df <- full_join(pheno_df, germ_tc, by = "ID")
pheno_df <- full_join(pheno_df, growth, by = "ID")
pheno_df <- full_join(pheno_df, sporulation, by = "ID")
pheno_df <- full_join(pheno_df, cfe, by = "ID")
pheno_df <- full_join(pheno_df, sev, by = "ID")

phenotype <- pheno_df
phenotype$ID <- as.character(phenotype$ID)
phenotype$log_toxin <- as.numeric(as.character(phenotype$log_toxin))
phenotype$log_germ_tc_and_gly <- as.numeric(as.character(phenotype$log_germ_tc_and_gly))
phenotype$log_germ_tc <- as.numeric(as.character(phenotype$log_germ_tc))
phenotype$log_growth <- as.numeric(as.character(phenotype$log_growth))
phenotype$log_sporulation <- as.numeric(as.character(phenotype$log_sporulation))
phenotype$log_cfe <- as.numeric(as.character(phenotype$log_cfe))
phenotype$severity <- as.numeric(as.character(phenotype$severity))

write_tsv(phenotype, "../../data/4_phenotypes/all_phenos.tsv")
