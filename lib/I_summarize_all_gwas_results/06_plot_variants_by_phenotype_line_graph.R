# Plot the number of variants: 
#   - given to hogwash (left),
#   - tested (middle),
#   - and high interest (right
#   in a line graph, one line per phenotype
source("gwas_analysis_lib.R")
phenotypes <- read_tsv("../../data/4_phenotypes/pheno_names.tsv")

variant_indiv_df <- read_tsv("../../data/12_variant_types/variant_type_indiv.tsv")
variant_group_df <- read_tsv("../../data/12_variant_types/variant_type_group.tsv")

variant_df <- rbind(variant_indiv_df, variant_group_df)

variant_df$`Input or output`[variant_df$`Input or output` == "input"] <- "1. Fed to hogwash"
variant_df$`Input or output`[variant_df$`Input or output` == "output"] <- "2. Tested by hogwash"
variant_df$`Input or output`[variant_df$`Input or output` == "sig. output"] <- "3. Sig. P & Hi. Epsilon"

summary_indiv <- variant_df %>% 
  filter(`Indiv. or Group.` == "Individual") %>% 
  group_by(shorthand, `Input or output`, Test) %>% 
  summarize("Variants (#)" = sum(n)) %>% 
  mutate("Variant type" = `Input or output`) 

summary_indiv <- left_join(summary_indiv, phenotypes, by = "shorthand") %>% 
  mutate("Phenotype" = longhand)

summary_indiv %>% 
  filter(Test == "phyc") %>% 
  ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype, shape = Test)) + 
  geom_line() + 
  theme_bw() + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../../figures/gwas_results/variant_funnel_through_hogwash_indiv_phyc.jpg", units = "in", width = 6, height = 5)

summary_indiv %>% 
  filter(Test == "synchronous") %>% 
  ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype, shape = Test)) + 
  geom_line() + 
  theme_bw() + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../../figures/gwas_results/variant_funnel_through_hogwash_indiv_synchronous.jpg", units = "in", width = 6, height = 5)


summary_indiv %>% 
  filter(Test == "continuous") %>% 
  ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype, shape = Test)) + 
  geom_line() + 
  theme_bw() + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../../figures/gwas_results/variant_funnel_through_hogwash_indiv_continuous.jpg", units = "in", width = 6, height = 5)


summary_group <- variant_df %>% 
  filter(`Indiv. or Group.` == "Grouped") %>% 
  group_by(shorthand, `Input or output`, Test) %>% 
  summarize("Variants (#)" = sum(n)) %>% 
  mutate("Variant type" = `Input or output`) 

summary_group <- left_join(summary_group, phenotypes, by = "shorthand") %>% 
  mutate("Phenotype" = longhand)

summary_group %>% 
  filter(Test == "phyc") %>% 
  ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype)) + 
  geom_line() + 
  theme_bw() + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../../figures/gwas_results/variant_funnel_through_hogwash_group_phyc.jpg", units = "in", width = 6, height = 5)

summary_group %>% 
  filter(Test == "synchronous") %>% 
  ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype)) + 
  geom_line() + 
  theme_bw() + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../../figures/gwas_results/variant_funnel_through_hogwash_group_synchronous.jpg", units = "in", width = 6, height = 5)

summary_group %>% 
  filter(Test == "continuous") %>% 
  ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype)) + 
  geom_line() + 
  theme_bw() + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../../figures/gwas_results/variant_funnel_through_hogwash_group_continuous.jpg", units = "in", width = 6, height = 5)


for (i in 1:nrow(phenotypes)) {
  name <- phenotypes$shorthand[i]
  longhand <- phenotypes$longhand[i]
  
  if (phenotypes$bin_cont[i] == "continuous") {
    variant_df %>% 
      filter(shorthand == name, 
             `Indiv. or Group.` == "Individual") %>% 
      ggplot(aes(y = n, x = `Input or output`, group = `Variant type`, col = `Variant type`)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      ylab("Variant (#)") +
      xlab("") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_variant_funnel_through_hogwash_indiv_by_type.jpg"), 
           units = "in", width = 6, height = 5)  
    
    variant_df %>% 
      filter(shorthand == name, 
             `Indiv. or Group.` == "Grouped") %>% 
      ggplot(aes(y = n, x = `Input or output`, group = `Variant type`, col = `Variant type`)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      ylab("Variant (#)") +
      xlab("") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_variant_funnel_through_hogwash_group_by_type.jpg"), 
           units = "in", width = 6, height = 5)  
    
    summary_indiv %>% 
      filter(shorthand == name) %>% 
      ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_variant_funnel_through_hogwash_indiv.jpg"), 
           units = "in", width = 6, height = 5)
    
    summary_group %>% 
      filter(shorthand == name) %>% 
      ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_variant_funnel_through_hogwash_group.jpg"), 
           units = "in", width = 6, height = 5)
    
  } else {# binary have both sync and phyC 
    # PhyC  ---- 
    variant_df %>% 
      filter(shorthand == name, 
             `Indiv. or Group.` == "Individual", 
             Test == "phyc") %>% 
      ggplot(aes(y = n, x = `Input or output`, group = `Variant type`, col = `Variant type`)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      ylab("Variant (#)") +
      xlab("") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_phyc_variant_funnel_through_hogwash_indiv_by_type.jpg"), 
           units = "in", width = 6, height = 5)  
    
    variant_df %>% 
      filter(shorthand == name, 
             `Indiv. or Group.` == "Grouped", 
             Test == "phyc") %>% 
      ggplot(aes(y = n, x = `Input or output`, group = `Variant type`, col = `Variant type`)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      ylab("Variant (#)") +
      xlab("") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_phyc_variant_funnel_through_hogwash_group_by_type.jpg"), 
           units = "in", width = 6, height = 5)  
    
    summary_indiv %>% 
      filter(shorthand == name, 
             Test == "phyc") %>%
      ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_phyc_variant_funnel_through_hogwash_indiv.jpg"), 
           units = "in", width = 6, height = 5)
    
    summary_group %>% 
      filter(shorthand == name, 
             Test == "phyc") %>% 
      ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_phyc_variant_funnel_through_hogwash_group.jpg"), 
           units = "in", width = 6, height = 5)
    
    # synchronous  ---- 
    variant_df %>% 
      filter(shorthand == name, 
             `Indiv. or Group.` == "Individual", 
             Test == "synchronous") %>% 
      ggplot(aes(y = n, x = `Input or output`, group = `Variant type`, col = `Variant type`)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      ylab("Variant (#)") +
      xlab("") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_synchronous_variant_funnel_through_hogwash_indiv_by_type.jpg"), 
           units = "in", width = 6, height = 5)  
    
    variant_df %>% 
      filter(shorthand == name, 
             `Indiv. or Group.` == "Grouped", 
             Test == "synchronous") %>% 
      ggplot(aes(y = n, x = `Input or output`, group = `Variant type`, col = `Variant type`)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      ylab("Variant (#)") +
      xlab("") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_synchronous_variant_funnel_through_hogwash_group_by_type.jpg"), 
           units = "in", width = 6, height = 5)  
    
    summary_indiv %>% 
      filter(shorthand == name, 
             Test == "synchronous") %>%
      ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_synchronous_variant_funnel_through_hogwash_indiv.jpg"), 
           units = "in", width = 6, height = 5)
    
    summary_group %>% 
      filter(shorthand == name, 
             Test == "synchronous") %>% 
      ggplot(aes(y = `Variants (#)`, x = `Variant type`, group = Phenotype, col = Phenotype)) + 
      geom_line() + 
      theme_bw() + 
      geom_point() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("../../figures/gwas_results/",
                  name, 
                  "_synchronous_variant_funnel_through_hogwash_group.jpg"), 
           units = "in", width = 6, height = 5)
  }
}
