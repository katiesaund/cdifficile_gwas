library(tidyverse)
library(magrittr)

eggnog_path <- "../../data/18_eggnog/eggnog_run_2_vs_all/"
s6_tbl <- read_tsv(file = paste0(eggnog_path, "file_s6_pangenome_eggnog_key.tsv"))

num_na <- s6_tbl %>% pull(Eggnog_preferred_name) %>% is.na(.) %>% sum(.)
num_dash <- s6_tbl %>% filter(Eggnog_preferred_name == "-") %>% select(Eggnog_preferred_name) %>% nrow()
num_not_na_dash <- s6_tbl %>% filter(Eggnog_preferred_name != "-") %>% filter(!is.na(Eggnog_preferred_name)) %>%  nrow()
stopifnot(identical(num_na + num_dash + num_not_na_dash, nrow(s6_tbl)))

num_group <- sum(grepl("group", s6_tbl$Roary_ID))
num_not_group <- sum(!grepl("group", s6_tbl$Roary_ID))
stopifnot(identical(num_group + num_not_group, nrow(s6_tbl)))

# How many gene names can eggnog assign that roary didn't already
# This is rows where eggnog !- "-" or na and roary id grepl("group")
num_eggnog_assigns_but_roary_didnt <- 
  s6_tbl %>% 
  filter(Eggnog_preferred_name != "-") %>% 
  filter(!is.na(Eggnog_preferred_name)) %>% 
  filter(grepl("group", Roary_ID)) %>% 
  nrow()

num_neither_eggnog_nor_roary_assigns <- 
  s6_tbl %>% 
  filter(Eggnog_preferred_name == "-" | is.na(Eggnog_preferred_name)) %>% 
  filter(grepl("group", Roary_ID)) %>% 
  nrow()

# top 8
table_3_tbl <- read_tsv(paste0(eggnog_path, "table_3_cytotoxicity_with_eggnog_key.tsv"))

# Summary
num_not_group # Roary assigned it a gene name 
# 2608 
num_eggnog_assigns_but_roary_didnt # 
# 401

num_neither_eggnog_nor_roary_assigns # group by roary and dash or na by eggnog
# 7160 

stopifnot(identical(nrow(s6_tbl), 
          num_not_group + 
            num_eggnog_assigns_but_roary_didnt + 
            num_neither_eggnog_nor_roary_assigns))

# How many at least get a seed ortholog? 
sum(!is.na(s6_tbl$Eggnog_seed_ortholog)) # 7444
sum(is.na(s6_tbl$Eggnog_seed_ortholog)) # 2725

# How many got a seed ortholog but not a preferred name
num_eggnog_assigns_ortho_not_name <- 
  s6_tbl %>% 
  filter(Eggnog_preferred_name == "-" | is.na(Eggnog_preferred_name)) %>% 
  filter(!is.na(Eggnog_seed_ortholog))  %>% 
  nrow()
num_eggnog_assigns_ortho_not_name # 5790

num_eggnog_assigns_ortho_not_name_and_no_roary_id <- 
  s6_tbl %>% 
  filter(Eggnog_preferred_name == "-" | is.na(Eggnog_preferred_name)) %>% 
  filter(!is.na(Eggnog_seed_ortholog))  %>% 
  filter(grepl("group", Roary_ID)) %>% 
  nrow()
num_eggnog_assigns_ortho_not_name_and_no_roary_id # 4470


num_neither_eggnog_nor_roary_assigns_but_has_ortho <- 
  s6_tbl %>% 
  filter(Eggnog_preferred_name == "-" | is.na(Eggnog_preferred_name)) %>% 
  filter(grepl("group", Roary_ID)) %>% 
  filter(!is.na(Eggnog_seed_ortholog))  %>% 
  nrow()
