# Set up color palettes: 
cbf_1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

big_palette <- c('#e6194B',  '#3cb44b', '#a9a9a9','#ffe119', '#4363d8', '#f58231',
                 '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', 
                 '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075')

# skip the lightest yellow one: 13
big_palette_rainbow_order <- c(big_palette[c(14, 1, 6, 4, 15, 2, 10, 7, 5, 16, 11, 9, 8, 12, 3)], "snow2")


indiv_color_palette <- cbf_1[c(3, 4, 8)] # blue, green, pink
names(indiv_color_palette) <- c("Accessory", "Indel", "SNP")

group_color_palette <- cbf_1[c(2, 6)] # orange and navy
names(group_color_palette) <- c("Gene", "Intergenic region")

# group_big_color_pal <- c(big_palette[c(1, 2, 8, 4, 5, 6, 7, 9)], "snow2")
  
gene_and_intergenic_color_pal <- c("#e6194B", # 1 red
                                   "#800000", # 2 maroon
                                   "#f58231", # 3 orange
                                   "#9A6324", # 4 brown
                                   "#ffe119", # 5 yellow
                                   "#808000", # 6 olive
                                   "#3cb44b", # 7 green
                                   "#aaffc3", # 8 mint
                                   "#42d4f4", # 9 cyan
                                   "#469990", # 10 teal
                                   "#4363d8", # 11 blue
                                   "#000075", # 12 navy
                                   "#911eb4", # 13 purple, 
                                   "#dcbeff", # 14 lavender
                                   "#f032e6", # 15 magenta
                                   "snow2", # 16 grey
                                   "snow2") # 17 grey again
                                   
  #https://sashamaps.net/docs/resources/20-colors/
                       
names(gene_and_intergenic_color_pal) <- 
  c("cdu1",
    "cdu1-tcdR",
    "tcdR", 
    "tcdR-tcdB",
    "tcdB", 
    "tcdB-tcdE",
    "tcdE", 
    "tcdE-uviB pseudogene", 
    "uviB pseudogene", 
    "uviB pseudogene-tcdA",
    "tcdA", 
    "tcdA-tcdC", 
    "tcdC", 
    "tcdC-cdd1", 
    "cdd1", 
    "other", 
    "intergenic")

just_gene_palette <- gene_and_intergenic_color_pal[!grepl("-", names(gene_and_intergenic_color_pal)) & !grepl("other", names(gene_and_intergenic_color_pal))]

toxin_domain_palette <- c("#fabed4", # 16 pink 
                          "#ffd8b1", # apricot
                          "#aaffc3", # 8 mint
                          "#dcbeff", # 14 lavender
                          "snow2") # 5 grey
names(toxin_domain_palette) <- c("Glucosyltransferase", 
                                 "Autoprotease",
                                 "Pore-forming", 
                                 "CROPs", 
                                 "NA_domain")

promoter_palette <- c("#000075", # 1 navy
                      # "#469990", # 2 teal 
                      "#9A6324", # 3 brown
                      "snow2") # 4 grey
names(promoter_palette) <- c("Promoter", 
                             # "Overlapping promoters", 
                             "sigD promoter", 
                             "Not promoter")
