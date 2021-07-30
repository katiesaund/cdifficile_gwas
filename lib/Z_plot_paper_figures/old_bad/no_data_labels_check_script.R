# Goal
# Make sure ND labels are being assigned to correct positions
# Pairs with fig1.R


cfe_nd
# DA00315 DA00484 VPI

sev_nd
# "Cd093" "Cd072" "Cd182" "Cd054" "VPI"  
temp_tr_plot <- ggtree(tree) %<+% ribotype + 
  geom_tippoint(aes(color = Ribotype), size = .9) + 
  geom_tiplab(aes(angle=90)) + 
  scale_color_manual(values = ribotype_colors) + 
  coord_flip() + 
  theme(legend.position = "none")

temp_tr_plot 

toxin_plot %>% 
  insert_top(cfe_plot) %>% 
  insert_bottom(temp_tr_plot, height = 4) %>% 
  insert_bottom(ggpubr::as_ggplot(ribo_legend))

toxin_plot %>% 
  insert_top(sporulation_plot) %>% 
  insert_bottom(temp_tr_plot, height = 4) %>% 
  insert_bottom(ggpubr::as_ggplot(ribo_legend))

