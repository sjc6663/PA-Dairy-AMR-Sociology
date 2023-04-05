# BARPLOTS FOR DESCRIPTIVE STATS
# SAB 4-4-2023

# color scheme: Viridis mako / microshades micro_cvd_blue
color_palette <- c("#BCE1FF", "#7DCCFF", "#56B4E9", "#098BD9", "#4292C6")

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(viridis)
library(BiocManager)
library(stringr)

# read in phyloseq object
ps <- readRDS("data/ransom/decontam-ps.RDS")

dist_mat <- phyloseq::distance(ps, method = "euclidean")

meltdf <- psmelt(ps)

write.csv(meltdf, file = "tables/melted-phyloseq.csv")

#plot 
ps %>% 
 # subset_taxa(
   # Broadclass == "Drugs") %>% 
  psmelt() %>% 
ggplot(aes(x=Broadclass, fill = Male.Female)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette) +
  facet_wrap(~Group, scales = "free",nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = "Male or Female") +
  theme(axis.text= element_text(size = 12,face = "bold"),
        axis.title= element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold", hjust = 0),
        legend.position = "right",
        strip.background = element_rect(fill = color_palette[2]),
        strip.text.x = element_text(size = 12, color = "black",face = "bold")
  ) +
  #guides(fill="none")+
  #scale_y_continuous(expand = c(0,0.1), limits = c(0,14))+
  coord_flip()

ggsave(filename = "plots/barplots/MF-calves-cows.pdf", dpi = 600)








# save work
save.image("data/beta-div.RData")