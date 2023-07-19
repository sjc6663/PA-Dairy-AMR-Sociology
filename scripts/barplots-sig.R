# Barplot of Clinically Significant Taxa
# SAB 07/17/2023
# for colors use viridis mako 

color_palette <- c("#0070FF", "#D75CE0", "#FFC55A", "#FF8C76", "#F9F871", "#FF5EAA")

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(viridis)
library(BiocManager)
library(stringr)
library(microshades)
library(microbiome)

# read in phyloseq object
ps <- readRDS("data/full-run/sig-decontam-ps.RDS")

# filter out taxa that are not clinically significant and remove them
ps <- subset_taxa(ps, 
                  sig == "TRUE")

ps

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# melt phyloseq to work with it in barplots
meltdf <- psmelt(ps)
meltdf2 <- psmelt(psrel)

#filter out rows with 0 relabun 
out2 <- meltdf[meltdf$Abundance != 0.0000000000, ]
out3 <- meltdf2[meltdf2$Abundance != 0.0000000000, ]

# BARPLOTS ---------------------------------------------------------------------------------------------------------

# Male Female - Cows and Calves ----
ggplot(out3, aes(x=Broadclass, fill = Group)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette) +
  facet_wrap(~Male.Female, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = "Age Group",
       title = "C") +
  theme(axis.text= element_text(size = 12,face = "bold"),
        axis.title= element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold", hjust = 0),
        legend.position = "right",
        strip.background = element_rect(fill = color_palette[2]),
        strip.text.x = element_text(size = 12, color = "black",face = "bold")
  ) +
  #guides(fill="none")+
  scale_y_continuous(expand = c(0,0.1), limits = c(0,10000))+
  coord_flip()

ggsave(filename = "plots/full-run/MF-calves-cows.pdf", dpi = 600)


# RELATIVE ABUNDANCE / NORMALIZED PLOTS -----------------------------------------------------------------------------

psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

psbclass %>% plot_composition(group_by = "Farm", x.label = "Group") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "G", discrete = TRUE) +
  ggtitle("")

ggsave(filename = "plots/presentation/relabund-sig-bclass-farm-age.pdf", dpi = 600, width = 14, height = 10)

psbclass %>% plot_composition(group_by = "Male.Female", sample.sort = "Group", x.label = "Group") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("")

ggsave(filename = "plots/presentation/relabund-sig-bclass-gender-age.pdf", dpi = 600, width = 14, height = 10)

psbclass %>% plot_composition(average_by = "Male.Female", sample.sort = "Group", x.label = "Group") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("")

ggsave(filename = "plots/presentation/relabund-sig-avg-gender-age.pdf", dpi = 600, width = 14, height = 10)

psbclass %>% plot_composition(average_by = "Farm", sample.sort = "Group", x.label = "Group") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("")

ggsave(filename = "plots/presentation/relabund-sig-avg-farm-age.pdf", dpi = 600, width = 14, height = 10)

psbclass %>% plot_composition(group_by = "Housing.Style", sample.sort = "Group", x.label = "Group") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "G", discrete = TRUE) +
  ggtitle("")

ggsave(filename = "plots/presentation/relabund-sig-bclass-gender-age.pdf", dpi = 600, width = 14, height = 10)
