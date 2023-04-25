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
library(microshades)

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

meltdf <- psmelt(ps)
meltdf2 <- psmelt(psrel)

#filter out rows with 0 relabun 
out2 <- meltdf[meltdf$Abundance != 0.0000000000, ]
out3 <- meltdf2[meltdf2$Abundance != 0.0000000000, ]

write.csv(meltdf, file = "tables/melted-phyloseq.csv")
write.csv(meltdf2, file = "tables/melted-phyloseq-relabund.csv")


### ------------------------------------------------------------------------
#plot 

# Male Female - Cows and Calves ----
C <- ggplot(out2, aes(x=Broadclass, fill = Group)) + 
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
  #scale_y_continuous(expand = c(0,0.1), limits = c(0,7500))+
  coord_flip()

ggsave(filename = "plots/full-run/MF-calves-cows.pdf", dpi = 600)

# Male Female ----
A <- ggplot(out2, aes(x=Broadclass, fill = Sample.Type)) + 
  geom_bar()+
  theme_bw() +
  scale_fill_manual(values = color_palette[3:5]) +
  facet_wrap(~Male.Female, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = " ",
       title = " ") +
  theme(axis.text= element_text(size = 12,face = "bold"),
        axis.title= element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold", hjust = 0),
        legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 12, color = "black",face = "bold")+
  ) +
  #guides(fill="none")+
  scale_y_continuous(expand = c(0,0.1), limits = c(0,20000))+
  coord_flip()
A
ggsave(filename = "plots/full-run/MF.pdf", dpi = 600)

A|B|C

ggsave(filename = "plots/full-run/MF-cow-calf-combo.pdf", dpi = 600, height = 12, width = 18)

# Cows Calves ----
B <- ggplot(out2, aes(x=Broadclass, fill = Sample.Type)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette[3:5]) +
  facet_wrap(~Group, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = " ",
       title = "A") +
  theme(axis.text= element_text(size = 12,face = "bold"),
        axis.title= element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold", hjust = 0),
        legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 12, color = "black",face = "bold")
  ) +
  #guides(fill="none")+
  #scale_y_continuous(expand = c(0,0.1), limits = c(0,7500))+
  coord_flip()

ggsave(filename = "plots/full-run/calves-cows.pdf", dpi = 600)

# Herd Size ----

out2$Herd.Size <- factor(out2$Herd.Size,
       levels = c("0-50", "50-100", "100-150", "150-200", "200+"))

HS <- ggplot(out2, aes(x=Broadclass, fill = Herd.Size)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette) +
  facet_wrap(~Group, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = "Herd Size",
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
  #scale_y_continuous(expand = c(0,0.1), limits = c(0,7500))+
  coord_flip()

ggsave(filename = "plots/full-run/herd-size.pdf", dpi = 600)

# Organic Conventional ----
OC <- ggplot(out2, aes(x=Broadclass, fill = Conventional.Organic)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette) +
  facet_wrap(~Group, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = "Farm Type",
       title = "B") +
  theme(axis.text= element_text(size = 12,face = "bold"),
        axis.title= element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold", hjust = 0),
        legend.position = "right",
        strip.background = element_rect(fill = color_palette[2]),
        strip.text.x = element_text(size = 12, color = "black",face = "bold")
  ) +
  #guides(fill="none")+
  #scale_y_continuous(expand = c(0,0.1), limits = c(0,7500))+
  coord_flip()

ggsave(filename = "plots/full-run/con-org.pdf", dpi = 600)
library(patchwork)
B | OC | HS

ggsave(filename = "plots/full-run/soc-barplot.pdf", dpi = 600, height = 12, width = 18)


# Team Meetings ----

out2$Formal.Team.Meetings.Frequency <- factor(out2$Formal.Team.Meetings.Frequency,
                                                           levels = c("Never", "1 or 2 times/year", "Quarterly", "Once a month", "At least twice a month"))


TM <- ggplot(out2, aes(x=Broadclass, fill = Formal.Team.Meetings.Frequency)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette) +
  facet_wrap(~Group, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = "Formal Team Meeting Frequency",
       title = "A") +
  theme(axis.text= element_text(size = 12,face = "bold"),
        axis.title= element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold", hjust = 0),
        legend.position = "right",
        strip.background = element_rect(fill = color_palette[2]),
        strip.text.x = element_text(size = 12, color = "black",face = "bold")
  ) +
  #guides(fill="none")+
  #scale_y_continuous(expand = c(0,0.1), limits = c(0,7500))+
  coord_flip()

ggsave(filename = "plots/full-run/team-meetings.pdf", dpi = 600)


# Language Barriers ----
LB <- ggplot(out2, aes(x=Broadclass, fill = Cultural.Language.Barriers)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette) +
  facet_wrap(~Group, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = "Language Barriers Present",
       title = "B") +
  theme(axis.text= element_text(size = 12,face = "bold"),
        axis.title= element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold", hjust = 0),
        legend.position = "right",
        strip.background = element_rect(fill = color_palette[2]),
        strip.text.x = element_text(size = 12, color = "black",face = "bold")
  ) +
  #guides(fill="none")+
  #scale_y_continuous(expand = c(0,0.1), limits = c(0,7500))+
  coord_flip()

ggsave(filename = "plots/full-run/lang-barriers.pdf", dpi = 600)

TM | LB

ggsave(filename = "plots/full-run/soc-barplots-2.pdf", dpi = 600, width = 18, height = 12)


### Relative Abundance by Percent (Normalized Data) ----
library(microbiome)

psclass <- aggregate_taxa(psrel, level = "Class")

# fix taxa for aesthetics
# find and substitute
taxa_names(psclass) <- gsub(taxa_names(psclass), pattern = "_", replacement = " ") 

psclass %>% plot_composition(average_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE)

ggsave(filename = "plots/full-run/relabund-averaged.pdf", dpi = 600, width = 24, height = 12)


psclass %>% plot_composition(group_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE)

ggsave(filename = "plots/full-run/relabund-all-samples.pdf", dpi = 600, width = 18, height = 12)


psclass %>% plot_composition(group_by = "Herd.Size") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("Herd Size")

ggsave(filename = "plots/full-run/relabund-all-samples-herd-size.pdf", dpi = 600, width = 20, height = 12)


psclass %>% plot_composition(group_by = "Group") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("Age Group")

ggsave(filename = "plots/full-run/relabund-all-samples-group.pdf", dpi = 600, width = 18, height = 12)


psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

psbclass %>% plot_composition(average_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("Primary Manager Gender")

ggsave(filename = "plots/full-run/relabund-averaged-bclass-gender.pdf", dpi = 600, width = 12, height = 14)


psbclass %>% plot_composition(average_by = "Herd.Size") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE)

ggsave(filename = "plots/full-run/relabund-averaged-bclass-herd-size.pdf", dpi = 600, width = 12, height = 14)

psbclass %>% plot_composition(average_by = "Group") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE)

ggsave(filename = "plots/full-run/relabund-averaged-bclass-agegroup.pdf", dpi = 600, width = 12, height = 14)



# save work ----
save.image("data/beta-div.RData")






