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
library(microbiome)

# read in phyloseq object
ps <- readRDS("data/full-run/sig-decontam-ps.RDS")
ps2 <- readRDS("data/full-run/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")
psrel2 <- microbiome::transform(ps2, "compositional")

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
C <- ggplot(out3, aes(x=sig, fill = Group)) + 
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
  scale_y_continuous(expand = c(0,0.1), limits = c(0,40000))+
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
library(dplyr)

psclass <- aggregate_taxa(psrel, level = "sig")
psbclass <- aggregate_taxa(psrel2, level = "Broadclass")

psrel@tax_table@sig[psrel@tax_table@sig == "TRUE"] <- "sig"
tax_table(psrel)$sig[tax_table(psrel)$sig == "FALSE"] <- "ns"


pssig <- tax_glom(psrel, taxrank = "sig")

# fix taxa for aesthetics
# find and substitute
taxa_names(psclass) <- gsub(taxa_names(psclass), pattern = "_", replacement = " ") 

gender <- psclass %>% plot_composition(group_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
 # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  ggtitle("A")

ggsave(filename = "plots/full-run/relabund-averaged.pdf", dpi = 600, width = 24, height = 12)

age <- psbclass %>% plot_composition(group_by = "Group") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("A")

ggsave(filename = "plots/full-run/relabund-total-bclass-group.pdf", dpi = 600, width = 24, height = 12)

calves <- subset_samples(
  psbclass,
  Group == "Calves"
)

cows <- subset_samples(
  psbclass,
  Group == "Cows"
)

A <- calves %>% plot_composition(group_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
  theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("B Calves")

B <- cows %>% plot_composition(group_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("Cows")

MFAG <- A|B

ggsave(filename = "plots/full-run/relabund-bclass-all-samples-group-gender.pdf", dpi = 600, width = 18, height = 12)

sample_data(psbclass)$Herd.Size <- factor(sample_data(psbclass)$Herd.Size,
                         levels = c("0-50", "50-100", "100-150", "150-200", "200+"))

size <- psbclass %>% plot_composition(group_by = "Herd.Size") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("C")

ggsave(filename = "plots/full-run/relabund-all-samples-bclass-herd-size.pdf", dpi = 600, width = 20, height = 12)


psbclass %>% plot_composition(group_by = "Group") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("Age Group")

ggsave(filename = "plots/full-run/relabund-all-samples-group.pdf", dpi = 600, width = 18, height = 12)

type <- psbclass %>% plot_composition(group_by = "Conventional.Organic") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("B")

ggsave(filename = "plots/full-run/relabund-all-samples-bclass-farm-type.pdf", dpi = 600, width = 20, height = 12)

language <- psbclass %>% plot_composition(group_by = "Cultural.Language.Barriers") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("D")

ggsave(filename = "plots/full-run/relabund-all-samples-bclass-langbarrier.pdf", dpi = 600, width = 20, height = 12)


sample_data(psbclass)$Formal.Team.Meetings.Frequency <- factor(sample_data(psbclass)$Formal.Team.Meetings.Frequency,
                                              levels = c("Never", "1 or 2 times/year", "Quarterly", "Once a month", "At least twice a month"))


meetings <- psbclass %>% plot_composition(group_by = "Formal.Team.Meetings.Frequency") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("Team Meeting Frequency") +
  ggtitle("C")

ggsave(filename = "plots/full-run/relabund-all-samples-bclass-teammeetings.pdf", dpi = 600, width = 20, height = 12)


age / type / size

ggsave(filename = "plots/full-run/relabund-farm-details.pdf", dpi = 600, width = 20, height = 24)


gender / MFAG / meetings / language

ggsave(filename = "plots/full-run/relabund-soc-details.pdf", dpi = 600, width = 20, height = 24)


## BROADCLASS GROUP 


psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

gender2 <- psbclass %>% plot_composition(average_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("D")

ggsave(filename = "plots/full-run/relabund-averaged-bclass-gender.pdf", dpi = 600, width = 12, height = 14)


size2 <- psbclass %>% plot_composition(average_by = "Herd.Size") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("C")

ggsave(filename = "plots/full-run/relabund-averaged-bclass-herd-size.pdf", dpi = 600, width = 12, height = 14)

group <- psbclass %>% plot_composition(average_by = "Group") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("A")

ggsave(filename = "plots/full-run/relabund-averaged-bclass-agegroup.pdf", dpi = 600, width = 12, height = 14)

farmT <- psbclass %>% plot_composition(average_by = "Conventional.Organic") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("B")

ggsave(filename = "plots/full-run/relabund-averaged-bclass-farmtype.pdf", dpi = 600, width = 12, height = 14)

sample_data(psbclass)$Formal.Team.Meetings.Frequency <- factor(sample_data(psbclass)$Formal.Team.Meetings.Frequency,
                                                              levels = c("Never", "1 or 2 times/year", "Quarterly", "Once a month", "At least twice a month"))

meetings2 <- psbclass %>% plot_composition(average_by = "Formal.Team.Meetings.Frequency") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("F")

ggsave(filename = "plots/full-run/relabund-averaged-bclass-teammeetings.pdf", dpi = 600, width = 12, height = 14)


language2 <- psbclass %>% plot_composition(average_by = "Cultural.Language.Barriers") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("G")

ggsave(filename = "plots/full-run/relabund-averaged-bclass-langbarrier.pdf", dpi = 600, width = 12, height = 14)

group / farmT / size2

ggsave(filename = "plots/full-run/averaged-relabund-farm-details.pdf", dpi = 600, width = 12, height = 16)

A2 <- calves %>% plot_composition(average_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
  theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("E Calves")

B2 <- cows %>% plot_composition(average_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("Cows")

MFAG2 <- A2|B2

ggsave(filename = "plots/full-run/averaged-relabund-cows-calves-gender.pdf", dpi = 600, width = 12, height = 14)

gender2 / MFAG2 / meetings2 / language2

ggsave(filename = "plots/full-run/averaged-relabund-soc-details.pdf", dpi = 600, width = 12, height = 14)

(group / farmT / size2)  | (gender2 / MFAG2 / meetings2 / language2)

ggsave(filename = "plots/full-run/averaged-all-details.pdf", dpi = 600, width = 24, height = 30)

# save work ----
save.image("data/beta-div.RData")

## MECHANISM
test <- psrel@tax_table %>% 
  mutate(beta = if_else(str_detect(Class, "betalactams"), true = "Yes", false = "No"))

psmclass <- aggregate_taxa(psrel, level = "Mechanism")

psgclass <- aggregate_taxa(psrel, level = "Gene")

# find and substitute
taxa_names(psmclass) <- gsub(taxa_names(psmclass), pattern = "_", replacement = " ") 

tax_table(psgclass)

psmclass %>% plot_composition(average_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
  theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE)

ggsave(filename = "plots/full-run/relabund-averaged-mechanism-gender.pdf", dpi = 600, width = 12, height = 14)

psbclass %>% plot_composition(average_by = "Farm", group_by = "Cultural.Language.Barriers") +
  # scale_y_continuous(labels = percent) +
  #theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE)

beta <- subset_taxa(
  psgclass,
  Class == "betalactams")

amino <- subset_taxa(
  psgclass,
  Class == "Aminoglycosides")

fluoro <- subset_taxa(
  psgclass,
  Class == "Fluoroquinolones")

fosfo <- subset_taxa(
  psgclass,
  Class == "Fosfomycin")

glyco <- subset_taxa(
  psgclass,
  Class == "Glycopeptides")

lipo <- subset_taxa(
  psgclass,
  Class == "Lipopeptides")

mls <- subset_taxa(
  psgclass,
  Class == "MLS")

rif <- subset_taxa(
  psgclass,
  Class == "Rifampin")

tax_table(beta)

beta %>% plot_composition(average_by = "Farm", group_by = "Conventional.Organic") +
   #scale_y_continuous(limits = c(0, 1)) +
  #theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE)

test <- mutate_tax_table(psmclass, dplyr::mutate(beta = if_else(str_detect(Class, "betalactams"), true = "Yes", false = "No")))

# clinically significant

psg <- aggregate_taxa(psrel, level = "sig", verbose = TRUE)

tax_table(ps) <- tax_table(ps)[,2:5]


psmclass %>% plot_composition(average_by = "Male.Female") +
  # scale_y_continuous(labels = percent) +
  theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE)

