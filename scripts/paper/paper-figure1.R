# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(viridis)
library(BiocManager)
library(stringr)
library(microshades)
library(microbiome)
library(patchwork)
library(ggpubr)

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

### Relative Abundance by Percent (Normalized Data) ----
library(microbiome)
library(dplyr)

psbclass <- aggregate_taxa(psrel, level = "Broadclass")
# fix taxa for aesthetics
# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ")

# subset phyloseq object by calves/cows
calves <- subset_samples(psbclass,
                         Group == "Calves")
cows <- subset_samples(psbclass, 
                       Group == "Cows")


A <- calves %>% plot_composition(average_by = "Farm", sample.sort = "Farm", x.label = " ") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position = "top") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  ggtitle("A")

B <- cows %>% plot_composition(average_by = "Farm", sample.sort = "Farm", x.label = " ") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position = "top") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  ggtitle("B")

C <- calves %>% plot_composition(average_by = "Farm", sample.sort = "Conventional.Organic", x.label = "Farm", group_by = "Conventional.Organic") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position = "top") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  ggtitle("C")

D <- cows %>% plot_composition(average_by = "Farm", sample.sort = "Conventional.Organic", x.label = "Farm", group_by = "Conventional.Organic") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position = "top") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  ggtitle("D")

(A|B)/(C|D)

all <- ggarrange(A, B, C, D,
                 ncol = 2, nrow = 2, 
                 common.legend = TRUE, legend = "bottom")

all

ggsave(filename = "plots/PCA-all-samples.pdf", dpi = 600, width = 12, height = 16)

