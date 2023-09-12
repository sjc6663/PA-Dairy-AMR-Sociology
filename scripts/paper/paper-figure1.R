# Figure 2 - Barplot of Relative Abundance (small multiples) Calves/Cows + Org/Conv
# SAB 09/12/2023

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
#install.packages("gridExtra")
library(gridExtra)

set.seed(81299)

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
  theme(legend.position = "none") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("A")

B <- cows %>% plot_composition(average_by = "Farm", sample.sort = "Farm", x.label = " ") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position = "none") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("B")

C <- calves %>% plot_composition(average_by = "Farm", sample.sort = "Conventional.Organic", x.label = "Farm", group_by = "Conventional.Organic") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position = "none") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("C")

D <- cows %>% plot_composition(average_by = "Farm", sample.sort = "Conventional.Organic", x.label = "Farm", group_by = "Conventional.Organic") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position = "none") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("D")

legend <- cows %>% plot_composition(average_by = "Farm", sample.sort = "Conventional.Organic", x.label = "Farm", group_by = "Conventional.Organic") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 40)) + 
  theme(legend.position = "top") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("D")

#(A|B)/(C|D)

#all <- ggarrange(A, B, C, D,
#                 ncol = 2, nrow = 2, 
#                 common.legend = TRUE, legend = "bottom")

#all <- A + B + C + D + plot_layout(ncol = 2, nrow = 2, 
#                                   guides = "collect") +
#  theme(legend.position = "top")

# all

# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legend
shared_legend <- extract_legend(legend)

# Draw plots with shared legend
all <- grid.arrange(arrangeGrob(A, B, C, D, ncol = 2, nrow = 2),
             shared_legend, nrow = 2, heights = c(10, 1))

ggsave(all, filename = "plots/paper/figure2.pdf", dpi = 600, width = 24, height = 20)

