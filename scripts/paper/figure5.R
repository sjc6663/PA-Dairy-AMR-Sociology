# Paper Figure 5 - Employees: Relative Abund, Alpha, Beta
# SAB 09/13/2023

# setup ----
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
library(knitr)
library(dplyr)
library(car)
library(misty)
library(vegan)

set.seed(81299)

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

ps2 <- rarefy_even_depth(ps)

# fix metadata for aesthetics
sample_data(ps2)$employees <- gsub(sample_data(ps2)$employees, pattern = "No", replacement = "Family")
sample_data(ps2)$employees <- gsub(sample_data(ps2)$employees, pattern = "Yes", replacement = "Non-Family")

# transform to relative abundance
psrel <- microbiome::transform(ps2, "compositional")

# relative abundance ----
psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# fix taxa for aesthetics
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

A <- psbclass %>% plot_composition(average_by = "employees", sample.sort = "employees", x.label = "employees") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position = "top") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  ggtitle("A")

# alpha diversity ----
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps2, measures = "Shannon"),
  "employees" = phyloseq::sample_data(ps2)$employees,
  "batch" = phyloseq::sample_data(ps2)$Batch,
  "run" = phyloseq::sample_data(ps2)$Run,
  "aff" = phyloseq::sample_data(ps2)$affiliation
)

# test variance
varMF <- var.test(Shannon ~ employees, data = adiv, 
                  alternative = "two.sided")
varMF # p = 0.83, ns

# statistical test
model <- lm(Shannon ~ employees + run + batch + aff + employees*run*batch*aff, data = adiv)
summary(model) # p = 0.51, ns

# violin plot
ps.meta <- meta(ps2)
ps.meta$Shannon <- phyloseq::estimate_richness(ps2, measures = "Shannon")

ps.meta$'' <- alpha(ps2, index = 'shannon')

B <- ggviolin(ps.meta, x = "employees", y = "Shannon$Shannon",
              add = "boxplot", fill = "employees", palette = c("#38aaac", "#40498d"), title = "B", ylab = "Shannon's Diversity Index", xlab = " ") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(1, 6))

# beta diversity ----

# clr transform phyloseq objects
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()

dist_mat <- phyloseq::distance(transps, method = "euclidean")

vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$employees*phyloseq::sample_data(transps)$Run*phyloseq::sample_data(transps)$Batch*phyloseq::sample_data(transps)$affiliation)
# R2 = 0.02, F(1, 57) = 1.27, P = 0.036*

C <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "employees", size = 6, axes = c(2,3)) +
  scale_color_manual(values = c("#40498d", "#38aaac")) +
  stat_ellipse(aes(group = employees, color = employees)) + 
  theme_classic() +
  labs(color = "Employee Status") +
  ggtitle("C") + 
  labs(caption = "R2 = 0.02, F(1, 57) = 1.27, P = 0.036*") +
  theme(text = element_text(size = 20)) 

(A|B)/C

ggsave(filename = "plots/paper/figure5.pdf", dpi = 600, width = 20, height = 18)
