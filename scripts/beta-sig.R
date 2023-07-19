# Beta diversity of clinically significant genes
# SAB 07/17/2023

library(phyloseq)
library(microViz)
library(vegan)
library(plyr)
library(stringr)
library(dplyr)

color_palette <- c("#367aa1", "#348fa7", "#40b7ad", "#8ad9b1", "#def4e5")
small_color <- c("#367aa1", "#8ad9b1")

set.seed(81299)

# read in phyloseq object
ps <- readRDS("data/full-run/sig-decontam-ps.RDS")

ps2 <- rarefy_even_depth(ps)

# filter out taxa that are not clinically significant and remove them
ps <- subset_taxa(ps, 
                  sig == "TRUE")

ps

ps <- ps %>% 
  ps_mutate(
    employees = if_else(str_detect(Non.Family.Milkers, "0"), true = "No", false = "Yes")
  ) 

ps <- ps %>% 
  ps_mutate(
    meetings = if_else(str_detect(Formal.Team.Meetings.Frequency, "Never"), true = "No", false = "Yes")
  ) 
# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")


# clr transform phyloseq objects
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()

dist_mat <- phyloseq::distance(transps, method = "euclidean")

vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female) # p = 0.004**, SIGNIFICANT
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Group) # p = 0.001***, SIGNIFICANT
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female*phyloseq::sample_data(transps)$Group) # p = 0.454, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Formal.Team.Meetings.Frequency) # p = 0.934, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Cultural.Language.Barriers) # p = 0.988, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Herd.Size) # p = 0.958, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$employees) # p = 0.591, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$meetings) # p = 0.722, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Conventional.Organic) # p = 0.454, ns


# PCA Plots ---------------------------------------------------------------------------------------------------------

A <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Male.Female") +
  scale_color_manual(values = small_color) +
  stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  ggtitle("A") + 
  labs(caption = "R2 = 0.019, F(1,70) = 1.34, P = 0.004**")
A

B <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Group") +
  scale_color_manual(values = small_color) +
  stat_ellipse(aes(group = Group, color = Group)) + 
  theme_classic() +
  ggtitle("B") + 
  labs(caption = "R2 = 0.026, F(1,70) = 1.88, P = 0.001**")
B

A|B

ggsave(filename = "plots/presentation/PCA-crg-beta.pdf", dpi = 600, width = 10, height = 8)
