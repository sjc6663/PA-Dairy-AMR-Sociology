# Figure 7 - clinically significant genes, gender (rel abund, alpha, beta)
# SAB 09/12/2023

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

ps2 <- readRDS("data/full-run/sig-decontam-ps.RDS")

# filter out taxa that are not clinically significant and remove them
sig <- subset_taxa(ps2, 
                   sig == "TRUE")

# barplot (A) ----
psrel <- microbiome::transform(sig, "compositional")

psbclass <- aggregate_taxa(psrel, level = "Broadclass")
# fix taxa for aesthetics
# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

A <- psbclass %>% plot_composition(average_by = "Male.Female", sample.sort = "Male.Female", x.label = "Male.Female") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 15)) + 
  theme(legend.position = "top") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  ggtitle("A")
A

# relative abundance plot percentages ----
# merge samples by comparative categories (metadata) we want
mps <- merge_samples(psrel, "Male.Female")

# create a dataframe with the aggregated abundances for desired tax rank
phy <- mps %>% tax_glom(taxrank = "Broadclass") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

# select only relevant columns
phy <- select(phy, c("Sample", "Broadclass", "Abundance"))

# get abundance as a percent and round to whole numbers
phy$percent <- phy$Abundance * 100
phy$percent <- round(phy$percent, digits = 0)

phy <- select(phy, "Sample", "Broadclass", "percent")

phy

# alpha diversity (B) ----
# create data frame with relevant metadata for comparison
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(sig, measures = "Shannon"),
  "MF" = phyloseq::sample_data(sig)$Male.Female
)

# test variance
# Male Female 
varMF <- var.test(Shannon ~ MF, data = adiv, 
                  alternative = "two.sided")
varMF # p = 0.55, not sig

# # male female
wtestMF <- t.test(Shannon ~ MF, data = adiv, var.equal = FALSE)
wtestMF # p = 0.0005**, SIGNIFICANT

# violin plot 
ps.meta <- meta(sig)
ps.meta$Shannon <- phyloseq::estimate_richness(sig, measures = "Shannon")

ps.meta$'' <- alpha(sig, index = 'shannon')

B <- ggviolin(ps.meta, x = "Male.Female", y = "Shannon$Shannon",
              add = "boxplot", fill = "Group", palette = c("#367aa1", "#def4e5"), title = "B", ylab = "Shannon's Diversity Index", xlab = " ") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 15), axis.title = element_text(size = 15)) +
  scale_y_continuous(limits = c(1, 7))
B

# beta diversity (C) ----
# transform to relative abundance
psrel <- microbiome::transform(sig, "compositional")

# clr transform
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()

# create distance matrix 
dist_mat <- phyloseq::distance(transps, method = "euclidean")

# test
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female) # R2 = 0.019, F(1,70) = 1.34, P = 0.007**

# plot

C <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Male.Female", size = 6, axes = c(1,2)) +
  scale_color_manual(values = c("#367aa1", "#def4e5")) +
  stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  labs(color = "Gender") +
  ggtitle("C") + 
  labs(caption = "R2 = 0.019, F(1,70) = 1.34, P = 0.007**") +
  theme(text = element_text(size = 15)) 
C

(A|B)/C

ggsave(filename = "plots/paper/figure7.pdf", dpi = 600, width = 12, height = 10)
