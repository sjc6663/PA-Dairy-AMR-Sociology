# Paper Figure 6 - Language Barrier: Relative Abund, Alpha, Beta
# SAB 11/16/2023

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
library(tibble)

set.seed(81299)

# read in phyloseq object
ps <- readRDS("bovine-host-resistome/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# relative abundance ----
psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# fix taxa for aesthetics
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

A <- psbclass %>% plot_composition(average_by = "Cultural.Language.Barriers", sample.sort = "Cultural.Language.Barriers", x.label = "Cultural.Language.Barriers") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position = "top") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  ggtitle("A")

# relative abundance plot percentages ----
# merge samples by comparative categories (metadata) we want
mps <- merge_samples(psrel, "Cultural.Language.Barriers")

# create a dataframe with the aggregated abundances for desired tax rank
phy <- mps %>% tax_glom(taxrank = "Broadclass") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

# select only relevant columns
phy <- select(phy, c("Sample", "Broadclass", "Abundance"))

# get abundance as a percent and round to whole numbers
phy$percent <- phy$Abundance * 100
phy$percent <- round(phy$percent, digits = 0)

phy <- select(phy, "Sample", "Broadclass", "percent")

phy

# alpha diversity ----
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "lang" = phyloseq::sample_data(ps)$Cultural.Language.Barriers,
  "batch" = phyloseq::sample_data(ps)$Batch,
  "run" = phyloseq::sample_data(ps)$Run,
  "aff" = phyloseq::sample_data(ps)$affiliation,
  "type" = phyloseq::sample_data(ps)$Conventional.Organic
)

# statistical test
model <- lm(Shannon ~ run + batch, data = adiv)
residuals <- resid(model)

res <- as.data.frame(residuals)
res <- rownames_to_column(res, "samp")

adiv <- rownames_to_column(adiv, "samp")
dat1 <- merge(adiv, res, by = "samp")

dat1 <- column_to_rownames(dat1, "samp")

model2 <- lm(residuals ~ lang, data = dat1)
summary(model2) # P = 0.94, ns

# violin plot
ps.meta <- meta(ps)
ps.meta$Shannon <- phyloseq::estimate_richness(ps, measures = "Shannon")

ps.meta$'' <- alpha(ps, index = 'shannon')

B <- ggviolin(ps.meta, x = "Cultural.Language.Barriers", y = "Shannon$Shannon",
              add = "boxplot", fill = "Cultural.Language.Barriers", palette = c("#38aaac", "#40498d"), title = "B", ylab = "Shannon's Diversity Index", xlab = " ") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(1, 7))

# beta diversity ----

ait <- ps %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "Cultural.Language.Barriers") %>% bdisp_get() # p=0.31

# test with PERMANOVA
mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "Cultural.Language.Barriers + Run + Batch",
    n_perms = 9999
  )

mod1 # R2 = 0.04, F(1, 24) = 0.94, P = 0.52

C <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Cultural.Language.Barriers", size = 6, axes = c(1,2)) +
  scale_color_manual(values = c("#40498d", "#38aaac")) +
  stat_ellipse(aes(group = Cultural.Language.Barriers, color = Cultural.Language.Barriers)) + 
  theme_classic() +
  labs(color = "Language Barriers") +
  ggtitle("C") + 
  labs(caption = "R2 = 0.01, F(1, 65) = 0.91, P = 0.61") +
  theme(text = element_text(size = 20)) 

(A|B)/C

A|B

ggsave(filename = "plots/paper/figure6.pdf", dpi = 600, width = 20, height = 14)
