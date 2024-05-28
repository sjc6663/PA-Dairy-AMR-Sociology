# Figure 7 - clinically significant genes, gender (rel abund, alpha, beta)
# SAB 09/12/2023

# null statistical hypothesis:
# alpha diversity - Shannon's Diversity Index of females = Shannon's Diversity Index of males
# beta diversity - sample distribution/makeup of females = sample distribution/makeup of males

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

ps2 <- readRDS("bovine-host-resistome/sig-decontam-ps.RDS")

# filter out taxa that are not clinically significant and remove them
sig <- subset_taxa(ps2, 
                   sig == "TRUE")

sig <- subset_samples(sig, 
                   Conventional.Organic == "Conventional")

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
phy <- dplyr::select(phy, c("Sample", "Broadclass", "Abundance"))

# get abundance as a percent and round to whole numbers
phy$percent <- phy$Abundance * 100
phy$percent <- round(phy$percent, digits = 0)

phy <- dplyr::select(phy, "Sample", "Broadclass", "percent")

phy

# alpha diversity (B) ----
# create data frame with relevant metadata for comparison
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(sig, measures = "Shannon"),
  "MF" = phyloseq::sample_data(sig)$Male.Female,
  "batch" = phyloseq::sample_data(sig)$Batch,
  "run" = phyloseq::sample_data(sig)$Run,
  "aff" = phyloseq::sample_data(sig)$affiliation, 
  "type" = phyloseq::sample_data(sig)$Conventional.Organic
)

# # male female
model <- lm(Shannon ~ run + batch, data = adiv)
residuals <- resid(model)

res <- as.data.frame(residuals)
res <- rownames_to_column(res, "samp")

adiv <- rownames_to_column(adiv, "samp")
dat1 <- merge(adiv, res, by = "samp")

dat1 <- column_to_rownames(dat1, "samp")

model2 <- lm(residuals ~ MF, data = dat1)
summary(model2) # P = 0.00


# violin plot 
ps.meta <- meta(sig)
ps.meta$Shannon <- phyloseq::estimate_richness(sig, measures = "Shannon")

ps.meta$'' <- alpha(sig, index = 'shannon')

B <- ggviolin(ps.meta, x = "Male.Female", y = "Shannon$Shannon",
              add = "boxplot", fill = "Male.Female", palette = c("#367aa1", "#def4e5"), title = "B", ylab = "Shannon's Diversity Index", xlab = " ") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 15), axis.title = element_text(size = 15)) +
  scale_y_continuous(limits = c(1, 7))
B

# beta diversity (C) ----

ait <- sig %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "Male.Female") %>% bdisp_get() # p=0.00015

# test with PERMANOVA
mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "Male.Female + Run + Batch",
    n_perms = 9999
  )

mod1 # R2 = 0.01, F(1, 65) = 1.06, P = 0.34
# plot

C <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Male.Female", size = 6, axes = c(2,3)) +
  scale_color_manual(values = c("#def4e5", "#367aa1")) +
  stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  labs(color = "Gender") +
  ggtitle("C") + 
  labs(caption = "R2 = 0.01, F(1, 65) = 1.06, P = 0.34") +
  theme(text = element_text(size = 15)) 
C

(A|B)/C

ggsave(filename = "plots/paper/figure7.pdf", dpi = 600, width = 12, height = 10)
