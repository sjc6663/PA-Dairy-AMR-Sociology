# Getting Percentages for Relative Abundance

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
library(tibble)
library(dplyr)

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# add column about non family employees for that comparison
psrel <- psrel %>% 
  ps_mutate(
    employees = if_else(str_detect(Non.Family.Milkers, "0"), true = "Family", false = "Non-Family")
  ) 

# merge samples by comparative categories (metadata) we want
mps <- merge_samples(psrel, "Cultural.Language.Barriers")

# create a dataframe with the aggregated abundances for desired tax rank
phy <- mps %>% tax_glom(taxrank = "Broadclass") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

# select only relevant columns
phy <- select(phy, c("OTU", "Sample", "Abundance"))

# create a dataframe to match unique row/gene names to Broadclass since that is desired name
#DNAK - Metals
#SME - Multi-Compound
#TBTB - Biocides
#TETQ - Drugs
ids <- data.frame(OTU = c("DNAK", "SME", "TBTB", "TETQ"),
                  Broadclass = c("Metals", "Multicompound", "Biocides", "Drugs"))

# merge tables by gene name to get Broadclass information matched with comparative category and abundance
abund <- merge(phy, ids, by = "OTU")

# get abundance as a percent and round to whole numbers
abund$percent <- abund$Abundance * 100
abund$percent <- round(abund$percent, digits = 0)


