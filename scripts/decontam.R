# Decontam
## SAB 03/22/2023 

# load packages ----
library(BiocManager)
library(stringr)
library(tidyverse)
library(tidylog)
library(phyloseq)
library(vegan) #version 2.6-4
library(ggplot2)
library(decontam)
library(readxl)
library(microViz)

# load ps objects and other data ----
countsDF <- read.delim("data/ransom/ransom-countmatrix-cleanedall.txt", sep = "\t") 
met <- read_excel("/Users/stephanieclouser/OneDrive - The Pennsylvania State University/Shared-Projects/Ransom-AMR/3.Sample-Collection/Ransom-AMR-Metadata.xlsx")
genes <- read.delim("data/ransom/ransom-geneinfo-all.txt", sep = "\t") %>%
  select(-c(MEG_ID)) %>%
  unique()

ps <- readRDS("data/ransom/rawps.rds")

# dummy code samples vs controls ----
ps <- ps %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Sample.Type,"Control"), true = "Control", false = "Sample")
  ) 

# inspect # of counts per sample
sam <- as.data.frame(sample_data(ps))
sam$GeneCounts <- sample_sums(ps)
sam <- sam[order(sam$GeneCounts), ]
sam$Index <- seq(nrow(sam))


# plot to look at it 
ggplot(sam, aes(x = Index, y = GeneCounts, color = SampleBinary)) +
  geom_point()


sample_data(ps)$is.neg <- sample_data(ps)$SampleBinary == "Control"
contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg")

# make a phyloseq object of only the contaminants and controls
con <- rownames(contamdf.prev[contamdf.prev$contaminant == "TRUE",])
conps <- prune_taxa(rownames(ps@tax_table) %in% con, ps)
negps <- subset_samples(conps, is.neg == "TRUE")

# remove contaminant sequences
contams <- rownames(contamdf.prev[contamdf.prev$contaminant == "TRUE",])
contams

# remove contaminant sequences
nocontam <- prune_taxa(!rownames(ps@tax_table) %in% contams, ps)
  
# Sophia will figure this out and get it to me
ps_filter(ps, )

# remove all controls
nocontrol <- prune_samples(!str_detect(PUT SOMETHING HERE FROM ABOVE^^, rownames(nocontam@sam_data)), nocontam)


# save work
save.image("data/decontam.RData")