#Integration of AMR++ output into phyloseq object for downstream analysis 
#Stephanie Bierly - 11/16/23 

#laod packges
library(BiocManager)
library(stringr)
library(tidyverse)
# install.packages("tidylog")
library(tidylog)
# BiocManager::install("phyloseq")
library(phyloseq)
# BiocManager::install("vegan", force = TRUE)
library(vegan) # version 2.6-4
library(ggplot2)
library(readxl)

#read in necessary files: count matrix, gene info, metadata---- 
countsDF <- read.delim("bovine-host-resistome/ransom-countmatrix-cleanedall.txt", sep = "\t") 
met <- read_excel("/Users/stephanieclouser/OneDrive - The Pennsylvania State University/Shared-Projects/Ransom-AMR/3.Sample-Collection/Ransom-AMR-Metadata.xlsx")
genes <- read.delim("bovine-host-resistome/ransom-geneinfo-all.txt", sep = "\t") %>%
  select(-c(MEG_ID)) %>%
  unique()

counts <- as.matrix(countsDF)
otutab <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)

tax <- as.matrix(genes)
rownames(tax) <- genes$Gene
taxtab <- phyloseq::tax_table(tax)
taxa_names(taxtab)

samp <- phyloseq::sample_data(met)
rownames(samp) <- met$`Sample-ID`

ps <- phyloseq::phyloseq(otutab, taxtab, samp)


#save the ps object 
saveRDS(ps, "bovine-host-resistome/rawps.rds")

#save work ----
save.image("data/making-psobj.RData")
