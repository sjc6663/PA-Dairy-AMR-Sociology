#Integration of AMR++ output into phyloseq object for downstream analysis 
#Stephanie Bierly - 3/15/23 

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

#my filtering mess that yall dont need to see---- 
countsDF <- read.delim("tables/countmatrix-cleanedall.txt", sep = "\t") 
met <- read.delim("batchinfo/meta/meta-both.txt")
genes <- read.delim("tables/geneinfo-all.txt")
​
​
samples <- colnames(countsDF)
metfilt <- met %>%
  filter(srr %in% samples) #there are 8 extras 
​
dropped <- met %>%
  filter(!srr %in% samples) #these are expected 
​
dupecheck <- data.frame(table(met$srr)) %>%
  filter(Freq > 1)  #find duplicates
​
dupelist <- dupes$srr %>% unique()
​
dupes <- met %>%
  filter(srr %in% dupecheck$Var1) #find batch no associated and save this
​
#write.table(dupes, file = "tables/srr-dupes-batchno.txt", sep = "\t") #save this matrix
#for sake of tutorial fix metfilt to remove the dupes and go back later to sort which is appropriate to drop
​
duperm <- metfilt %>%
  filter(!srr %in% dupecheck$Var1)
​
test <- countsDF[, -grep(paste(dupelist[1]), colnames(countsDF))]
test <- test[, -grep(paste(dupelist[2]), colnames(test))]
test <- test[, -grep(paste(dupelist[3]), colnames(test))]
test <- test[, -grep(paste(dupelist[4]), colnames(test))]
test <- test[, -grep(paste(dupelist[5]), colnames(test))]
test <- test[, -grep(paste(dupelist[6]), colnames(test))]
test <- test[, -grep(paste(dupelist[7]), colnames(test))]
test <- test[, -grep(paste(dupelist[8]), colnames(test))]
​
#save these for working through amr++ tutorial 
#write.table(duperm, "batchinfo/meta/meta-both-duperm.txt", sep = "\t")
#write.table(test, "tables/countmatrix-duperm.txt", sep = "\t")
​
#read in necessary files: count matrix, gene info, metadata---- 
countsDF <- read.delim("/Users/stephanieclouser/OneDrive - The Pennsylvania State University/Shared-Projects/Ransom-AMR/7.Data-Analysis/SK_Prelim/countmatrix-cleanedall.txt", sep = "\t") 
met <- read_excel("/Users/stephanieclouser/OneDrive - The Pennsylvania State University/Shared-Projects/Ransom-AMR/3.Sample-Collection/Ransom-AMR-Metadata.xlsx")
genes <- read.delim("/Users/stephanieclouser/OneDrive - The Pennsylvania State University/Shared-Projects/Ransom-AMR/7.Data-Analysis/SK_Prelim/geneinfo-all.txt", sep = "\t") %>%
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
saveRDS(ps, "data/full-run/rawps.rds")

#save work ----
save.image("data/making-psobj.RData")