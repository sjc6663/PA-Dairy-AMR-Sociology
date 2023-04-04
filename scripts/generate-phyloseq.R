#Integration of AMR++ output into phyloseq object for downstream analysis 
#Stephanie Bierly - 3/15/23 
​
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
​
​
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
countsDF <- read.delim("data/ransom/ransom-countmatrix-cleanedall.txt", sep = "\t") 
met <- read_excel("/Users/stephanieclouser/OneDrive - The Pennsylvania State University/Shared-Projects/Ransom-AMR/3.Sample-Collection/Ransom-AMR-Metadata.xlsx")
genes <- read.delim("data/ransom/ransom-geneinfo-all.txt", sep = "\t") %>%
  select(-c(MEG_ID)) %>%
  unique()
​
counts <- as.matrix(countsDF)
otutab <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
​
tax <- as.matrix(genes)
rownames(tax) <- genes$Gene
taxtab <- phyloseq::tax_table(tax)
taxa_names(taxtab)
​
samp <- phyloseq::sample_data(met)
rownames(samp) <- met$`Sample-ID`
​
ps <- phyloseq::phyloseq(otutab, taxtab, samp)

ps <- phyloseq::phyloseq(samp)

#save the ps object 
saveRDS(ps, "data/ransom/rawps.rds")
​
​
​
​
#save work ----
save.image("data/making-psobj.RData")
​
​
​
​
#probably next week or wed----
df <- as.data.frame(sample_data(ps))
df$GeneCounts <- sample_sums(ps)
df <- df[order(df$GeneCounts), ]
df$Index <- seq(nrow(df))
p_genecount <- ggplot(data = df, aes(x = Index, y = GeneCounts, color = Male.Female)) +
  geom_point() +
  ggtitle("Gene Counts")
p_genecount

# test to see what samples came through and if we need to do decontam----

plot_bar(ps, x = "Male.Female")
sample_data(ps)$Sample.Type # Sample, Extra-Sample, Positive Control, and Negative Control 

pc <- subset_samples(
  ps,
  Sample.Type == "Positive-Control"
)
otu_table(pc) # there are a bunch of genes that have counts 
plot_bar(pc, x = "Sample.ID")

nc <- subset_samples(
  ps,
  Sample.Type == "Negative-Control"
) 
otu_table(nc) # there is one gene that I can see that has counts A16S|RequiresSNPConfirmation R80 - 31, R79 - 2

MINREADS = 1

ncdf <- as.data.frame(sample_data(nc))
plot_bar(nc, x = "Sample.ID")
