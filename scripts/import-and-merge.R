## Importing & Merging Files
## SAB last updated 3/1/2023

#load packages----
library(BiocManager)
library(dplyr)
library(stringr)
library(plyr)
library(tibble)

#read in and clean data----
b1 <- read.csv("amr++outputs/ransom/batch1_AMR_analytic_matrix.csv")
b2 <- read.csv("amr++outputs/ransom/batch2_AMR_analytic_matrix.csv")
b3 <- read.csv("amr++outputs/ransom/batch3_AMR_analytic_matrix.csv")
b4 <- read.csv("amr++outputs/ransom/batch4_AMR_analytic_matrix.csv")
b5 <- read.csv("amr++outputs/ransom_shotgun_723_AMR_analytic_matrix.csv")

b1 <- column_to_rownames(b1, var = "X")
b2 <- column_to_rownames(b2, var = "X")
b3 <- column_to_rownames(b3, var = "X")
b4 <- column_to_rownames(b4, var = "X")
b5 <- column_to_rownames(b5, var = "X")

#do this for all batches, n being the batch n 
df1 <- cbind(b1, as.data.frame(str_split_fixed(rownames(b1), "\\|", 5))[5]) %>%
  unique() %>%
  aggregate(.~V5, FUN = sum)
df2 <- cbind(b2, as.data.frame(str_split_fixed(rownames(b2), "\\|", 5))[5]) %>%
  unique() %>%
  aggregate(.~V5, FUN= sum)
df3 <- cbind(b3, as.data.frame(str_split_fixed(rownames(b3), "\\|", 5))[5]) %>%
  unique() %>%
  aggregate(.~V5, FUN= sum)
df4 <- cbind(b4, as.data.frame(str_split_fixed(rownames(b4), "\\|", 5))[5]) %>%
  unique() %>%
  aggregate(.~V5, FUN= sum)
df5 <- cbind(b5, as.data.frame(str_split_fixed(rownames(b5), "\\|", 5))[5]) %>%
  unique() %>%
  aggregate(.~V5, FUN= sum)


install.packages("plyr")
library(plyr)
merged <- join_all(list(df1, df2, df3, df4, df5), by = "V5", type = "full")
merged[is.na(merged)] <- 0 #replace NAs 
rownames(merged) <- merged$V5 #fix rownames
colnames(merged) <- trimws(colnames(merged), whitespace = "_.*") #fix colnames
merged_clean <- merged[2:ncol(merged)]

# write.table(merged_clean, file = "data/ransom/ransom-countmatrix-cleanedall.txt", sep = "\t") #save this matrix


#break up gene information for alllll batches - this will be useful for later 
#install.packages("rlist")
library("rlist")
genes <- as.data.frame(list.append(rownames(b1),
                                   rownames(b2),
                                   rownames(b3),
                                   rownames(b4),
                                   rownames(b5))) %>%
  unique()

colnames(genes)[1] <- "allmeg"
genes <- as.data.frame(str_split_fixed(genes$allmeg, "\\|", 5))
colnames(genes) <- c("MEG_ID", "Broadclass", "Class", "Mechanism", "Gene")

# write.table(genes, file = "data/ransom/ransom-geneinfo-all.txt", sep = "\t") 

