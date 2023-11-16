## Importing & Merging Files
## SAB last updated 11/15/2023

#load packages----
library(BiocManager)
library(dplyr)
library(stringr)
library(plyr)
library(tibble)

#read in and clean data----
b1 <- read.csv("amr++outputs/bovine-host/AMR_analytic_matrix_b1.csv")
b2 <- read.csv("amr++outputs/bovine-host/AMR_analytic_matrix_b2.csv")
b3 <- read.csv("amr++outputs/bovine-host/AMR_analytic_matrix_b3.csv")
b4 <- read.csv("amr++outputs/bovine-host/AMR_analytic_matrix_b4.csv")

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

# install.packages("plyr")
library(plyr)
merged <- join_all(list(df1, df2, df3, df4), by = "V5", type = "full")
merged[is.na(merged)] <- 0 #replace NAs 
rownames(merged) <- merged$V5 #fix rownames
colnames(merged) <- trimws(colnames(merged), whitespace = "_.*") #fix colnames
merged_clean <- merged[2:ncol(merged)]

# write.table(merged_clean, file = "bovine-host-resistome/ransom-countmatrix-cleanedall.txt", sep = "\t") #save this matrix


#break up gene information for alllll batches - this will be useful for later 
#install.packages("rlist")
library("rlist")
genes <- as.data.frame(list.append(rownames(b1),
                                   rownames(b2),
                                   rownames(b3),
                                   rownames(b4))) %>%
  unique()

colnames(genes)[1] <- "allmeg"
genes <- as.data.frame(str_split_fixed(genes$allmeg, "\\|", 5))
colnames(genes) <- c("MEG_ID", "Broadclass", "Class", "Mechanism", "Gene")

# write.table(genes, file = "bovine-host-resistome/ransom-geneinfo-all.txt", sep = "\t") 

