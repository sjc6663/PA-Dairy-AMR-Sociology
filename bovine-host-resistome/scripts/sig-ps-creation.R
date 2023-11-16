## Edit phyloseq to have clinically significant genes
## SAB 7/5/23

# load packages
library(phyloseq)
library(stringr)
library(microViz)
library(ALDEx2)
library(dplyr)
library(tibble)

siggene <- read.csv("CIA_genes.csv")

siggene

genes <- siggene$SNPConfirmation

gene <- as.data.frame(genes)

ps <- readRDS("bovine-host-resistome/decontam-ps.RDS")

tax <- as.data.frame(ps@tax_table)

tax <- tax %>% 
  dplyr::select(Gene)

colnames(gene) <- c("Gene")

tax$sig <- tax$Gene %in% gene$Gene

tax <- rownames_to_column(tax, "gene2")

tax <- dplyr::select(tax, Gene, sig)

tax <- as.data.frame(tax)

tax3 <- as.data.frame(ps@tax_table)

tax4 <- merge(tax, tax3, by = "Gene")

tax4 <- as.data.frame(tax4)

tax5 <- tax4[, c("sig", "Broadclass", "Class", "Mechanism", "Gene")]

tax5 <- as.matrix(tax5)
rownames(tax5) <- tax$Gene
taxtab <- phyloseq::tax_table(tax5)
taxa_names(taxtab)

tax_table(ps) <- phyloseq(tax_table(taxtab))

tax_table(ps)

saveRDS(ps, "bovine-host-resistome/sig-decontam-ps.rds")


