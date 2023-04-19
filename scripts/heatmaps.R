## TRYING TO DO HEATMAPS
## 4/4/23 SAB

library(phyloseq)
library(tidyr)
library(dplyr)
library(data.table)
library(microViz)
library(pheatmap)
library(AggregateR)
library(stringr)
library(readxl)
library(tibble)
library(viridis)


### GOOD ##############################################################################################################################
# https://towardsdatascience.com/pheatmap-draws-pretty-heatmaps-483dab9a3cc

ps <- readRDS("data/full-run/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

#save the ps object 
saveRDS(psrel, "data/full-run/relabund-ps.rds")

## WRANGLING GENE METADATA FOR HEATMAP ----
# create data frame from phyloseq object
meltdf2 <- psmelt(psrel)

# remove any genes that have 0 relative abundance for samples
meltdf3 <- meltdf2[meltdf2$Abundance != 0.000000000000, ]

# from the data table, select only the columns we want for gene info and abundance
genemeta <- dplyr::select(meltdf3, Sample.ID, Broadclass, Class, Abundance) 

# convert to data frame
genemeta <- as.data.frame(genemeta)

# combine all columns that are the same resistance type for each sample
genemeta <- setDT(genemeta)[, .(relabund = sum(Abundance)), by = list(Sample.ID, Broadclass, Class)]

# select only the gene metadata tables we want (Broadclass and Class)
genemeta <- dplyr::select(genemeta, Broadclass, Class)

# combine all groups that are the same
test <- Aggregate(genemeta, c("Class", "Broadclass"))
genemeta <- dplyr::select(test, Broadclass, Class)

genemeta$Broadclass <- as.factor(genemeta$Broadclass) # change broadclass to a factor
genemeta <- as.data.frame(genemeta) # make into a data frame so it can be called for annotation_row =
rownames(genemeta) <- genemeta$Class
genemeta2 <- genemeta[order(genemeta$Broadclass, genemeta$Class),] # sort this so rows will cluster by broadclass and not alphabetical

#fix row names for plot
gm2names <- rownames(genemeta2)
gm2names <- str_replace_all(gm2names, "_", " ") #sub underscores with spaces
gm2names <- str_replace(gm2names, "resistance", "")#drop resistance because that's a given
gm2names <- str_replace_all(gm2names, "and", "+") #simplfy this to make shorter

rownames(genemeta2) <- gm2names

colnames(genemeta2)[1] <- "Resistance Type"
genemeta2

# write.table(genemeta2, file = "tables/genemetadata_cleaned.txt", sep = "\t")

## WRANGLING SAMPLE METADATA FOR HEATMAP ----

met <- read_excel("/Users/stephanieclouser/OneDrive - The Pennsylvania State University/Shared-Projects/Ransom-AMR/3.Sample-Collection/Ransom-AMR-Metadata.xlsx")

metfilt <- met %>% 
  filter(`Sample-Type` == "Sample")

# edit metadata so that it is numeric so that it can be used in the heatmap
is.character(metfilt$`Male/Female`)
metfilt$`Male/Female`[metfilt$`Male/Female` == "Male"] <- 1
metfilt$`Male/Female`[metfilt$`Male/Female` == "Female"] <- 2
is.numeric(metfilt$`Male/Female`)
metfilt$`Male/Female` <- as.numeric(metfilt$`Male/Female`)
is.numeric(metfilt$`Male/Female`)


is.character(metfilt$`Herd-Size`)
metfilt$`Herd-Size`[metfilt$`Herd-Size` == "0-50"] <- 0
metfilt$`Herd-Size`[metfilt$`Herd-Size` == "50-100"] <- 1
metfilt$`Herd-Size`[metfilt$`Herd-Size` == "100-150"] <- 2
metfilt$`Herd-Size`[metfilt$`Herd-Size` == "150-200"] <- 3
metfilt$`Herd-Size`[metfilt$`Herd-Size` == "200+"] <- 4
is.numeric(metfilt$`Herd-Size`)
metfilt$`Herd-Size` <- as.numeric(metfilt$`Herd-Size`)
is.numeric(metfilt$`Herd-Size`)


# column names for aesthetics
colnames(metfilt)[6] <- "Male v. Female Operator"
colnames(metfilt)[54] <- "Milking Herd Size"

rownames(metfilt) <- metfilt$`Sample-ID`

# make into dataframe so it can be called in annotation_col = 
dfmetfilt <- as.data.frame(metfilt)

metfilt$`Milking Herd Size` <- factor(metfilt$`Milking Herd Size`,
                         levels = c("0-50", "50-100", "100-150", "150-200", "200+"))


# add Herd Size and Male Female annotation
annoHS <-data.frame(row.names=metfilt$`Sample-ID`, HerdSize=metfilt$`Milking Herd Size`, MaleFemale=metfilt$`Male v. Female Operator`)

# add colors for each group
newCols <- colorRampPalette(grDevices::rainbow(length(unique(annoHS$HerdSize))))
newCols2 <- colorRampPalette(grDevices::rainbow(length(unique(annoHS$MaleFemale))))
annoHSC <- newCols(length(unique(annoHS$HerdSize)))
annoHSC2 <- newCols2(length(unique(annoHS$MaleFemale)))
names(annoHSC) <- unique(annoHS$HerdSize)
names(annoHSC2) <- unique(annoHS$MaleFemale)
annoHSC <- list(category = annoHSC, mf = annoHSC2)


## COUNTS FOR HEATMAP FILL ----
# get counts and genes
counts <- as.data.frame(psrel@otu_table)
gene <- as.data.frame(psrel@tax_table)

gene2 <- dplyr::select(gene, Broadclass, Class)

counts2 <- rownames_to_column(counts, var = "Gene")

# combine counts and gene info
merge <- merge(counts2, gene)
# write.table(merge, file = "tables/merged-counts-and-gene-clean.txt", sep = "\t")

# select only the columns we want (samples and class)
merge2 <- merge %>% dplyr::select(R23:Class)
merge3 <- merge2 %>% dplyr::select(!Broadclass)

# collapse the classes so that they are only one row
out3 <- aggregate(.~Class, data = merge3, FUN = sum)
# write.table(out3, file = "tables/clr-relabund_byclass.txt", sep = "\t")

counts2 <- column_to_rownames(out3, "Class")
# write.table(counts2, file = "tables/class-count-data-only.txt", sep = "\t")

# fix rowname issue with the data
out3 <- rownames_to_column(out3, "XXX")
out3 <- column_to_rownames(out3, "Class")
out3 <- dplyr::select(out3, !XXX)
out3

# Aesthetics
#fix rownames like with genemeta2
rownames(out3) <- str_replace_all(rownames(out3), "_", " ")#sub underscores with spaces
rownames(out3) <- str_replace(rownames(out3), "resistance", "")#drop resistance because that's a given
rownames(out3) <- str_replace_all(rownames(out3), "and", "+") #simplfy this to make shorter
out3

#now that our rownames match lets sort them to match order of gm2names
mat <- out3[match(gm2names, rownames(out3)),]
mat
mat <- as.matrix(mat)

#set color vars
colors <- list(
  "Male v. Female Operator" = c(Male = "#088F8F", Female = "#5F9EA0"),
  #"Group" = c(Calves = "#1F51FF", Cows = "#1434A4"),
 "Milking Herd Size" = c("0-50" = "#BCE1FF", "50-100" = "#7DCCFF", "100-150" = "#56B4E9", "150-200" = "#098BD9", "200+" = "#4292C6"),
  "Resistance Type" = c("Biocides" = "#552F7A",
                   "Drugs" = "#7C5F98",
                   "Metals" = "#B09FC1",
                   "Multi-compound" = "#CABED6"))

# Plot ----

test <- pheatmap(mat, 
         annotation_col = annoHS,
         annotation_row = subset(x=genemeta2, select = "Resistance Type"),
         cluster_cols = hclust(dist(t(out3), method = "euclidean")),
         cluster_rows = F,
         annotate_colors = colors,
         color = mako(15),
         gaps_row = c(10,27,43),
         fontsize = 15,
         annotation_legend = T,
         legend = T,
         show_rownames = TRUE,
         border_color = NA)

pdf("test2.pdf", width = 18, height = 24)
test
dev.off()

