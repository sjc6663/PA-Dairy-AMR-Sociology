## Exploratory Plots
## SAB - last updated 2/15/2023

# load packages
library(BiocManager)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)


# load data from last time
load("data/wrangling.RData")

# load metadata
meta <- read.delim("/Users/stephanieclouser/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Shared-Projects/rhAMR_preliminary/rhAMR_miseq-analysis/sampleinfo/rhamr-samples.txt", sep = "\t")


# merge metadata with melted df
# first we need to correct the column names in meltdf
colnames(meltdf) <- c("Gene", "rh_ID", "counts")

# then merge them together  
meltdf2 <- merge(meltdf, meta, by = "rh_ID")

# add gene info as well, but we need to remove duplicates of genes otherwise it gets very confused because the rows are the same except for the megID column
genefilt <- genes[2:5] %>% 
  unique() # will drop duplicate rows

# merge genes table with other data
meltdf2 <- merge(meltdf2, genefilt, by = "Gene")




# save work
save.image("data/prelim-plots.RData")