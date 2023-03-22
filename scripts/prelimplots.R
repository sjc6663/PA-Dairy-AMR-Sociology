## Exploratory Plots
## SAB - last updated 2/15/2023

## ---- setup ----
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

## ---- merge data ----

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


## ---- OPTIONS FOR EXPLORATORY PLOTS ----

# make a heatmap
ggplot(meltdf2 %>% 
         filter(counts != "0"), aes(x = Gene, fill = Class, y = rh_ID)) +
  geom_tile() +
  theme_bw() +
  coord_fixed() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))

# make a barplot
ggplot(meltdf2 %>% 
         filter(counts != "0") %>% 
         filter(rh_ID == "rh02"), aes(x = rh_ID, fill = Class, y = counts)) +
  geom_bar(stat = "identity") + 
  scale_y_log10() +
  facet_wrap(~Broadclass) + # can facet wrap by male vs female
  theme(legend.position = "none")



# save work
save.image("data/prelim-plots.RData")
