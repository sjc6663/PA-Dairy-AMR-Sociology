## Wrangling AMR++ Output
## SAB 02/13/2023

# load packages
library(BiocManager)
library(dplyr)
library(stringr)

# read in and clean data
out <- read.csv("AMR_analytic_matrix.csv")

# drop extension from file name
colnames(out) <- str_split_i(colnames(out), "_", 1)

# save work
save.image("data/wrangling.RData")

# break up gene information
genes <- as.data.frame(rownames(out))
genes <- as.data.frame(str_split_fixed(genes$`rownames(out)`, "\\|", 5))
# give the column names
colnames(genes) <- c("MEG_ID", "Broadclass", "Class", "Mechanism", "Gene")

out2 <- cbind(out, genes[5])

# collapse the genes so that they are only one row
out3 <- aggregate(.~Gene, data = out2, FUN = sum)
# make the table super pretty
out4 <- out3[2:ncol(out3)]
# change the rownames to be the genes
rownames(out4) <- out3$Gene
# save the file 
write.csv(out4, file = "data/countmatrix-cleaned.csv")

## ---- Prelim Plots ----
# load packages
library(ggplot2)
library(data.table)
library(microViz)

# put into a data table to make it easier to plot
# you can add the other information to your out3 table (meaning drug class, class, sample metadata, etc.) and then melt and it will include all of the other info
meltdf <- melt(as.data.table(out3))

# plot
ggplot(meltdf %>% 
        filter(value != "0"), aes(x = variable, fill = variable, y = value)) +
  geom_bar(stat = "identity") + 
  scale_y_log10() +
  facet_wrap(~Gene) +
  theme(legend.position = "none")

# make a heatmap
ggplot(meltdf %>% 
         filter(value != "0"), aes(x = Gene, fill = log(value), y = variable)) +
  geom_tile() +
  theme_bw() +
  coord_fixed() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  




# random line of code to keep the space at the bottom so I'm not typing tight against the console.