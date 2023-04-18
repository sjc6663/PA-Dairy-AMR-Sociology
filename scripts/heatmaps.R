## TRYING TO DO HEATMAPS
## 4/4/23 SAB


plot_heatmap(ps, sample.label = "Group", taxa.label = "Broadclass")

heatmap(ps, scale = "column")

sample_data(psrel)
df.matrix <- data.matrix(meltdf)

heatmap(df.matrix, Rowv = NA,
        Colv = NA, col = color_palette, scale = "column")


meltdf2 <- psmelt(psrel)

new <- meltdf %>% select(c("OTU", "Sample", "Abundance"))
new.mat <- data.matrix(meltdf3)


ggplot(meltdf, aes(x = Sample, y = OTU, fill = Abundance)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# BASE R PACKAGE

heatmap(out2,
        xlab = "Sample",
        ylab = "Abundance")
otu_table(ps)

out <- as.data.frame(otu_table(ps))
out2 <- data.matrix(out)




all_site_v4

#filter out rows with 0 relabun 
meltdf3 <- meltdf2[meltdf2$Abundance != 0.0000000000, ]

#set color vars
colors <- list(
  "Male.Female" = c(Male = "#000004FF", Female = "#24868EFF"),
  "Group" = c(Calves = "#7DCCFF", Cows = "#098BD9"),
  "Broadclass" = c("Biocides" = "#552F7A",
                        "Drugs" = "#7C5F98",
                        "Metals" = "#B09FC1",
                        "Multi-compound" = "#CABED6"))


pheatmap(new.mat, 
           annotation_col = subset(x=meltdf3, select = c("Male.Female", "Group")),
           annotation_row = subset(x=meltdf3, select = "Broadclass"),
           #cluster_cols = hclust(dist(t(out), method = "euclidean")), 
           cluster_rows = F,
           annotation_colors = colors,
           color = inferno(15),
           gaps_row = c(8,31,43),
           fontsize = 15,
           annotation_legend = F,
           legend = T,
           show_rownames = TRUE,
           border_color = "black",
          show_colnames = FALSE)










##############################################################################################################################

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

meltdf <- psmelt(ps)

#filter out rows with 0 relabun 
out2 <- meltdf[meltdf$Abundance != 0.0000000000, ]


#first row annotations
genemeta <- out2 %>% dplyr::select(Sample, Broadclass, Class, Abundance) %>% 
  group_by(Sample, Class)%>%
  pivot_wider(values_from = Abundance) %>% 
  dplyr::select(Broadclass, Class)

genemeta %>% pivot_wider(values_from = Abundance)


  genemeta$Broadclass <- as.factor(genemeta$Broadclass) #change broadclass to factor
  genemeta <- as.data.frame(genemeta) #make into df so it can be called for annotation_row =
  rownames(genemeta) <- genemeta$Class #convert rownames to match o6
  genemeta2 <- genemeta[order(genemeta$Broadclass, genemeta$Class),] #sort this so rows will cluster by broadclass and not just be plotted alphabetically

  
#fix row names for plot
gm2names <- rownames(genemeta2)
gm2names <- str_replace_all(gm2names, "_", " ") #sub underscores with spaces
gm2names <- str_replace(gm2names, "resistance", "")#drop resistance because that's a given
gm2names <- str_replace_all(gm2names, "and", "+") #simplfy this to make shorter

rownames(genemeta2) <- gm2names #assign back
#fix colnames for aesthetic purposes
colnames(genemeta2)[1] <- "Sample"
colnames(genemeta2)[3] <- "Resistance Type"    
genemeta2

metfilt <- met %>% 
  filter(`Sample-Type` == "Sample")

metfilt

#fix rownames like with genemeta2
rownames(out2) <- str_replace_all(rownames(out2), "_", " ")#sub underscores with spaces
rownames(o6) <- str_replace(rownames(o6), "resistance", "")#drop resistance because that's a given
rownames(o6) <- str_replace_all(rownames(o6), "and", "+") #simplfy this to make shorter


#set color vars
colors <- list(
  "Male.Female" = c(Male = "#000004FF", Female = "#24868EFF"),
  "Group" = c(Calves = "#7DCCFF", Cows = "#098BD9"),
  "Broadclass" = c("Biocides" = "#552F7A",
                   "Drugs" = "#7C5F98",
                   "Metals" = "#B09FC1",
                   "Multi-compound" = "#CABED6"))


mat <- out2[match(gm2names, rownames(out2)),]

  pheatmap(meltdf, 
           annotation_col = subset(x=meltdf, select = c("Male.Female", "Herd.Size")),
           annotation_row = subset(x=meltdf, select = "Class"),
           #cluster_cols = hclust(dist(t(out2), method = "euclidean")), 
          # cluster_rows = F,
           annotation_colors = colors,
           color = inferno(15),
           gaps_row = c(8,31,43),
           fontsize = 15,
           annotation_legend = F,
           legend = T,
           show_rownames = TRUE,
           border_color = NA)


df_num <- as.matrix(genemeta2[4])

### GOOD ##############################################################################################################################

# get counts and genes
counts <- as.data.frame(ps@otu_table)
gene <- as.data.frame(ps@tax_table)

counts2 <- rownames_to_column(counts, var = "Gene")

merge <- merge(counts2, gene)
write.table(merge, file = "tables/merged-counts-and-gene-clean.txt", sep = "\t")

merge2 <- merge %>% dplyr::select(R23:Class)
merge3 <- merge2 %>% dplyr::select(!Broadclass)

# collapse the classes so that they are only one row
out3 <- aggregate(.~Class, data = merge3, FUN = sum)

counts2 <- column_to_rownames(out3, "Class")
write.table(counts2, file = "tables/class-count-data-only.txt", sep = "\t")

df_num_scale <- scale(counts2)
pheatmap(df_num_scale, 
         #annotation_col = subset(x = mat, select = c("Male.Female", "Herd.Size")),
         annotation_row = subset(x = merge, select = "Broadclass"),
         cluster_rows = F,
         color = mako(15),
         annotation_legend = F,
        # gaps_row = c(8, 31, 43),
         main = "pheatmap test")




