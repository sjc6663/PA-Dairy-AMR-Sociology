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

