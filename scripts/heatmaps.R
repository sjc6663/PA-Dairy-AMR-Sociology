## TRYING TO DO HEATMAPS
## 4/4/23 SAB


plot_heatmap(ps, sample.label = "Group", taxa.label = "Broadclass")

heatmap(ps, scale = "column")


df.matrix <- data.matrix(meltdf)

heatmap(df.matrix, Rowv = NA,
        Colv = NA, col = color_palette, scale = "column")


ggplot(data = meltdf, aes(x = meltdf$Male.Female, y = meltdf$Conventional.Organic)) +
  geom_tile(aes(fill = meltdf$Number.Milking.Cows)) +
  coord_flip() +
  scale_fill_gradient(low = "#BCE1FF", high = "#56B4E9") +
  theme_minimal()
