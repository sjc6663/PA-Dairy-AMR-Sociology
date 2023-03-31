## ---- create data frames for each breed and type ----
adiv <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "Conventional.Organic" = phyloseq::sample_data(ps)$Conventional.Organic,
)
head(adiv)

## ---- Shannon ANOVA ----
test <- aov(Shannon ~ Conventional.Organic, data = adiv)
summary(test) # p = 0.0247*

plot_richness(ps, x="Male.Female", measures=c("Shannon"), title = "A", color = "Towel.Type") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1", "#999999")) +
  theme_classic()
  # theme(legend.position = "none")
A
ggsave(plot = A, filename = "plots/alpha-diversity-holstein-fecal.pdf", dpi = 600)


# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")
otu_table(psrel)

# clr transform phyloseq objects at Genus level
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()

dist_mat <- phyloseq::distance(transps, method = "euclidean")

vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female) # p = 0.072
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Conventional.Organic) # p = 0.15
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Group) # p = 0.001***
##  PCA plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Group", shape = "Group") +
  #scale_colour_brewer(palette = "BrBG") +
  stat_ellipse(aes(group = Group, color = Group)) +
  theme_classic() +
  ggtitle("A") +
  #theme(legend.position = "none") +
  labs(caption = "")
