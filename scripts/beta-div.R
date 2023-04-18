## BETA DIVERSITY

ps <- readRDS("data/ransom/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")
otu_table(psrel)

# clr transform phyloseq objects
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()

OTU <- as(sample_data(transps), "matrix")
OTU2 <- as.data.frame(OTU)
head(OTU2)
OTU3 <- rownames_to_column(OTU2, var = "sample.id.2")


dist_mat <- phyloseq::distance(transps, method = "euclidean")

vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female) # p = 0.065, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Group) # p = 0.001***
vegan::adonis2(dist_mat ~ Group*Male.Female, data = OTU3) # p (Group:Male.Female) = 0.489, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Herd.Size) # p = 0.632, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Conventional.Organic) # p = 0.18, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Formal.Team.Meetings.Frequency) # p = 0.96, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Cultural.Language.Barriers) # p = 0.967, not sig




##  PCA plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Conventional.Organic", shape = "Conventional.Organic") +
  scale_color_manual(values = color_palette) +
  stat_ellipse(aes(group = Conventional.Organic, color = Conventional.Organic)) + 
  theme_classic() +
  ggtitle("B") + 
  labs(caption = "")

ggsave(filename = "plots/PCA-microViz-angus-fecal.pdf", dpi = 600)


