ps <- readRDS("data/full-run/decontam-ps.rds")

sample_data(ps)

english <- subset_samples(ps,
                          affiliation == "English")

sample_data(english)

# create data frame with relevant metadata for comparison
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(english, measures = "Shannon"),
  "MF" = phyloseq::sample_data(english)$Male.Female, 
  "run" = phyloseq::sample_data(english)$Run, 
  "batch" = phyloseq::sample_data(english)$Batch
)

# # male female
wtestMF <- lm(Shannon ~ MF + run + batch + MF*run*batch, data = adiv)
summary(wtestMF) # p = 0.91, ns

# transform to relative abundance
psrel <- microbiome::transform(english, "compositional")

# clr transform
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()

# create distance matrix 
dist_mat <- phyloseq::distance(transps, method = "euclidean")

# test
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female*phyloseq::sample_data(transps)$Batch*phyloseq::sample_data(transps)$Run) # R2 = 0.021, F(1,70) = 1.49, P = 0.006**

psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Male.Female", size = 6, axes = c(1,3)) +
  scale_color_manual(values = c("#367aa1", "#def4e5")) +
  stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  labs(color = "Gender") +
  ggtitle("C") + 
  labs(caption = "R2 = 0.021, F(1,70) = 1.49, P = 0.006**") +
  theme(text = element_text(size = 15)) 

