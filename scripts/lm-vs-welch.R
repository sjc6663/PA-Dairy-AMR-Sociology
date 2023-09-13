ps <- readRDS("data/full-run/decontam-ps.rds")

# create data frame with relevant metadata for comparison
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "MF" = phyloseq::sample_data(ps)$Male.Female, 
  "batch" = phyloseq::sample_data(ps)$Batch,
  "run" = phyloseq::sample_data(ps)$Run,
  "aff" = phyloseq::sample_data(ps)$affiliation
)

# test variance
# Male Female 
varMF <- var.test(Shannon ~ MF, data = adiv, 
                  alternative = "two.sided")
varMF # p = 0.28, not sig

# # male female
wtestMF <- lm(Shannon ~ MF + run + batch + MF*run*batch, data = adiv)
out1 <- summary(wtestMF) # p = 0.91, ns

wt <- t.test(Shannon ~ MF, data = adiv, var.equal = FALSE)
wt

write.csv(out1, file = "data/full-run/lm-alpha.csv")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# clr transform
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()

# create distance matrix 
dist_mat <- phyloseq::distance(transps, method = "euclidean")

# test
out <- vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female*phyloseq::sample_data(transps)$Batch*phyloseq::sample_data(transps)$Run*phyloseq::sample_data(transps)$affiliation) # R2 = 0.021, F(1,70) = 1.49, P = 0.006**

write.csv(out, file = "data/full-run/perm-inter.csv")
