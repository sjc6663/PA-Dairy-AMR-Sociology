
set.seed(81299)

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

# alpha diversity ----
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "MF" = phyloseq::sample_data(ps)$Male.Female,
  "batch" = phyloseq::sample_data(ps)$Batch,
  "run" = phyloseq::sample_data(ps)$Run,
  "aff" = phyloseq::sample_data(ps)$affiliation, 
  "type" = phyloseq::sample_data(ps)$Conventional.Organic
)


# get residuals
model <- lm(Shannon ~ aff + run + batch + type, data = adiv)
residuals <- resid(model)

res <- as.data.frame(residuals)
res <- rownames_to_column(res, "samp")

adiv <- rownames_to_column(adiv, "samp")
dat1 <- merge(adiv, res, by = "samp")

dat1 <- column_to_rownames(dat1, "samp")

model2 <- lm(residuals ~ MF, data = dat1)
summary(model2) # P = 0.0216*