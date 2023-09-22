# calculate the lm for x + y interaction and use the residuals from that to do a comparison of x + z interaction

# a residual is the difference of the actual minus the estimated for a y value at a given x
library(tibble)

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

# alpha diversity ----
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "MF" = phyloseq::sample_data(ps)$Male.Female,
  "batch" = phyloseq::sample_data(ps)$Batch,
  "run" = phyloseq::sample_data(ps)$Run,
  "aff" = phyloseq::sample_data(ps)$affiliation,
  "type" = phyloseq::sample_data(ps)$Conventional.Organic,
  "age" = phyloseq::sample_data(ps)$Group,
  "language" = phyloseq::sample_data(ps)$Cultural.Language.Barriers,
  "employees" = phyloseq::sample_data(ps)$employees
)

model <- lm(Shannon ~ aff + run + batch + type, data = adiv)
residuals <- resid(model)

res <- as.data.frame(residuals)
res <- rownames_to_column(res, "samp")

adiv <- rownames_to_column(adiv, "samp")
dat1 <- merge(adiv, res, by = "samp")

dat1 <- column_to_rownames(dat1, "samp")

model2 <- lm(residuals ~ MF, data = dat1)
summary(model2)

model3 <- lm(residuals ~ age, data = dat1)
summary(model3)

model4 <- lm(residuals ~ language, data = dat1)
summary(model4)

model5 <- lm(residuals ~ employees, data = dat1)
summary(model5)
