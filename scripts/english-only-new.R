# English only gender

set.seed(81299)
library(phyloseq)

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

eng <- subset_samples(ps,
               affiliation == "English")

# alpha diversity ----
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(eng, measures = "Shannon"),
  "MF" = phyloseq::sample_data(eng)$Male.Female,
  "batch" = phyloseq::sample_data(eng)$Batch,
  "run" = phyloseq::sample_data(eng)$Run,
  "type" = phyloseq::sample_data(eng)$Conventional.Organic
)


# get residuals
model <- lm(Shannon ~ run + batch + type, data = adiv)
residuals <- resid(model)

res <- as.data.frame(residuals)
res <- rownames_to_column(res, "samp")

adiv <- rownames_to_column(adiv, "samp")
dat1 <- merge(adiv, res, by = "samp")

dat1 <- column_to_rownames(dat1, "samp")

model2 <- lm(residuals ~ MF, data = dat1)
summary(model2) # P = 0.0517, Padj = 0.517, ns



ait <- eng %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "Male.Female") %>% bdisp_get() # p=0.0007068

# test with PERMANOVA
mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "Male.Female + Run + Batch + Conventional.Organic",
    n_perms = 9999
  )

mod1 # R2 = 0.04, F(1, 25) = 1.22, P = 0.0145, Padj = 0.145, ns



# get residuals
model2 <- lm(Shannon ~ run + batch + MF, data = adiv)
residuals2 <- resid(model2)

res2 <- as.data.frame(residuals2)
res2 <- rownames_to_column(res2, "samp")

adiv2 <- rownames_to_column(adiv, "samp")
dat12 <- merge(adiv2, res2, by = "samp")

dat12 <- column_to_rownames(dat12, "samp")

model3 <- lm(residuals ~ type, data = dat12)
summary(model3) # P = 0.0517, Padj = 0.517, ns

