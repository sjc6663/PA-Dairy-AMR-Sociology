broad <- tax_glom(ps, taxrank = "Broadclass")

# alpha diversity ----
adiv2 <- data.frame(
  "Shannon" = phyloseq::estimate_richness(broad, measures = "Shannon"),
  "MF" = phyloseq::sample_data(broad)$Male.Female,
  "batch" = phyloseq::sample_data(broad)$Batch,
  "run" = phyloseq::sample_data(broad)$Run,
  "aff" = phyloseq::sample_data(broad)$affiliation, 
  "type" = phyloseq::sample_data(broad)$Conventional.Organic
)


# get residuals
model3 <- lm(Shannon ~ aff + run + batch + aff + type, data = adiv2)
residuals2 <- resid(model3)

res2 <- as.data.frame(residuals2)
res2 <- rownames_to_column(res2, "samp")

adiv2 <- rownames_to_column(adiv2, "samp")
dat2 <- merge(adiv2, res2, by = "samp")

dat2 <- column_to_rownames(dat2, "samp")

model4 <- lm(residuals2 ~ MF, data = dat2)
summary(model4) # P = 0.0205*



ait2 <- broad %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait2 %>% dist_bdisp(variables = "Male.Female") %>% bdisp_get() # p=0.501

# test with PERMANOVA
mod2 <- ait2 %>%
  dist_permanova(
    seed = 81299,
    variables = "Male.Female + Run + Batch + affiliation + Conventional.Organic",
    n_perms = 9999
  )

mod2 # R2 = 0.02, F(1, 65) = 1.16, P = 0.0287
