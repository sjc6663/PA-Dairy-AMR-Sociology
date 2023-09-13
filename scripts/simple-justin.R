
set.seed(81299)

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

ps2 <- rarefy_even_depth(ps)

# alpha diversity ----
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps2, measures = "Shannon"),
  "MF" = phyloseq::sample_data(ps2)$Male.Female,
  "batch" = phyloseq::sample_data(ps2)$Batch,
  "run" = phyloseq::sample_data(ps2)$Run,
  "aff" = phyloseq::sample_data(ps2)$affiliation,
  "type" = phyloseq::sample_data(ps2)$Conventional.Organic
)

# test variance
varMF <- var.test(Shannon ~ MF, data = adiv, 
                  alternative = "two.sided")
varMF # p = 0.5, ns

# statistical test - WHERE WE ARE HAVING THE ISSUE/QUESTION

# what I was running (t-test with assumed unequal variances): 
wtestMF <- t.test(Shannon ~ MF, data = adiv, var.equal = FALSE)
wtestMF # p = 0.0007**, SIGNIFICANT

# what I think we need to change it to to account for run, batch, relgious affiliation, and farm production type:
model <- lm(Shannon ~ MF + run + batch + aff + type + MF*run*batch*aff*type, data = adiv)
out <- summary(model) 
# significance at MFMale:affEnglish (P = 0.02*)
# MFMale, p = 0.70, ns
