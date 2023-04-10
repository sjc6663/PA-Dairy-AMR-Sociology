## BETA DIVERSITY

ps <- readRDS("data/ransom/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")
otu_table(psrel)

# clr transform phyloseq objects at Genus level
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


