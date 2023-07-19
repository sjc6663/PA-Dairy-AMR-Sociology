# Beta diversity of clinically significant genes
# SAB 07/17/2023

color_palette <- c("#367aa1", "#348fa7", "#40b7ad", "#8ad9b1", "#def4e5")
small_color <- c("#367aa1", "#def4e5")

set.seed(81299)

# read in phyloseq object
ps <- readRDS("data/full-run/sig-decontam-ps.RDS")

ps2 <- rarefy_even_depth(ps)

# filter out taxa that are not clinically significant and remove them
ps <- subset_taxa(ps, 
                  sig == "TRUE")

ps

ps <- ps %>% 
  ps_mutate(
    employees = if_else(str_detect(Non.Family.Milkers, "0"), true = "No", false = "Yes")
  ) 

ps <- ps %>% 
  ps_mutate(
    meetings = if_else(str_detect(Formal.Team.Meetings.Frequency, "Never"), true = "No", false = "Yes")
  ) 
# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")


# clr transform phyloseq objects
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()

dist_mat <- phyloseq::distance(transps, method = "euclidean")

vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female) # p = 0.225, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Group) # p = 0.001 ***
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female*phyloseq::sample_data(transps)$Group) # p = 0.62, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Formal.Team.Meetings.Frequency) # p = 0.31, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Cultural.Language.Barriers) # p = 0.2, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Herd.Size) # p = 0.117, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$employees) # p = 0.176, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$meetings) # p = 0.694, ns


# PCA Plots ---------------------------------------------------------------------------------------------------------

psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Conventional.Organic") +
  scale_color_viridis(option = "viridis", discrete = TRUE) +
 # stat_ellipse(aes(group = meetings, color = meetings)) + 
  theme_classic() +
  ggtitle("A") + 
  labs(caption = "R2 = 0.016, F(1,70) = 1.07, P = 0.30")

