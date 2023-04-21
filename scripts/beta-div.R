## BETA DIVERSITY

ps <- readRDS("data/full-run/decontam-ps.RDS")

colors <- c("#BCF5F9", "#89C5FD", "#3A80EC", "#0229BF", "#080B6C")

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

vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female) # p = 0.303, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Group) # p = 0.001***
vegan::adonis2(dist_mat ~ Group*Male.Female, data = OTU3) # p (Group:Male.Female) = 0.559, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Herd.Size) # p = 0.639, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Conventional.Organic) # p = 0.19, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Formal.Team.Meetings.Frequency) # p = 0.96, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Cultural.Language.Barriers) # p = 0.98, not sig


##  PCA plot - Male Female ----
A_soc <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Male.Female") +
  scale_color_manual(values = colors) +
  stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  ggtitle("A") + 
  labs(caption = "R2 = 0.016, F(1,70) = 1.07, P = 0.30")
A_soc

ggsave(filename = "plots/full-run/PCA-male-female.pdf", dpi = 600)

##  PCA plot - Male Female Age Group ----
B_soc <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Group", shape = "Male.Female") +
  scale_color_manual(values = colors) +
  stat_ellipse(aes(group = Group, color = Group)) + 
  theme_classic() +
  ggtitle("B") + 
  labs(caption = "R2 = 0.011, F(1, 68) = 0.85, P = 0.63")

ggsave(filename = "plots/full-run/PCA-male-female-age-group.pdf", dpi = 600)

##  PCA plot - Formal Team Meetings ----

sample_data(psrel)$Formal.Team.Meetings.Frequency <- factor(sample_data(psrel)$Formal.Team.Meetings.Frequency,
                                                            levels = c("Never", "1 or 2 times/year", "Quarterly", "Once a month", "At least twice a month"))


C_soc <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Formal.Team.Meetings.Frequency") +
  scale_color_manual(values = colors) +
  stat_ellipse(aes(group = Formal.Team.Meetings.Frequency, color = Formal.Team.Meetings.Frequency)) + 
  theme_classic() +
  ggtitle("C") + 
  labs(caption = "R2 = 0.057, F(4, 67) = 1.01, P = 0.40")

C_soc

ggsave(filename = "plots/full-run/PCA-team-meetings.pdf", dpi = 600)

##  PCA plot - Cultural Language Barriers ----
D_soc <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Cultural.Language.Barriers") +
  scale_color_manual(values = colors) +
  stat_ellipse(aes(group = Cultural.Language.Barriers, color = Cultural.Language.Barriers)) + 
  theme_classic() +
  ggtitle("D") + 
  labs(caption = "R2 = 0.016, F(1, 70) = 1.14, P = 0.22")

D_soc

ggsave(filename = "plots/full-run/PCA-lang-barriers.pdf", dpi = 600)

## Combine Soc Plots ----
library(patchwork)

(A_soc | B_soc) / (C_soc | D_soc)

ggsave(filename = "plots/full-run/PCA-soc-plots.pdf", dpi = 600, height = 10, width = 12)


##  PCA plot - Farm Type  ----
A_farm <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Conventional.Organic") +
  scale_color_manual(values = colors) +
  stat_ellipse(aes(group = Conventional.Organic, color = Conventional.Organic)) + 
  theme_classic() +
  ggtitle("A") + 
  labs(caption = "R2 = 0.017, F(1, 70) = 1.20, P = 0.17") 

A_farm

ggsave(filename = "plots/full-run/PCA-farm-type.pdf", dpi = 600)

##  PCA plot - Age Group  ----
B_farm <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Group") +
  scale_color_manual(values = colors) +
  stat_ellipse(aes(group = Group, color = Group)) + 
  theme_classic() +
  ggtitle("B") + 
  labs(caption = "R2 = 0.068, F(1, 70) = 5.08, P = 0.001***, ") 

B_farm

ggsave(filename = "plots/full-run/PCA-age-group.pdf", dpi = 600)

##  PCA plot - Herd Size ----

sample_data(psrel)$Herd.Size <- factor(sample_data(psrel)$Herd.Size,
                                       levels = c("0-50", "50-100", "100-150", "150-200", "200+"))

C_farm <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Herd.Size") +
  scale_color_manual(values = colors) +
  stat_ellipse(aes(group = Herd.Size, color = Herd.Size)) + 
  theme_classic() +
  ggtitle("C") + 
  labs(caption = "R2 = 0.078, F(4, 67) = 1.42, P = 0.11") +
  theme(legend.position = "bottom")

C_farm

ggsave(filename = "plots/full-run/PCA-herd-size.pdf", dpi = 600)

## Combine Farm Plots ----

(A_farm | B_farm) / C_farm

ggsave(filename = "plots/full-run/PCA-farm-stats.pdf", dpi = 600, width = 12, height = 10)
