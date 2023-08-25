## BETA DIVERSITY

library(phyloseq)
library(vegan)
library(microViz)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)

set.seed(81299)

ps <- readRDS("data/full-run/decontam-ps.RDS")

ps2 <- rarefy_even_depth(ps)

color_palette <- c("#8ad9b1", "#40b7ad", "#348fa7", "#37659e", "#423d7b", "black")
small_color <- c("#348fa7", "#423d7b")

# transform to relative abundance
psrel <- microbiome::transform(ps2, "compositional")

psrel <- psrel %>% 
  ps_mutate(
    orgtime = if_else(str_detect(How.Long.Organic, "Not Organic"), true = "Not Organic", false = "0-5 years")
  ) 

# factor how long farms have been organic into buckets, took a lot of confusing code -----
psrel <- psrel %>% 
  ps_mutate(orgtime = case_when(
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "1 year") ~ "0-5 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "2 years") ~ "0-5 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "3 years") ~ "0-5 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "5 years") ~ "0-5 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "10 years") ~ "6-10 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "12 years") ~ "11-15 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "13 years") ~ "11-15 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "14 years") ~ "11-15 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "16 years") ~ "16-20 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "17 years") ~ "16-20 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "20 years") ~ "16-20 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "20+ years") ~ "> 21 years",
    orgtime == "Not Organic" ~ "Not Organic"
  ))


psrel <- psrel %>% 
  ps_mutate(orgtime = case_when(
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "12 years") ~ "11-15 years",
    orgtime == "0-5 years" & str_detect(How.Long.Organic, "13 years") ~ "11-15 years",
    orgtime == NA ~ "20+ years",
    orgtime == "Not Organic" ~ "Not Organic",
    orgtime == "0-5 years" ~ "0-5 years",
    orgtime == "6-10 years" ~ "6-10 years",
    orgtime == "11-15 years" ~ "11-15 years",
    orgtime == "16-20 years" ~ "16-20 years",
  ))


psrel <- psrel %>% 
  ps_mutate(orgtime = case_when(
   orgtime == NA ~ "20+ years",
    orgtime == "Not Organic" ~ "Not Organic",
    orgtime == "0-5 years" ~ "0-5 years",
    orgtime == "6-10 years" ~ "6-10 years",
    orgtime == "11-15 years" ~ "11-15 years",
    orgtime == "16-20 years" ~ "16-20 years",
  ))

sample_data(psrel)$orgtime[is.na(sample_data(psrel)$orgtime)] <- "20+ years"

sample_data(psrel)$orgtime
sample_data(psrel)$How.Long.Organic

time <- as.data.frame(sample_data(psrel)$orgtime, sample_data(psrel)$How.Long.Organic)


org <- subset_samples(psrel,
                      Conventional.Organic == "Organic")

# clr transform phyloseq objects ----
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()

orgtransps <- org %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()


dist_mat <- phyloseq::distance(transps, method = "euclidean")
org_dist_mat <- phyloseq::distance(orgtransps, method = "euclidean")


vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female) # p = 0.009**, SIGNIFICANT
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Group) # p = 0.001***, SIGNIFICANT
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Group*phyloseq::sample_data(transps)$Male.Female) # p (Group:Male.Female) = 0.399, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Herd.Size) # p = 0.577, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Conventional.Organic) # p = 0.017*, SIGNIFICANT
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Formal.Team.Meetings.Frequency) # p = 0.987, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Cultural.Language.Barriers) # p = 0.792, not sig
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$employees) # p = 0.095, not sig
vegan::adonis2(org_dist_mat ~ phyloseq::sample_data(orgtransps)$orgtime) # p = 0.812, ns


vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Run) # p = 0.008***, SIGNIFICANT

##  PCA plot - Male Female ----

sample_data(psrel)$Run <- as.character(sample_data(psrel)$Run)

A <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Male.Female", size = 6, axes = c(2,3)) +
  scale_color_manual(values = c("#40498d", "#38aaac")) +
  stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  labs(color = "Gender") +
  ggtitle("A") + 
  labs(caption = "R2 = 0.018, F(1,70) = 1.30, P = 0.032*") +
  theme(text = element_text(size = 30)) 
C2

(A2|B2)/C2

ggsave(filename = "plots/presentation/figure2.pdf", dpi = 600, width = 18, height = 12)

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
B <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Conventional.Organic") +
  scale_color_manual(values = small_color) +
  stat_ellipse(aes(group = Conventional.Organic, color = Conventional.Organic)) + 
  theme_classic() +
  ggtitle("B") + 
  labs(caption = "R2 = 0.018, F(1, 70) = 1.31, P = 0.026*") 

B

ggsave(filename = "plots/full-run/PCA-farm-type.pdf", dpi = 600)

##  PCA plot - Age Group  ----
B <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Group", size = 6) +
  scale_color_manual(values = c("#38aaac", "#40498d")) +
  stat_ellipse(aes(group = Group, color = Group)) + 
  theme_classic() +
  ggtitle("B") + 
  labs(color = "Age Group") +
  labs(caption = "R2 = 0.035, F(1, 70) = 2.53, P = 0.001***") +
  theme(text = element_text(size = 30)) 

C3

(A3|B3)/C3

ggsave(filename = "plots/presentation/figure3.pdf", dpi = 600, width = 18, height = 12)

(A|B)

ggsave(filename = "plots/presentation/PCA-beta-multiple.pdf", dpi = 600, width = 16, height = 10)
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

# How long organic
org %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "orgtime") +
  scale_color_manual(values = color_palette) +
  stat_ellipse(aes(group = orgtime, color = orgtime)) + 
  theme_classic() +
  ggtitle("C") + 
  labs(caption = "R2 = 0.078, F(4, 67) = 1.42, P = 0.11") +
  theme(legend.position = "bottom")


## Combine Farm Plots ----

(A_farm | B_farm) / C_farm

ggsave(filename = "plots/full-run/PCA-farm-stats.pdf", dpi = 600, width = 12, height = 10)
