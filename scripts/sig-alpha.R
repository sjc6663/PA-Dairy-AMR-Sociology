# alpha diversity, clinically significant genes
# SAB 07/17/2023

# set seed
set.seed(81299)

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(viridis)
library(BiocManager)
library(stringr)
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
library(car)
library(misty)
library(patchwork)

color_palette <- scale_fill_viridis(option = "G")

# read in phyloseq object
ps <- readRDS("data/full-run/sig-decontam-ps.RDS")

# filter out taxa that are not clinically significant and remove them
ps <- subset_taxa(ps, 
                  sig == "TRUE")

ps

# create data frame
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "MF" = phyloseq::sample_data(ps)$Male.Female,
  "AgeGroup" = phyloseq::sample_data(ps)$Group,
  "Size" = phyloseq::sample_data(ps)$Number.Milking.Cows,
  "OrgCon" = phyloseq::sample_data(ps)$Conventional.Organic,
  "TeamMeetings" = phyloseq::sample_data(ps)$Formal.Team.Meetings.Frequency,
  "LangBarrier" = phyloseq::sample_data(ps)$Cultural.Language.Barriers,
  "HerdSize" = phyloseq::sample_data(ps)$Herd.Size
)

# test variances ----

# Male Female 
varMF <- var.test(Shannon ~ MF, data = adiv, 
                  alternative = "two.sided")
varMF # p = 0.008, SIGNIFICANT DIFFERENCE

varAG <- var.test(Shannon ~ AgeGroup, data = adiv,
                  alternative = "two.sided")
varAG # p = 9.402e-12, SIGNIFICANT DIFFERENCE

varOC <- var.test(Shannon ~ OrgCon, data = adiv,
                  alternative = "two.sided")
varOC # p = 0.037, SIGNIFICANT DIFFERENCE

varTM <- leveneTest(Shannon ~ TeamMeetings, data = adiv) # have to do levene's test because its more than 2 groups
varTM # p = 0.5, not sig

varLB <- var.test(Shannon ~ LangBarrier, data = adiv,
                  alternative = "two.sided")
varLB # p = 0.017, SIGNIFICANT DIFFERENCE

varHS <- leveneTest(Shannon ~ HerdSize, data = adiv)
varHS # p = 0.086, not sig

# Welch's t-test ----
# https://www.statology.org/welchs-anova-in-r/ 

# male female
wtestMF <- oneway.test(Shannon ~ MF, data = adiv, var.equal = FALSE)
wtestMF # p = 0.19, ns

# age group
wtestAG <- oneway.test(Shannon ~ AgeGroup, data = adiv, var.equal = FALSE)
wtestAG # p = 1.473e-07, SIGNFICANT

# age group as a subgroup of male/female
wtestMFAG <- oneway.test(Shannon ~ MF*AgeGroup, data = adiv, var.equal = FALSE)
wtestMFAG # p = 1.92e-05, SIGNIFICANT

# organic conventional
wtestOC <- oneway.test(Shannon ~ OrgCon, data = adiv, var.equal = FALSE)
wtestOC # p = 0.098, ns

# formal team meetings
wtestTM <- oneway.test(Shannon ~ TeamMeetings, data = adiv, var.equal = FALSE)
wtestTM # p = 0.67, ns

# cultural language barrier
wtestLB <- oneway.test(Shannon ~ LangBarrier, data = adiv, var.equal = FALSE)
wtestLB # p = 0.26, ns

# herd size
wtestHS <- oneway.test(Shannon ~ HerdSize, data = adiv, var.equal = FALSE)
wtestHS # p = 0.57, ns

# Violin Plots ---------------------------------------------

ps.meta <- meta(ps)
ps.meta$Shannon <- phyloseq::estimate_richness(ps, measures = "Shannon")

ps.meta$'' <- alpha(ps, index = 'shannon')

# Formal Team Meetings ------------------------------------

v_FTM <- ggviolin(ps.meta, x = "Formal.Team.Meetings.Frequency", y = "Shannon$Shannon",
                  add = "boxplot", fill = "Formal.Team.Meetings.Frequency", palette = color_palette, title = "C", ylab = "Shannon's Diversity Index", xlab = "Formal Team Meeting Frequency") +
  scale_x_discrete(limits = c("Never", "1 or 2 times/year", "Quarterly", "Once a month", "At least twice a month")) 
theme(legend.position = "none")
v_FTM


v_CO <- ggviolin(ps.meta, x = "Conventional.Organic", y = "Shannon$Shannon",
                 add = "boxplot", fill = "Conventional.Organic", palette = "RdBu", title = "CRG", ylab = "Shannon's Diversity Index", xlab = "Farm Type") +
  theme(legend.position = "none")
v_CO
