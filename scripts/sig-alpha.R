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

ps <- ps %>% 
  ps_mutate(
    employees = if_else(str_detect(Non.Family.Milkers, "0"), true = "No", false = "Yes")
  ) 

ps <- rarefy_even_depth(ps)

# create data frame
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "MF" = phyloseq::sample_data(ps)$Male.Female,
  "AgeGroup" = phyloseq::sample_data(ps)$Group,
  "Size" = phyloseq::sample_data(ps)$Number.Milking.Cows,
  "OrgCon" = phyloseq::sample_data(ps)$Conventional.Organic,
  "TeamMeetings" = phyloseq::sample_data(ps)$Formal.Team.Meetings.Frequency,
  "LangBarrier" = phyloseq::sample_data(ps)$Cultural.Language.Barriers,
  "HerdSize" = phyloseq::sample_data(ps)$Herd.Size, 
  "employees" = phyloseq::sample_data(ps)$employees
)

# test variances ----

# Male Female 
varMF <- var.test(Shannon ~ MF, data = adiv, 
                  alternative = "two.sided")
varMF # p = 0.67, ns

varAG <- var.test(Shannon ~ AgeGroup, data = adiv,
                  alternative = "two.sided")
varAG # p = 0.065, ns

varOC <- var.test(Shannon ~ OrgCon, data = adiv,
                  alternative = "two.sided")
varOC # p = 0.38, ns

varTM <- leveneTest(Shannon ~ TeamMeetings, data = adiv) # have to do levene's test because its more than 2 groups
varTM # p = 0.93, not sig

varLB <- var.test(Shannon ~ LangBarrier, data = adiv,
                  alternative = "two.sided")
varLB # p = 0.18, ns

varHS <- leveneTest(Shannon ~ HerdSize, data = adiv)
varHS # p = 0.53, not sig

# since none of the variances are significantly different, we can run one-way ANOVA

# One-Way ANOVA ----
# https://www.statology.org/welchs-anova-in-r/ 

# male female
wtestMF <- aov(Shannon ~ MF, data = adiv)
summary(wtestMF) # P = 0.00177**, SIGNIFICANT

# age group
wtestAG <- aov(Shannon ~ AgeGroup, data = adiv)
summary(wtestAG) # p = 0.00693**, SIGNIFICANT

# age group as a subgroup of male/female
wtestMFAG <- aov(Shannon ~ MF*AgeGroup, data = adiv)
summary(wtestMFAG) # MF:AgeGroup, P = 0.04*, SIGNIFICANT

# organic conventional
wtestOC <- aov(Shannon ~ OrgCon, data = adiv)
summary(wtestOC) # p = 0.018*, SIGNIFICANT

# formal team meetings
wtestTM <- aov(Shannon ~ TeamMeetings, data = adiv)
summary(wtestTM) # p = 0.83, ns

# cultural language barrier
wtestLB <- aov(Shannon ~ LangBarrier, data = adiv)
summary(wtestLB) # p = 0.063, ns

# herd size
wtestHS <- aov(Shannon ~ HerdSize, data = adiv)
summary(wtestHS) # p = 0.0128*, SIGNIFICANT

# non-family employees
wtestNFE <- aov(Shannon ~ employees, data = adiv)
summary(wtestNFE) # p = 0.0102*, SIGNIFICANT

# Violin Plots ---------------------------------------------

sample_data(ps)$Herd.Size <- factor(sample_data(ps)$Herd.Size,
                                       levels = c("0-50", "50-100", "100-150", "150-200", "200+"))

ps.meta <- meta(ps)
ps.meta$Shannon <- phyloseq::estimate_richness(ps, measures = "Shannon")

ps.meta$'' <- alpha(ps, index = 'shannon')

# Male v. Female ------------------------------------

A <- ggviolin(ps.meta, x = "Male.Female", y = "Shannon$Shannon",
              add = "boxplot", fill = "Male.Female", palette = c("#367aa1", "#def4e5"), title = "A", ylab = "Shannon's Diversity Index", xlab = "Primary Farm Operator Gender") +
  theme(legend.position = "none")
A

# create a list of pairwise comaprisons
A2 <- unique(adiv$MF) # get the variables

A_pair <- combn(seq_along(A2), 2, simplify = FALSE, FUN = function(i)A2[i])


v_A <- A + stat_compare_means(comparisons = A_pair, label = "p.signif",
                              symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_A

# Cow v. Calves

B <- ggviolin(ps.meta, x = "Group", y = "Shannon$Shannon",
              add = "boxplot", fill = "Group", palette = c("#367aa1", "#def4e5"), title = "B", ylab = "Shannon's Diversity Index", xlab = "Age Group") +
  theme(legend.position = "none")
B

# create a list of pairwise comaprisons
group <- unique(adiv$AgeGroup) # get the variables

group_pair <- combn(seq_along(group), 2, simplify = FALSE, FUN = function(i)group[i])


v_B <- B + stat_compare_means(comparisons = group_pair, label = "p.signif",
                              symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_B

# Organic v. Conventional

v_C <- ggviolin(ps.meta, x = "Conventional.Organic", y = "Shannon$Shannon",
                add = "boxplot", fill = "Conventional.Organic", palette = c("#367aa1", "#def4e5"), title = "C", ylab = "Shannon's Diversity Index", xlab = "Farm Type") +
  theme(legend.position = "none")


v_C

# create a list of pairwise comaprisons
CO <- unique(adiv$OrgCon) # get the variables

CO_pair <- combn(seq_along(CO), 2, simplify = FALSE, FUN = function(i)CO[i])


v_C <- v_C + stat_compare_means(comparisons = CO_pair, label = "p.signif",
                                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_C

# Non-Family Employees

v_D <- ggviolin(ps.meta, x = "employees", y = "Shannon$Shannon",
                add = "boxplot", fill = "employees", palette = c("#367aa1", "#def4e5"), title = "D", ylab = "Shannon's Diversity Index", xlab = "Non-Family Employees")  +
  theme(legend.position = "none")
v_D

# create a list of pairwise comaprisons
NFE <- unique(adiv$employees) # get the variables

NFE_pair <- combn(seq_along(NFE), 2, simplify = FALSE, FUN = function(i)NFE[i])


v_D <- v_D + stat_compare_means(comparisons = NFE_pair, label = "p.signif",
                                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_D

# Herd Size

v_E <- ggviolin(ps.meta, x = "Herd.Size", y = "Shannon$Shannon",
                 add = "boxplot", fill = "Herd.Size", palette = c("#367aa1", "#348fa7", "#40b7ad", "#8ad9b1", "#def4e5"), title = "E", ylab = "Shannon's Diversity Index", xlab = "Herd Size") +
  scale_x_discrete(limits = c("0-50", "50-100", "100-150", "150-200", "200+"))
# theme(legend.position = "none")
v_E

# create a list of pairwise comaprisons
HS <- unique(adiv$HerdSize) # get the variables

HS_pair <- combn(seq_along(HS), 2, simplify = FALSE, FUN = function(i)HS[i])


v_E <- v_E + stat_compare_means(comparisons = HS_pair, label = "p.signif", hide.ns = TRUE, 
                                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_E

(v_A/v_B/v_C)|(v_D/v_E)

ggsave(filename = "plots/presentation/alpha-crg.pdf", dpi = 600, width = 15, height = 20)
