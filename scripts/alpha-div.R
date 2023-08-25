# ALPHA DIVERSITY COMPARISONS
# SAB 4-4-2023

# setup ----
# color scheme: Viridis mako / microshades micro_cvd_blue

color_palette <- c("#8ad9b1", "#40b7ad", "#348fa7", "#37659e", "#423d7b")
small_color <- c("#40b7ad", "#423d7b")

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

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

ps2 <- rarefy_even_depth(ps)

ps2 <- ps2 %>% 
  ps_mutate(
    employees = if_else(str_detect(Non.Family.Milkers, "0"), true = "No", false = "Yes")
  ) 

# create data frame ----
# create data frame with relevant metadata for comparison
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps2, measures = "Shannon"),
  "MF" = phyloseq::sample_data(ps2)$Male.Female,
  "AgeGroup" = phyloseq::sample_data(ps2)$Group,
  "Size" = phyloseq::sample_data(ps2)$Number.Milking.Cows,
  "OrgCon" = phyloseq::sample_data(ps2)$Conventional.Organic,
  "TeamMeetings" = phyloseq::sample_data(ps2)$Formal.Team.Meetings.Frequency,
  "LangBarrier" = phyloseq::sample_data(ps2)$Cultural.Language.Barriers,
  "HerdSize" = phyloseq::sample_data(ps2)$Herd.Size,
  "employees" = phyloseq::sample_data(ps2)$employees
)
write.csv(adiv, "tables/full-run_alpha-div.csv", row.names = FALSE)

# test variances ----

# Male Female 
varMF <- var.test(Shannon ~ MF, data = adiv, 
         alternative = "two.sided")
varMF # p = 0.05, not sig

varAG <- var.test(Shannon ~ AgeGroup, data = adiv,
                  alternative = "two.sided")
varAG # p = 0.003***

varOC <- var.test(Shannon ~ OrgCon, data = adiv,
                  alternative = "two.sided")
varOC # p = 0.11, not sig

varTM <- leveneTest(Shannon ~ TeamMeetings, data = adiv) # have to do levene's test because its more than 2 groups
varTM # p = 0.51, not sig

varLB <- var.test(Shannon ~ LangBarrier, data = adiv,
                  alternative = "two.sided")
varLB # p = 0.63, not sig

varHS <- leveneTest(Shannon ~ HerdSize, data = adiv)
varHS # p = 0.37, not sig

varNFE <- var.test(Shannon ~ employees, data = adiv,
                   alternative = "two.sided")
varNFE # p = 0.78, not sig
# none of our groups have significant variance differences, so we can run standard t-tests.  

# Welch's t-test ----
# https://www.statology.org/welchs-anova-in-r/ 

# male female
wtestMF <- t.test(Shannon ~ MF, data = adiv, var.equal = FALSE)
wtestMF # p = 0.0009**, SIGNIFICANT

# age group
wtestAG <- t.test(Shannon ~ AgeGroup, data = adiv, var.equal = FALSE)
wtestAG # p = 0.00043***, SIGNFICANT

# age group as a subgroup of male/female
wtestMFAG <- test.welch(Shannon ~ AgeGroup*MF, data = adiv, var.equal = FALSE)
wtestMFAG # p = 1.49e-07, SIGNIFICANT

# organic conventional
wtestOC <- t.test(Shannon ~ OrgCon, data = adiv, var.equal = FALSE)
wtestOC # p = 0.001***, SIGNIFICANT

# formal team meetings
wtestTM <- test.welch(Shannon ~ TeamMeetings, data = adiv, var.equal = FALSE)
wtestTM # p = 0.90, ns

# cultural language barrier
wtestLB <- t.test(Shannon ~ LangBarrier, data = adiv, var.equal = FALSE)
wtestLB # p = 0.02*, SIGNIFICANT

# herd size
wtestHS <- test.welch(Shannon ~ HerdSize, data = adiv, var.equal = FALSE)
awtestHS # p = 0.12, ns

# non-family employees
wtestNFE <- t.test(Shannon ~ employees, data = adiv, var.equal = FALSE)
wtestNFE # p = 0.004***, SIGNIFICANT

# p = 0.000111

# boxplots ----
MF <- plot_richness(ps, x="Male.Female", measures=c("Shannon"), title = "Male vs. Female", color = "Male.Female") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
MF
ggsave(plot = MF, filename = "plots/full-run/male-vs-female.pdf", dpi = 600)

AG <- plot_richness(ps, x="Group", measures=c("Shannon"), title = "Age Group", color = "Group") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
AG
ggsave(plot = AG, filename = "plots/full-run/cows-vs-calves.pdf", dpi = 600)

MFAG <- plot_richness(ps, x="Male.Female", measures=c("Shannon"), title = "Male vs. Female by Age Group", color = "Group") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette[1:5]) +
  theme_classic()
  #theme(legend.position = "none")
MFAG
ggsave(plot = MFAG, filename = "plots/full-run/male-vs-female-age-group.pdf", dpi = 600)

OC <- plot_richness(ps, x="Conventional.Organic", measures=c("Shannon"), title = "Conventional vs. Organic", color = "Conventional.Organic") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
OC
ggsave(plot = OC, filename = "plots/full-run/conventional-organic.pdf", dpi = 600)

sample_data(ps)$"Formal.Team.Meetings.Frequency" <- factor(sample_data(ps)$"Formal.Team.Meetings.Frequency",
                                                           levels = c("Never", "1 or 2 times/year", "Quarterly", "Once a month", "At least twice a month"))

TM <- plot_richness(ps, x="Formal.Team.Meetings.Frequency", measures=c("Shannon"), title = "Team Meetings", color = "Formal.Team.Meetings.Frequency") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
TM
ggsave(plot = TM, filename = "plots/full-run/team-meeting-frequency.pdf", dpi = 600)

LB <- plot_richness(ps, x="Cultural.Language.Barriers", measures=c("Shannon"), title = "Language Barriers", color = "Cultural.Language.Barriers") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
LB
ggsave(plot = LB, filename = "plots//full-run/language-barriers.pdf", dpi = 600)

sample_data(ps)$"Herd.Size" <- factor(sample_data(ps)$"Herd.Size",
                                                           levels = c("0-50", "50-100", "100-150", "150-200", "200+"))

HS <- plot_richness(ps, x="Herd.Size", measures=c("Shannon"), title = "Herd Size", color = "Herd.Size") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
HS
ggsave(plot = HS, filename = "plots/full-run/herd-size.pdf", dpi = 600)

## Violin Plots ----
# https://microbiome.github.io/tutorials/PlotDiversity.html 
ps.meta <- meta(ps2)
ps.meta$Shannon <- phyloseq::estimate_richness(ps2, measures = "Shannon")

ps.meta$'' <- alpha(ps2, index = 'shannon')

# Formal Team Meetings ------------------------------------

ggviolin(ps.meta, x = "Formal.Team.Meetings.Frequency", y = "Shannon$Shannon",
                 add = "boxplot", fill = "Formal.Team.Meetings.Frequency", palette = color_palette, title = "C", ylab = "Shannon's Diversity Index", xlab = "Formal Team Meeting Frequency") +
  scale_x_discrete(limits = c("Never", "1 or 2 times/year", "Quarterly", "Once a month", "At least twice a month")) 
theme(legend.position = "none")
v_FTM


# create a list of pairwise comaprisons
FTMF <- unique(adiv$TeamMeetings) # get the variables

FTMF_pair <- combn(seq_along(FTMF), 2, simplify = FALSE, FUN = function(i)FTMF[i])

ggviolin(ps.meta, x = "Formal.Team.Meetings.Frequency", y = "Shannon$Shannon",
               add = "boxplot", fill = "Formal.Team.Meetings.Frequency", palette = color_palette) +
  scale_x_discrete(limits = c("Never", "1 or 2 times/year", "Quarterly", "Once a month", "At least twice a month")) +
  theme(legend.position = "none") + 
  stat_compare_means(comparisons = FTMF_pair, label = "p.signif",
                      symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                     )
p1



ggsave(filename = "plots/full-run/violin-plot-team-meetings-stats.pdf", dpi = 600)

# Conventional / Organic ------------------------------------

v_C <- ggviolin(ps.meta, x = "Conventional.Organic", y = "Shannon$Shannon",
               add = "boxplot", fill = "Conventional.Organic", palette = c("#38aaac", "#40498d"), title = "C", ylab = "Shannon's Diversity Index", xlab = "Farm Type") +
  theme(legend.position = "none")


v_C

# create a list of pairwise comaprisons
CO <- unique(adiv$OrgCon) # get the variables

CO_pair <- combn(seq_along(CO), 2, simplify = FALSE, FUN = function(i)CO[i])


v_C <- v_C + stat_compare_means(comparisons = CO_pair, label = "p.signif",
                                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                                  )
v_C


ggsave(filename = "plots/presentation/conventional-organic-stats.pdf", dpi = 600)

# Age Group ------------------------------------
  
B3 <- ggviolin(ps.meta, x = "Group", y = "Shannon$Shannon",
                 add = "boxplot", fill = "Group", palette = c("#38aaac", "#40498d"), title = "B", ylab = "Shannon's Diversity Index", xlab = "Age Group") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 30), axis.title = element_text(size = 20))
B3

# create a list of pairwise comaprisons
group <- unique(adiv$AgeGroup) # get the variables

group_pair <- combn(seq_along(group), 2, simplify = FALSE, FUN = function(i)group[i])


v_B <- B + stat_compare_means(comparisons = group_pair, label = "p.signif",
                                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_B


ggsave(filename = "plots/full-run/cow-calf-stats.pdf", dpi = 600)


# Male / Female ------------------------------------
  
B2 <- ggviolin(ps.meta, x = "Male.Female", y = "Shannon$Shannon",
                      add = "boxplot", fill = "Male.Female", palette = c("#38aaac", "#40498d"), title = "B", ylab = "Shannon's Diversity Index", xlab = "Farm Manager Gender") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 30), axis.title = element_text(size = 20))
B2

# create a list of pairwise comaprisons
A2 <- unique(adiv$MF) # get the variables

A_pair <- combn(seq_along(A2), 2, simplify = FALSE, FUN = function(i)A2[i])


B2 <- B2 + stat_compare_means(comparisons = A_pair, label = "p.signif",
                                        symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
B2


ggsave(filename = "plots/full-run/male-female-stats.pdf", dpi = 600)


# Male / Female Age Groups (PROBLEMS WITH THIS PLOT) ------------------------------------

v_MFAG <- ggviolin(ps.meta, x = "Male.Female", y = "Shannon$Shannon",
                 add = "boxplot", fill = "Group", palette = c("#0070FF", "#FFC55A"), title = "B", ylab = "Shannon's Diversity Index", xlab = "Primary Farm Operator Gender")
 # theme(legend.position = "none")
v_MFAG

ggsave(v_MFAG, filename = "plots/presentation/age-gender.pdf", dpi = 600)

# create a list of pairwise comaprisons
MFAG <- unique(adiv$AgeGroup) # get the variables

MFAG_pair <- combn(seq_along(MFAG), 2, simplify = FALSE, FUN = function(i)MFAG[i])


v_MFAG <- v_MFAG + stat_compare_means(comparisons = MFAG_pair, label = "p.signif",
                                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_MFAG


ggsave(filename = "plots/male-female-stats.pdf", dpi = 600)


# Herd Size ------------------------------------

v_HS <- ggviolin(ps.meta, x = "Herd.Size", y = "Shannon$Shannon",
                   add = "boxplot", fill = "Herd.Size", palette = color_palette, title = "C", ylab = "Shannon's Diversity Index", xlab = "Herd Size") +
  scale_x_discrete(limits = c("0-50", "50-100", "100-150", "150-200", "200+"))
# theme(legend.position = "none")
v_HS

# create a list of pairwise comaprisons
HS <- unique(adiv$HerdSize) # get the variables

HS_pair <- combn(seq_along(HS), 2, simplify = FALSE, FUN = function(i)HS[i])


v_HS <- v_HS + stat_compare_means(comparisons = HS_pair, label = "p.signif",
                                      symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_HS


ggsave(filename = "plots/full-run/herd-size-stats.pdf", dpi = 600)


# Language Barriers ------------------------------------

B5 <- ggviolin(ps.meta, x = "Cultural.Language.Barriers", y = "Shannon$Shannon",
                 add = "boxplot", fill = "Cultural.Language.Barriers", palette = c("#38aaac", "#40498d"), title = "B", ylab = "Shannon's Diversity Index", xlab = "Language Barriers") +
                theme(legend.position = "none") +
  theme(text = element_text(size = 30), axis.title = element_text(size = 20))
B5

A5|B5

ggsave(filename = "plots/presentation/figure5.pdf", dpi = 600, width = 10, height = 6)

# create a list of pairwise comaprisons
CLB <- unique(adiv$LangBarrier) # get the variables

CLB_pair <- combn(seq_along(CLB), 2, simplify = FALSE, FUN = function(i)CLB[i])


v_D <- v_D + stat_compare_means(comparisons = CLB_pair, label = "p.signif",
                                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_D


ggsave(filename = "plots/full-run/language-barriers-stats.pdf", dpi = 600)

# non-family employees --------------------
ps <- ps %>% 
  ps_mutate(
    employees = if_else(str_detect(Non.Family.Milkers, "0"), true = "No", false = "Yes")
  ) 

sample_data(ps)

B4 <- ggviolin(ps.meta, x = "employees", y = "Shannon$Shannon",
                   add = "boxplot", fill = "employees", palette = c("#38aaac", "#40498d"), title = "B", ylab = "Shannon's Diversity Index", xlab = "Non-Family Employees")  +
  theme(legend.position = "none") +
  theme(text = element_text(size = 30), axis.title = element_text(size = 20))
B4

(A4|B4)

ggsave(filename = "plots/presentation/figure4.pdf", dpi = 600, width = 10, height = 6)
# create a list of pairwise comaprisons
NFE <- unique(adiv$employees) # get the variables

NFE_pair <- combn(seq_along(NFE), 2, simplify = FALSE, FUN = function(i)NFE[i])


v_E <- v_E + stat_compare_means(comparisons = NFE_pair, label = "p.signif",
                                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_E

ggsave(v_emp, filename = "plots/presentation/nonfamemp.pdf", dpi = 600)

(v_A/v_B/v_C)|(v_D/v_E)

ggsave(filename = "plots/presentation/all-alpha-div.pdf", dpi = 600, width = 12, height = 20)

## COMBINED PLOTS ----
# Combined Plots Showing CO, AG, HS 
(v_CO | v_group) / v_HS

ggsave(filename = "plots/full-run/combined-farm-details.pdf", dpi = 600, width = 12, height = 18)

# Combined Plots Showing FTM, MF, LB, MFAG

(v_MF | v_MFAG ) / (v_FTM | v_CLB)

ggsave(filename = "plots/full-run/combined-soc-details.pdf", dpi = 600, width = 18, height = 18)

# save work ----
save.image("data/alpha-div.RData")


# dotplot

ps.meta$Shannon <- as.factor(ps.meta$Shannon)

ggplot(ps.meta, aes(x=Batch, y=Conventional.Organic, fill=Run)) +
  geom_dotplot(binaxis='x', stackdir='center')

sample_data(ps2)$Run <- as.character(sample_data(ps2)$Run)
sample_data(ps)$Batch <- as.character(sample_data(ps)$Batch)

plot_richness(ps2, x = "Run", measures = c("Shannon"), color = "Male.Female", shape = "Group") +
  theme_classic()

