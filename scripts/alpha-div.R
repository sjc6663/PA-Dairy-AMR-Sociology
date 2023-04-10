# ALPHA DIVERSITY COMPARISONS
# SAB 4-4-2023

# color scheme: Viridis mako / microshades micro_cvd_blue
color_palette <- c("#BCE1FF", "#7DCCFF", "#56B4E9", "#098BD9", "#4292C6")

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(viridis)
library(BiocManager)
library(stringr)

# read in phyloseq object
ps <- readRDS("data/ransom/decontam-ps.RDS")

# create data frame with relevant metadata for comparison
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
write.csv(adiv, "tables/alpha-div.csv", row.names = FALSE)

testMF <- aov(Shannon ~ MF, data = adiv)
summary(testMF) # p = 0.00206 **

testAG <- aov(Shannon ~ AgeGroup, data = adiv)
summary(testAG) # p = 0.00128 **

testMFAG <- aov(Shannon ~ AgeGroup*MF, data = adiv)
summary(testMFAG) # p = 0.030 *

testS <- aov(Shannon ~ Size, data = adiv)
summary(testS) # p = 0.229, not significant

testOC <- aov(Shannon ~ OrgCon, data = adiv)
summary(testOC) # p = 0.0322 *

testTM <- aov(Shannon ~ TeamMeetings, data = adiv)
summary(testTM) # p = 0.76, not significant

testLB <- aov(Shannon ~ LangBarrier, data = adiv)
summary(testLB) # p = 0.0847 . not significant

testHS <- aov(Shannon ~ HerdSize, data = adiv)
summary(testHS) # p = 0.000489 ***

MF <- plot_richness(ps, x="Male.Female", measures=c("Shannon"), title = "Male vs. Female", color = "Male.Female") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
MF
ggsave(plot = MF, filename = "plots/male-vs-female.pdf", dpi = 600)

AG <- plot_richness(ps, x="Group", measures=c("Shannon"), title = "Age Group", color = "Group") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
AG
ggsave(plot = AG, filename = "plots/cows-vs-calves.pdf", dpi = 600)

MFAG <- plot_richness(ps, x="Male.Female", measures=c("Shannon"), title = "Male vs. Female by Age Group", color = "Group") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic()
  #theme(legend.position = "none")
MFAG
ggsave(plot = MFAG, filename = "plots/male-vs-female-age-group.pdf", dpi = 600)

OC <- plot_richness(ps, x="Conventional.Organic", measures=c("Shannon"), title = "Conventional vs. Organic", color = "Conventional.Organic") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
OC
ggsave(plot = OC, filename = "plots/conventional-organic.pdf", dpi = 600)

sample_data(ps)$"Formal.Team.Meetings.Frequency" <- factor(sample_data(ps)$"Formal.Team.Meetings.Frequency",
                                                           levels = c("Never", "1 or 2 times/year", "Quarterly", "Once a month", "At least twice a month"))

TM <- plot_richness(ps, x="Formal.Team.Meetings.Frequency", measures=c("Shannon"), title = "Team Meetings", color = "Formal.Team.Meetings.Frequency") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
TM
ggsave(plot = TM, filename = "plots/team-meeting-frequency.pdf", dpi = 600)

LB <- plot_richness(ps, x="Cultural.Language.Barriers", measures=c("Shannon"), title = "Language Barriers", color = "Cultural.Language.Barriers") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
LB
ggsave(plot = LB, filename = "plots/language-barriers.pdf", dpi = 600)

sample_data(ps)$"Herd.Size" <- factor(sample_data(ps)$"Herd.Size",
                                                           levels = c("0-50", "50-100", "100-150", "150-200", "200+"))

HS <- plot_richness(ps, x="Herd.Size", measures=c("Shannon"), title = "Herd Size", color = "Herd.Size") + 
  geom_boxplot() + 
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = "none")
HS
ggsave(plot = HS, filename = "plots/herd-size.pdf", dpi = 600)

## Violin Plots ----

ps.meta <- meta(ps)
ps.meta$Shannon <- phyloseq::estimate_richness(ps, measures = "Shannon")

ps.meta$'' <- alpha(ps, index = 'shannon')

p1 <- ggviolin(ps.meta, x = "Formal.Team.Meetings.Frequency", y = "Shannon$Shannon",
               add = "boxplot", fill = "Formal.Team.Meetings.Frequency", palette = color_palette) +
  scale_x_discrete(limits = c("Never", "1 or 2 times/year", "Quarterly", "Once a month", "At least twice a month")) +
  theme(legend.position = "none") + 
  stat_compare_means(comparisons = FTMF_pair, label = "p.signif",
                      symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                     )
p1

# create a list of pairwise comaprisons
FTMF <- unique(adiv$TeamMeetings) # get the variables

FTMF_pair <- combn(seq_along(FTMF), 2, simplify = FALSE, FUN = function(i)FTMF[i])


p1 <- p1 + stat_compare_means(comparisons = FTMF_pair, label = "p.format",
                              symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                              )
p1 


ggsave(filename = "plots/violin-plot-team-meetings-stats.pdf", dpi = 600)

# ------------------------------------

v_CO <- ggviolin(ps.meta, x = "Conventional.Organic", y = "Shannon$Shannon",
               add = "boxplot", fill = "Conventional.Organic", palette = color_palette) +
  theme(legend.position = "none")
v_CO

# create a list of pairwise comaprisons
CO <- unique(adiv$OrgCon) # get the variables

CO_pair <- combn(seq_along(CO), 2, simplify = FALSE, FUN = function(i)CO[i])


v_CO <- v_CO + stat_compare_means(comparisons = CO_pair, label = "p.signif",
                                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                                  )
v_CO


ggsave(filename = "plots/conventional-organic-stats.pdf", dpi = 600)

# ------------------------------------
  
v_group <- ggviolin(ps.meta, x = "Group", y = "Shannon$Shannon",
                 add = "boxplot", fill = "Group", palette = color_palette) +
  theme(legend.position = "none")
v_group

# create a list of pairwise comaprisons
group <- unique(adiv$AgeGroup) # get the variables

group_pair <- combn(seq_along(group), 2, simplify = FALSE, FUN = function(i)group[i])


v_group <- v_group + stat_compare_means(comparisons = group_pair, label = "p.signif",
                                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_group


ggsave(filename = "plots/cow-calf-stats.pdf", dpi = 600)


# ------------------------------------
  
v_MF <- ggviolin(ps.meta, x = "Male.Female", y = "Shannon$Shannon",
                      add = "boxplot", fill = "Male.Female", palette = color_palette) +
  theme(legend.position = "none")
v_MF

# create a list of pairwise comaprisons
MF <- unique(adiv$MF) # get the variables

MF_pair <- combn(seq_along(MF), 2, simplify = FALSE, FUN = function(i)MF[i])


v_MF <- v_MF + stat_compare_means(comparisons = MF_pair, label = "p.signif",
                                        symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_MF


ggsave(filename = "plots/male-female-stats.pdf", dpi = 600)


# ------------------------------------

v_MFAG <- ggviolin(ps.meta, x = "Male.Female", y = "Shannon$Shannon",
                 add = "boxplot", fill = "Group", palette = color_palette)
 # theme(legend.position = "none")
v_MFAG

# create a list of pairwise comaprisons
MFAG <- unique(adiv$AgeGroup) # get the variables

MFAG_pair <- combn(seq_along(MFAG), 2, simplify = FALSE, FUN = function(i)MFAG[i])


v_MFAG <- v_MFAG + stat_compare_means(comparisons = MFAG_pair, label = "p.signif",
                                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_MFAG


ggsave(filename = "plots/male-female-stats.pdf", dpi = 600)


# ------------------------------------

v_HS <- ggviolin(ps.meta, x = "Herd.Size", y = "Shannon$Shannon",
                   add = "boxplot", fill = "Herd.Size", palette = color_palette) +
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


ggsave(filename = "plots/herd-size-stats.pdf", dpi = 600)


# ------------------------------------

v_CLB <- ggviolin(ps.meta, x = "Cultural.Language.Barriers", y = "Shannon$Shannon",
                 add = "boxplot", fill = "Cultural.Language.Barriers", palette = color_palette) +
                theme(legend.position = "none")
v_CLB

# create a list of pairwise comaprisons
CLB <- unique(adiv$LangBarrier) # get the variables

CLB_pair <- combn(seq_along(CLB), 2, simplify = FALSE, FUN = function(i)CLB[i])


v_CLB <- v_CLB + stat_compare_means(comparisons = CLB_pair, label = "p.signif",
                                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_CLB


ggsave(filename = "plots/language-barriers-stats.pdf", dpi = 600)

# save work
save.image("data/alpha-div.RData")
