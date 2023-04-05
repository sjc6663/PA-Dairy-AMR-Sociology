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
  theme(legend.position = "none")
p1

ggsave(filename = "plots/violin-plot-team-meetings.pdf", dpi = 600)


# save work
save.image("data/alpha-div.RData")
