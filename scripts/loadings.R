# Biplots with Loadings - PCA
# May 9, 2023 - SAB
# https://joey711.github.io/phyloseq/plot_ordination-examples.html 

ps <- readRDS("data/full-run/decontam-ps.RDS")

psrel <- microbiome::transform(ps, "compositional")
library(stats)
library(phyloseq)
library(microViz)
library(ggplot2)

##  PCA plot - Male Female ----
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Male.Female") +
  scale_color_manual(values = colors) +
  stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  ggtitle("A") + 
  labs(caption = "R2 = 0.016, F(1,70) = 1.07, P = 0.30")

rel.ord <- ordinate(psrel, "NMDS", "bray")

p1 = plot_ordination(psrel, rel.ord, type="taxa", color="Broadclass", title="broadclass")

p1 + facet_wrap(~Broadclass, 3)

p2 <- plot_ordination(psrel, rel.ord, type = "samples", color = "Farm", shape = "Male.Female", label = "Farm") + theme(legend.position = "none")
p2

p2 + geom_polygon(aes(fill="Group")) + geom_point(size = 5) + ggtitle("Farms")

p3 <- plot_ordination(psrel, rel.ord, type = "split", color = "Broadclass", shape = "Group", label = "Farm", title = "split")
p3

dist.bc <- distance(psrel, "bray")
kin.mds <- metaMDS(dist.bc, trace = 0)
## Vector fitting
ef <- envfit(kin.mds, sample_data(psrel)) ## Plot only most significant variables plot(kin.mds)
plot.new()
plot(ef)


psdf <- psmelt(psrel)






ps2 <- psrel %>%
  ps_mutate(
    Male.Female = as.numeric(Male.Female == "male"),
    Group = as.numeric(Group == "cows")
  )


ps3 <- ps2 %>% tax_transform("clr")
 
  
ps4 <- ps3 %>% ord_calc(
    constraints = c("Male.Female", "Group"),
    # method = "RDA", # Note: you can specify RDA explicitly, and it is good practice to do so, but microViz can guess automatically that you want an RDA here (helpful if you don't remember the name?)
    scale_cc = FALSE # doesn't make a difference
  )

  ord_plot(ps4)

  
  psrel %>% 
    # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
    tax_transform(trans = "clr") %>%
    ord_calc(method = "PCA") %>% 
    ord_plot(constraints = "Male.Female", color = "Male.Female") +
    scale_color_manual(values = colors) +
    stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
    theme_classic() +
    ggtitle("A") + 
    labs(caption = "R2 = 0.016, F(1,70) = 1.07, P = 0.30")
  
  
sample_data(psrel)$Male.Female[sample_data(psrel)$Male.Female == "Male"] <- 1
sample_data(psrel)$Male.Female[sample_data(psrel)$Male.Female == "Female"] <- 2

sample_data(psrel)$Male.Female <- as.numeric(sample_data(psrel)$Male.Female)
is.numeric(sample_data(psrel)$Male.Female)

sample_data(psrel)$Group[sample_data(psrel)$Group == "Cows"] <- 1
sample_data(psrel)$Group[sample_data(psrel)$Group == "Calves"] <- 2

sample_data(psrel)$Group <- as.numeric(sample_data(psrel)$Group)
is.numeric(sample_data(psrel)$Group)

sample_data(psrel)$Conventional.Organic[sample_data(psrel)$Conventional.Organic == "Conventional"] <- 1
sample_data(psrel)$Conventional.Organic[sample_data(psrel)$Conventional.Organic == "Organic"] <- 2

sample_data(psrel)$Conventional.Organic <- as.numeric(sample_data(psrel)$Conventional.Organic)
is.numeric(sample_data(psrel)$Conventional.Organic)
  
sample_data(psrel)$Herd.Size[sample_data(psrel)$Herd.Size == "0-50"] <- 1
sample_data(psrel)$Herd.Size[sample_data(psrel)$Herd.Size == "50-100"] <- 2
sample_data(psrel)$Herd.Size[sample_data(psrel)$Herd.Size == "100-150"] <- 3
sample_data(psrel)$Herd.Size[sample_data(psrel)$Herd.Size == "150-200"] <- 4
sample_data(psrel)$Herd.Size[sample_data(psrel)$Herd.Size == "200+"] <- 5

sample_data(psrel)$Herd.Size <- as.numeric(sample_data(psrel)$Herd.Size)
is.numeric(sample_data(psrel)$Herd.Size)

sample_data(psrel)$Formal.Team.Meetings.Frequency[sample_data(psrel)$Formal.Team.Meetings.Frequency == "Never"] <- 1
sample_data(psrel)$Formal.Team.Meetings.Frequency[sample_data(psrel)$Formal.Team.Meetings.Frequency == "1 or 2 times/year"] <- 2
sample_data(psrel)$Formal.Team.Meetings.Frequency[sample_data(psrel)$Formal.Team.Meetings.Frequency == "Quarterly"] <- 3
sample_data(psrel)$Formal.Team.Meetings.Frequency[sample_data(psrel)$Formal.Team.Meetings.Frequency == "Once a month"] <- 4
sample_data(psrel)$Formal.Team.Meetings.Frequency[sample_data(psrel)$Formal.Team.Meetings.Frequency == "At least twice a month"] <- 5

sample_data(psrel)$Formal.Team.Meetings.Frequency <- as.numeric(sample_data(psrel)$Formal.Team.Meetings.Frequency)
is.numeric(sample_data(psrel)$Formal.Team.Meetings.Frequency)

sample_data(psrel)$Cultural.Language.Barriers[sample_data(psrel)$Cultural.Language.Barriers == "Yes"] <- 1
sample_data(psrel)$Cultural.Language.Barriers[sample_data(psrel)$Cultural.Language.Barriers == "No"] <- 2

sample_data(psrel)$Cultural.Language.Barriers <- as.numeric(sample_data(psrel)$Cultural.Language.Barriers)
is.numeric(sample_data(psrel)$Cultural.Language.Barriers)

# alternatively, constrain variation on weight and female
  constrained_aitchison_rda <- psrel %>%
    tax_transform("clr") %>%
    ord_calc(constraints = c("Male.Female", "Group", "Conventional.Organic", "Herd.Size", "Formal.Team.Meetings.Frequency", "Cultural.Language.Barriers"), 
             method = "RDA") # constraints --> RDA
  
  constrained_aitchison_rda %>%
    ord_plot(colour = "Male.Female", constraint_vec_length = 2)
  
ggsave(filename = "plots/biplot.pdf", dpi = 600, width = 40, height = 40)

