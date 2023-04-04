
ps <- readRDS("data/ransom/decontam-ps.RDS")

ps <- ps0

ps0@sam_data <- sample_data(samp)

group <- c("Cows", "Calves")

ps0 <- subset_samples(
  ps0,
  Group %in% group
)

# save the decontaminated phyloseq object for downstrem analysis
saveRDS(ps0, file = "data/ransom/decontam-ps.rds")


## ---- create data frames for each breed and type ----
adiv <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "Conventional.Organic" = phyloseq::sample_data(ps)$Conventional.Organic,
)
head(adiv)

## ---- Shannon ANOVA ----
test <- aov(Shannon ~ Conventional.Organic, data = adiv)
summary(test) # p = 0.0247*

plot_richness(ps, x="Male.Female", measures=c("Shannon"), title = "Males vs Females", color = "Number.Milking.Cows") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1", "#999999")) +
  theme_classic() +
  theme(legend.position = "none")
A
ggsave(filename = "plots/alpha-diversity-mf.pdf", dpi = 600)


# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")
otu_table(psrel)

# clr transform phyloseq objects at Genus level
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()

dist_mat <- phyloseq::distance(transps, method = "euclidean")

vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Male.Female) # p = 0.072
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Conventional.Organic) # p = 0.15
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$Group) # p = 0.001***
##  PCA plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Group", shape = "Group") +
  #scale_colour_brewer(palette = "BrBG") +
  stat_ellipse(aes(group = Group, color = Group)) +
  theme_classic() +
  ggtitle("A") +
  #theme(legend.position = "none") +
  labs(caption = "")


## BAR PLOTS ----

male <- subset_samples(
  ps,
  Male.Female == "Male"
)

female <- subset_samples(
  ps,
  Male.Female == "Female"
)

df.male.otu <- as.data.frame(male@otu_table)
df.male.tax <- as.data.frame(male@tax_table)

df.male.otu <- rownames_to_column(df.male.otu, var = "Gene")

df.male <- merge.data.frame(df.male.otu, df.male.tax)
 
bmale <- ggplot(df.male, aes(x=Class, fill = Broadclass)) + 
  geom_bar()+
  theme_bw()+
  ggtitle("Males") + 
  scale_fill_viridis(discrete=TRUE, option = "A", begin = 0.1, end = 0.9) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
bmale

df.female.otu <- as.data.frame(female@otu_table)
df.female.tax <- as.data.frame(female@tax_table)

df.female.otu <- rownames_to_column(df.female.otu, var = "Gene")

df.female <- merge.data.frame(df.female.otu, df.female.tax)

bfemale <- ggplot(df.female, aes(x=Broadclass, fill = Class)) + 
  geom_bar()+
  theme_bw()+
  ggtitle("Females") + 
  scale_fill_viridis(discrete=TRUE, option = "A", begin = 0.1, end = 0.9) +
  theme(legend.position = "none")
bfemale

library(patchwork)
bmale/bfemale

df.otu <- as.data.frame(ps0@otu_table)
df.tax <- as.data.frame(ps0@tax_table)
df.otu <- rownames_to_column(df.otu, var = "Gene")
df <- merge.data.frame(df.otu, df.tax)
df2 <- t(df)
df2 <- as.data.frame(df2)

df.met <- as.data.frame(ps0@sam_data)
df2 <- rownames_to_column(df2, var = "Sample.ID")
df3 <- merge.data.frame(met, df2)

p <- ggplot(df, aes(x=Broadclass, fill = Class)) + 
  geom_bar()+
  theme_bw()+
  ggtitle("Males v. Females") + 
  scale_fill_viridis(discrete=TRUE, option = "A", begin = 0.1, end = 0.9) +
  theme(legend.position = "none")


p
