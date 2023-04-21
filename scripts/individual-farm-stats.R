## Individual Farmer AMR Results
## 4-21-2023, Stephanie Bierly

# color scheme: Viridis mako / microshades micro_cvd_blue
color_palette <- c("#BCE1FF", "#7DCCFF", "#56B4E9", "#098BD9", "#4292C6")

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(viridis)
library(BiocManager)
library(stringr)
library(microshades)
library(dplyr)

# read in phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

## Farm AF ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AF"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor") +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AF")

ggsave(filename = "plots/individual-farmers-letters/Farm-AF.pdf", dpi = 600)


## Farm AK ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AK"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor") +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AK")

ggsave(filename = "plots/individual-farmers-letters/Farm-AK.pdf", dpi = 600)

## Farm AG ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AG"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor") +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AG")

ggsave(filename = "plots/individual-farmers-letters/Farm-AG.pdf", dpi = 600)

# Farm S ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"S"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor") +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm S")

ggsave(filename = "plots/individual-farmers-letters/Farm-S.pdf", dpi = 600)

# Farm B ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"B"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor") +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm B")

ggsave(filename = "plots/individual-farmers-letters/Farm-B.pdf", dpi = 600)

# Farm AI ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AI"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor") +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AI")

ggsave(filename = "plots/individual-farmers-letters/Farm-AI.pdf", dpi = 600)

# Farm A ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bA\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor") +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm A")

ggsave(filename = "plots/individual-farmers-letters/Farm-A.pdf", dpi = 600)


# Barplot Random Playing Around ----

AF <- subset_samples(psrel,
                     Farm == "AF")

meltaf <- psmelt(AF)
#filter out rows with 0 relabun 
af2 <- meltaf[meltaf$Abundance != 0.0000000000, ]

ggplot(af2, aes(x=Broadclass, fill = Sample.Type)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette[3:5]) +
  facet_wrap(~Group, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = " ",
       title = "A") +
  theme(axis.text= element_text(size = 12,face = "bold"),
        axis.title= element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold", hjust = 0),
        legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 12, color = "black",face = "bold")
  ) +
  #guides(fill="none")+
  #scale_y_continuous(expand = c(0,0.1), limits = c(0,7500))+
  coord_flip()

psrel <- psrel %>% 
  ps_mutate(
    Samples = if_else(str_detect(Farm,"B"), true = "blue", false = "grey")
  ) 

psrel <- psrel %>% 
  ps_mutate(Samples = case_when(
    Samples == "blue" & str_detect(Group, "Calves") ~ "Calves",
    Samples == "blue" & !str_detect(Group, "Calves") ~ "Cows",
    Samples == "grey" ~ "Other Farms"
  ))



# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Samples") +
  scale_color_manual(values = c("pink", "blue", "grey")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  # theme(legend.position = "none") +
  ggtitle("Farm B")

ggsave(filename = "plots/individual-farmers-letters/Farm-AF.pdf", dpi = 600)
