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

ps2 <- rarefy_even_depth(ps)

# transform to relative abundance
psrel <- microbiome::transform(ps2, "compositional")

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
  ord_plot(color = "pointcolor", size = 6) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm A")

ggsave(filename = "plots/individual-farmers-letters/Farm-A.pdf", dpi = 600)

# Farm AA ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bAA\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 6) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AA")

ggsave(filename = "plots/individual-farmers-letters/Farm-AA.pdf", dpi = 600)


# Farm AE ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bAE\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 6) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AE")

ggsave(filename = "plots/individual-farmers-letters/Farm-AE.pdf", dpi = 600)

# Farm AD ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bAD\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 6) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AD")

ggsave(filename = "plots/individual-farmers-letters/Farm-AD.pdf", dpi = 600)

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
  ord_plot(color = "pointcolor", size = 6) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AF")

ggsave(filename = "plots/individual-farmers-letters/Farm-AF.pdf", dpi = 600)


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
  ord_plot(color = "pointcolor", size = 6) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AG")

ggsave(filename = "plots/individual-farmers-letters/Farm-AG.pdf", dpi = 600)


## Farm AH ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AH"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 3) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AH")

ggsave(filename = "plots/individual-farmers-letters/Farm-AH.pdf", dpi = 600)

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
  ord_plot(color = "pointcolor", size = 6) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AI")

ggsave(filename = "plots/individual-farmers-letters/Farm-AI.pdf", dpi = 600)

## Farm AJ ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AJ"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AJ")

ggsave(filename = "plots/individual-farmers-letters/Farm-AJ.pdf", dpi = 600)

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
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AK")

ggsave(filename = "plots/individual-farmers-letters/Farm-AK.pdf", dpi = 600)

## Farm AL ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AL"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AL")

ggsave(filename = "plots/individual-farmers-letters/Farm-AL.pdf", dpi = 600)

## Farm AM ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AM"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AM")

ggsave(filename = "plots/individual-farmers-letters/Farm-AM.pdf", dpi = 600)

## Farm AN ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AN"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AN")

ggsave(filename = "plots/individual-farmers-letters/Farm-AN.pdf", dpi = 600)

## Farm AO ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AO"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AO")

ggsave(filename = "plots/individual-farmers-letters/Farm-AO.pdf", dpi = 600)

## Farm AP ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AP"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AP")

ggsave(filename = "plots/individual-farmers-letters/Farm-AP.pdf", dpi = 600)

## Farm AQ ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AQ"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AQ")

ggsave(filename = "plots/individual-farmers-letters/Farm-AQ.pdf", dpi = 600)

## Farm AR ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AR"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AR")

ggsave(filename = "plots/individual-farmers-letters/Farm-AR.pdf", dpi = 600)

## Farm AS ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"AS"), true = "blue", false = "grey")
  ) 

# AF Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm AS")

ggsave(filename = "plots/individual-farmers-letters/Farm-AS.pdf", dpi = 600)

## Farm B ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bB\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm B")

ggsave(filename = "plots/individual-farmers-letters/Farm-B.pdf", dpi = 600)

## Farm E ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bE\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm E")

ggsave(filename = "plots/individual-farmers-letters/Farm-E.pdf", dpi = 600)

## Farm F ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bF\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm F")

ggsave(filename = "plots/individual-farmers-letters/Farm-F.pdf", dpi = 600)

## Farm H ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bH\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm H")

ggsave(filename = "plots/individual-farmers-letters/Farm-H.pdf", dpi = 600)

## Farm K ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bK\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm K")

ggsave(filename = "plots/individual-farmers-letters/Farm-K.pdf", dpi = 600)

## Farm L ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bL\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm L")

ggsave(filename = "plots/individual-farmers-letters/Farm-L.pdf", dpi = 600)

## Farm O ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bO\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm O")

ggsave(filename = "plots/individual-farmers-letters/Farm-O.pdf", dpi = 600)

## Farm P ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bP\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm P")

ggsave(filename = "plots/individual-farmers-letters/Farm-I.pdf", dpi = 600)

## Farm Q ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bQ\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm Q")

ggsave(filename = "plots/individual-farmers-letters/Farm-Q.pdf", dpi = 600)

## Farm W ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bW\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm W")

ggsave(filename = "plots/individual-farmers-letters/Farm-W.pdf", dpi = 600)

## Farm X ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bX\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm X")

ggsave(filename = "plots/individual-farmers-letters/Farm-X.pdf", dpi = 600)

## Farm Y ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bY\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm Y")

ggsave(filename = "plots/individual-farmers-letters/Farm-Y.pdf", dpi = 600)

## Farm N ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bC\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm C")

ggsave(filename = "plots/individual-farmers-letters/Farm-C.pdf", dpi = 600)

## Farm I ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bI\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm I")

ggsave(filename = "plots/individual-farmers-letters/Farm-I.pdf", dpi = 600)

## Farm J ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bJ\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm J")

ggsave(filename = "plots/individual-farmers-letters/Farm-J.pdf", dpi = 600)

## Farm R ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bR\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm R")

ggsave(filename = "plots/individual-farmers-letters/Farm-R.pdf", dpi = 600)

## Farm S ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bS\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm S")

ggsave(filename = "plots/individual-farmers-letters/Farm-S.pdf", dpi = 600)

# Farm V ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bV\\b"), true = "blue", false = "grey")
  ) 

# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 4) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm V")

ggsave(filename = "plots/individual-farmers-letters/Farm-V.pdf", dpi = 600)

# Farm P ----
psrel <- psrel %>% 
  ps_mutate(
    pointcolor = if_else(str_detect(Farm,"\\bP\\b"), true = "blue", false = "grey")
  ) 


# Plot
psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "pointcolor", size = 6) +
  scale_color_manual(values = c("blue", "gray")) +
  # stat_ellipse(aes(group = Male.Female, color = Male.Female)) + 
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Farm P")

ggsave(filename = "plots/individual-farmers-letters/Farm-P.pdf", dpi = 600)
