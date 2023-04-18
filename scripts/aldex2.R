# ALDEx2

# load packages
require(tidyverse)
require(phyloseq)
library(ALDEx2)
library(microViz)
library(ggpubr)
library(ggrepel)

# load phyloseq of counts
# load("ps-decontam-filtered-counts.RData")


# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps)), phyloseq::sample_data(ps)$Group, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "Gene_")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# make a table of significant corrected p-values
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "Gene_") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(Gene_, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
sig_aldex2 <- left_join(sig_aldex2, taxa_info)
sig_aldex2
write.csv(sig_aldex2, file = "tables/aldex-cow-calf.csv")

ggplot(data = sig_aldex2, aes(x = Broadclass, y = effect)) +
  geom_point()
  
ggplot(data=sig_aldex2, aes(x=Broadclass, y=(effect), col=Broadclass, label=Gene)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = color_palette) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black")

ggsave(filename = "plots/aldex2-cow-calf.pdf", dpi = 600)
