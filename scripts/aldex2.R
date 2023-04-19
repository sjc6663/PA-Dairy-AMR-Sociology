# ALDEx2
# Stephanie Bierly 4/19/2023

# setup ----
# load packages
require(tidyverse)
require(phyloseq)
library(ALDEx2)
library(microViz)
library(ggpubr)
library(ggrepel)

# load phyloseq of counts
ps <- readRDS("data/full-run/decontam-ps.RDS")

# color scheme: Viridis mako / microshades micro_cvd_blue
color_palette <- c("#BCE1FF", "#7DCCFF", "#56B4E9", "#098BD9", "#4292C6")

## Age Group ----
# run AlDEx2 function
aldex2_da_G <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps)), phyloseq::sample_data(ps)$Group, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
gsig_aldex2 <- aldex2_da_G %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "Gene_")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da_G, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# make a table of significant corrected p-values
gsig_aldex2 <- aldex2_da_G %>%
  rownames_to_column(var = "Gene_") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(Gene_, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
gsig_aldex2 <- left_join(gsig_aldex2, taxa_info)
gsig_aldex2
write.csv(gsig_aldex2, file = "tables/aldex-cow-calf.csv")

ggplot(data = gsig_aldex2, aes(x = Broadclass, y = effect)) +
  geom_point()
  
ggplot(data=gsig_aldex2, aes(x=Broadclass, y=(effect), col=Broadclass, label=Gene)) +
  geom_point() + 
  geom_jitter() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = color_palette) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black") +
  theme(legend.position = "none")

ggsave(filename = "plots/full-run/aldex2-cow-calf.pdf", dpi = 600)


## Male Female ----
# run AlDEx2 function
aldex2_da_MF <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps)), phyloseq::sample_data(ps)$Male.Female, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
mfsig_aldex2 <- aldex2_da_MF %>%
  filter(wi.eBH < 0.05)

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da_MF, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# make a table of significant corrected p-values
mfsig_aldex2 <- aldex2_da_MF %>%
  rownames_to_column(var = "Gene_") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(Gene_, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
mfsig_aldex2 <- left_join(mfsig_aldex2, taxa_info)
mfsig_aldex2
write.csv(sig_aldex2, file = "tables/aldex-cow-calf.csv")

ggplot(data=mfsig_aldex2, aes(x=Broadclass, y=(effect), col=Broadclass, label=Gene)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = color_palette) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black") +
  theme(legend.position = "none")

## Organic Conventional ----
# run AlDEx2 function
aldex2_da_OC <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps)), phyloseq::sample_data(ps)$Conventional.Organic, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
ocsig_aldex2 <- aldex2_da_OC %>%
  filter(wi.eBH < 0.05)

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da_OC, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# make a table of significant corrected p-values
ocsig_aldex2 <- aldex2_da_OC %>%
  rownames_to_column(var = "Gene_") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(Gene_, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
ocsig_aldex2 <- left_join(ocsig_aldex2, taxa_info)
ocsig_aldex2
write.csv(ocsig_aldex2, file = "tables/aldex-farm-type.csv")

ggplot(data = ocsig_aldex2, aes(x = Broadclass, y = effect)) +
  geom_point()

ggplot(data=ocsig_aldex2, aes(x=Broadclass, y=(effect), col=Broadclass, label=Gene)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = color_palette[4]) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black") +
  theme(legend.position = "none")

ggsave(filename = "plots/full-run/aldex2-farm-type.pdf", dpi = 600)




## Cultural Language Barriers ----
# run AlDEx2 function
aldex2_da_LB <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps)), phyloseq::sample_data(ps)$Cultural.Language.Barriers, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
lbsig_aldex2 <- aldex2_da_LB %>%
  filter(wi.eBH < 0.05)

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da_LB, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# make a table of significant corrected p-values
lbsig_aldex2 <- aldex2_da_LB %>%
  rownames_to_column(var = "Gene_") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(Gene_, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
lbsig_aldex2 <- left_join(lbsig_aldex2, taxa_info)
lbsig_aldex2
write.csv(sig_aldex2, file = "tables/aldex-cow-calf.csv")

ggplot(data = lbsig_aldex2, aes(x = Broadclass, y = effect)) +
  geom_point()

ggplot(data=lbsig_aldex2, aes(x=Broadclass, y=(effect), col=Broadclass, label=Gene)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = color_palette) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black")

ggsave(filename = "plots/aldex2-cow-calf.pdf", dpi = 600)



### TESTING FOR FUTURE ANALYSIS IF I CAN GET IT TO WORK ----
test <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps)), phyloseq::sample_data(ps)$Formal.Team.Meetings.Frequency, taxa_rank = "all", norm = "CLR", method = "kruskal", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")



# the output is in the S3 object 'x'
x <- aldex.clr(data.frame(phyloseq::otu_table(ps)), phyloseq::sample_data(ps)$Male.Female, mc.samples=128, denom="all", verbose=F)

x.kw <- aldex.kw(x)

x.effect <- aldex.effect(x, CI=T, verbose=FALSE, paired.test=FALSE)

x.all <- data.frame(x.kw,x.effect)

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(x.all, type="MW", test="kruskal", called.cex = 1, cutoff = 0.05)

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
x2 <- x.kw %>%
  filter(kw.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "Gene_")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(x.kw, type="MW", test="kruskal", called.cex = 1, cutoff = 0.05) ## ERROR

# make a table of significant corrected p-values
x3 <- x.all %>%
  rownames_to_column(var = "Gene_") %>%
  filter(kw.eBH < 0.05) %>%
  arrange(effect, kw.eBH) %>%
  dplyr::select(Gene_, effect, kw.ep, kw.eBH)

# add in previously formed taxa information to complete the table
x4 <- left_join(x3, taxa_info)
x4
write.csv(gsig_aldex2, file = "tables/aldex-cow-calf.csv")

ggplot(data=x4, aes(x=Broadclass, y=(kw.ep), col=Broadclass, label=Gene)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = color_palette) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black")

ggsave(filename = "plots/aldex2-cow-calf.pdf", dpi = 600)