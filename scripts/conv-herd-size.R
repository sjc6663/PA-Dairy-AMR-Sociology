# Herd Size ----

ps2 <- psrel %>% 
  ps_mutate(Herd.Size = case_when(
    Herd.Size == "0-50" ~ "<150",
    Herd.Size == "50-100" ~ "<150",
    Herd.Size == "100-150" ~ "<150",
    Herd.Size == "150-200" ~ ">150",
    Herd.Size == "200+" ~ ">150"
  ))

sample_data(ps2)$Herd.Size

ps2 <- ps2 %>% subset_samples(
  Conventional.Organic == "Conventional"
)

meltdf <- psmelt(ps2)

#filter out rows with 0 relabun 
out2 <- meltdf[meltdf$Abundance != 0.0000000000, ]

out2$Herd.Size <- factor(out2$Herd.Size,
                         levels = c("<150", ">150"))

ggplot(out2, aes(x=Broadclass, fill = Herd.Size)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette) +
  facet_wrap(~Group, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = "Herd Size",
       title = "C") +
  theme(axis.text= element_text(size = 12,face = "bold"),
        axis.title= element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold", hjust = 0),
        legend.position = "right",
        strip.background = element_rect(fill = color_palette[2]),
        strip.text.x = element_text(size = 12, color = "black",face = "bold")
  ) +
  #guides(fill="none")+
  #scale_y_continuous(expand = c(0,0.1), limits = c(0,7500))+
  coord_flip()

ggsave(filename = "plots/full-run/herd-size.pdf", dpi = 600)

psbclass <- aggregate_taxa(ps2, level = "Broadclass")

psbclass %>% plot_composition(group_by = "Herd.Size") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("C")

psbclass %>% plot_composition(average_by = "Herd.Size") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("Conventional Only, Herd Size")

ggsave(filename = "plots/deep-analysis/con-herd-size.pdf", dpi = 600)
