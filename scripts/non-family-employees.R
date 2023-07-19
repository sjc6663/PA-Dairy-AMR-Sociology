# Non-Family Employees ----

sample_data(psrel)

sample_data(psrel)$Non.Family.Milkers

is.character(sample_data(psrel)$Non.Family.Milkers)
as.character(sample_data(psrel)$Non.Family.Milkers)

psrel <- psrel %>% 
  ps_mutate(
    employees = if_else(str_detect(Non.Family.Milkers, "0"), true = "No", false = "Yes")
  ) 

sample_data(ps)

meltdf <- psmelt(psrel)

#filter out rows with 0 relabun 
out2 <- meltdf[meltdf$Abundance != 0.0000000000, ]


ggplot(out2, aes(x=Broadclass, fill = employees)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette) +
  facet_wrap(~Group, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = "Employees",
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

psbclass <- aggregate_taxa(psrel, level = "Broadclass")

psbclass %>% plot_composition(group_by = "employees") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  ggtitle("C")

psbclass %>% plot_composition(average_by = "employees") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("Employees")

ggsave(filename = "plots/deep-analysis/employees.pdf", dpi = 600)
