psbclass %>% plot_composition(group_by = "Conventional.Organic",  sample.sort = "affiliation", x.label = "affiliation") +
  # change y axis to be percentages instead of numbers
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position = "top") +
  labs(x = " ", fill=' ') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("A")


psbclass %>% plot_composition(average_by = "affiliation", sample.sort = "Conventional.Organic", x.label = "affiliation", group_by = "Conventional.Organic") +
  scale_y_continuous(labels = scales::percent) +
  #theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 30)) +
  labs(x = "Farm") 
#ggtitle("A")

org <- subset_samples(psbclass,
                      Conventional.Organic == "Organic")
A <- org %>% plot_composition(average_by = "affiliation", sample.sort = "Conventional.Organic", x.label = "affiliation", group_by = "Conventional.Organic") +
  scale_y_continuous(labels = scales::percent) +
  #theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 30)) +
  labs(x = "Farm") 
#ggtitle("A")

conv <- subset_samples(psbclass,
                       Conventional.Organic == "Conventional")
B <- conv %>% plot_composition(average_by = "affiliation", sample.sort = "Conventional.Organic", x.label = "affiliation", group_by = "Conventional.Organic") +
  scale_y_continuous(labels = scales::percent) +
  #theme(legend.position = "none") +
  scale_fill_viridis(option = "mako", discrete = TRUE) + 
  theme(text = element_text(size = 30)) +
  labs(x = "Farm") 
#ggtitle("A")

library(patchwork)

A/B
