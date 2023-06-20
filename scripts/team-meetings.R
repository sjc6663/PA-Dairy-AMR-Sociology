sample_data(psrel)$Formal.Team.Meetings.Frequency

is.character(sample_data(psrel)$Formal.Team.Meetings.Frequency)

library(dplyr)
library(ggplot2)

psrel <- psrel %>% 
  ps_mutate(
    meetings = if_else(str_detect(Formal.Team.Meetings.Frequency, "Never"), true = "No", false = "Yes")
  ) 

ggplot(out3, aes(x=Broadclass, fill = meetings)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette) +
  facet_wrap(~Male.Female, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = "Meetings",
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

