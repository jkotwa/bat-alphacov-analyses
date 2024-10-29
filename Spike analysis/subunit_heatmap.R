library("ggplot2")
library("tidyr")

setwd("")

data <- read.table(".csv", sep=",", header=T, stringsAsFactor=F)

data$Label <- as.character(data$Label)
data$Label <- factor(data$Label, levels=unique(data$Label))

s_heat <- ggplot(data=data, mapping = aes(x=identity,
                                          y=Label,
                                          fill=identity.1), fsize=1) +
  geom_tile(colour="white") +
  facet_grid(~subunit, ) +
  scale_x_discrete(position="top", guide=guide_axis(angle=45)) +
  xlab("") +
  ylab("") +
  scale_fill_gradient(name = "% Identity",
                      low = "#FFFFFF",
                      high = "#150c25") + 
  theme(strip.placement = "outside",
        axis.text.x = element_text(angle=90, size=5),
        axis.text.y = element_text(size=6.5)) +
  coord_fixed() 

# Save the plot
ggsave(".png", s_heat, dpi=600, width = 5, height = 6)

