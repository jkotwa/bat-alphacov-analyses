# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(ggforce)
library(reshape2)
library(gridExtra)
library(ggchicklet)
library(dplyr)

# Set working directory
setwd("")

# Load data
coverage_data <- read.csv(".csv")
gene_coords <- read.csv(".csv")

# Reshape coverage data for plotting
coverage_data_long <- reshape2::melt(coverage_data, id.vars = "position", variable.name = "sample", value.name = "similarity")

# Create the coverage plot
coverage_plot <- ggplot(coverage_data_long, aes(x = position, y = similarity, color = sample)) +
  geom_line() +
  scale_color_manual(values=c("#3B84a1","#172773","#377373", "#732737")) +
  scale_y_continuous(breaks = seq(0, 1, by=0.05), limits = c(0, 1), expand=c(0,0)) +
  scale_x_continuous(breaks = seq(0, max(gene_coords$end), by = 5000), limits = c(0, max(gene_coords$end)), expand=c(0,0)) +
  labs(x = "", y = "Similarity", color = "") +
  theme_minimal() + 
  theme(
    legend.position = "top",
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(size = 1, color = "black"),
    axis.ticks = element_line(size = 1, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 0, 10)
  ) 

# Create the genome organization plot
genome_plot <- ggplot(gene_coords) +
  ggchicklet:::geom_rrect(aes(xmin = start, xmax = end, ymin = 0.98, ymax = 1.02), r=unit(0.3,'npc'), fill = "#001e4e", alpha = 0.85) +
  geom_text(aes(x = (start + end) / 2, y = 1, label = gene), col="white", size = 2, fontface=2, hjust = 0.5, vjust = 0.5) + 
  scale_x_continuous(breaks = seq(0, max(gene_coords$end), by = 5000),  expand=c(0,0)) +
  labs(x = "Genome position") +  # Add x-axis label here
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = "white", ),
        plot.margin = margin(0, 5, 20, 10),
        axis.title.x = element_text(size = 10, margin = margin(t=5))) 

# Combine the two plots
combined_plot <- ggarrange(coverage_plot, genome_plot, ncol = 1, nrow=2,  heights = c(2, 0.3), align = "v")

# Save the plot
ggsave(".png", combined_plot, width = 8, height = 5)