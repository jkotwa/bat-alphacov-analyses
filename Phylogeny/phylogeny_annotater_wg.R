#load libraries
library(ggtree)
library(ggimage)
library(ggplot2)
library(rphylopic)
library(treeio)
library(phylotools)

# Set working directory
setwd("")

# Import tree
tree <- read.tree(".tre")

# Import metadata - new labels, phylopics
meta <- read.table(".csv", sep=",", header=T, comment.char = "")

# Extract tip data
tip_data <- data.frame(
  tip.label=c(meta$Label),
  image=c(meta$Phylo),
  colour=c(meta$Colour))

# Sub taxa label
treelab <- sub.taxa.label(tree, meta)

# Merge tree data with tip data
tree_data <- full_join(as_tibble(treelab), tip_data, by = c("label" = "tip.label"))

# Generate tree to get node data
p <- ggtree(treelab, ladderize = TRUE) %<+% tree_data

# Generate node support
node_data <- data.frame(
  node=c(p$data$node),
  node.label=c(as.numeric(p$data$label)),
               node.labelsub=c(100-as.numeric(p$data$label)) 
)

# Hide bootstrap values under 75%
node_data$node.label[node_data$node.label <= 74] = NA #hide bootstraps with <75%
node_data$node.labelsub[node_data$node.labelsub >= 26] = NA #hide bootstraps with <75%

# Generate colours
colours <- data.frame(colours=c(p$data$colour))
colours <- na.omit(colours)

# Pie charts
pies <- nodepie(node_data, cols=2:3)
pies <- lapply(pies, function(p) p+scale_fill_manual(values = c("#3F0954","white")))

# Plot the tree and add images
p <- ggtree(treelab, ladderize = TRUE) %<+% tree_data + 
  geom_tiplab(aes(label = label), geom = "label", fill = colours, alpha = 0.33, label.size=0, hjust=-0.01, size=2.5, offset=0, align=TRUE) +
  geom_tiplab(aes(image = image), geom = "phylopic", size = 0.03, offset=0, linetype = NULL, align=TRUE) +
  geom_inset(pies, x="node", width=0.023, height=0.023, hjust=0.02) +
  geom_treescale(x=0, fontsize = 2.5) +
  xlim(NA, 2.95) +
  theme_tree()

# Save the plot
ggsave(".png", p, dpi=600, width = 10, height = 8)