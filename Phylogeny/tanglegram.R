library(ggplot2)
library(ggtree)
library(phangorn)
library(dplyr)
library(phylotools)
library(phytools)
library(TreeTools)

setwd("")

#load metadata for each phylogeny - make sure tips match
meta <- read.table("alpha_meta_short_clades.csv", 
                   sep=",", header=T, comment.char = "")
meta2 <- read.table("alpha_meta_tangle_spike.csv", 
                    sep=",", header=T, comment.char = "")

#load tree 1
tree1 <- read.tree(".tre")


#load tree 2
tree2 <- read.tree(".tre")


# Create a named vector of colors
colour_map <- setNames(meta$Colour, meta$Label)

#replace tip labels with cleaned labels
tree1lab <- sub.taxa.label(tree1, meta)
tree2lab <- sub.taxa.label(tree2, meta2)


tree1lab$edge.length<-NULL
tree2lab$edge.length<-NULL

#generate cophylogeny
tangle <- cophylo(tree1lab, tree2lab, assoc = NULL, rotate=TRUE)

# Extract the association matrix from the cophylogeny object
assoct <- tangle$assoc


# Function to map colors based on taxon labels
get_link_colour <- function(assoct, colour_map) {
  sapply(assoct, function(label) {
    if (label %in% names(colour_map)) {
      colour_map[label]
    } else {
      "black"  # Default color if not specified
    }
  })
}

# Get link colors based on the tip labels
link_colours <- get_link_colour(assoct, colour_map)

#generate figure and export as PDF
pdf(".pdf",width=20, height=10, bg="white")
plot(tangle, show.tip.label=TRUE, link.type="curved", link.lwd=6, 
     link.lty="solid", link.col=make.transparent(link_colours, 0.95), 
     cex=1, fsize=1)

nodelabels.cophylo(node=1:tangle$trees[[1]]$Nnode+Ntip(tangle$trees[[1]]), 
                   pie=cbind(as.numeric(tangle$trees[[1]]$node.label),100- as.numeric(tangle$trees[[1]]$node.label)), 
                   piecol=c("black","white"),cex=0.15,which="left")
nodelabels.cophylo(node=1:tangle$trees[[2]]$Nnode+Ntip(tangle$trees[[2]]), 
                   pie=cbind(as.numeric(tangle$trees[[2]]$node.label),100- as.numeric(tangle$trees[[2]]$node.label)), 
                   piecol=c("black","white"),cex=0.15,which="right")
dev.off()