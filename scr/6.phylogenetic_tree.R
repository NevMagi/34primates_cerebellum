
### -------------------------------------- ###
# Create an annotated phylogenetic tree
# N. Magielse, 2022
### -------------------------------------- ###

#------------------------------------
# Set-up
#------------------------------------
setwd("~/Documents/PhD/EvolutionPrimates/Final_code/")
outputdir <- "6.phylogenetic_tree"
dir.create(outputdir)

library(ape)
library(phytools)

#------------------------------------
# Load Data
#------------------------------------
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
tree13 <- force.ultrametric(read.nexus("./input/AnsiformAreaTree/consensusTree_10kTrees_Primates_Version3_13.nex"), method = "extend")

## Plot the phylogenetic tree, annotated with the internal node number and archaelogical time.
tiff("./6.phylogenetic_tree/6.1.annotated_Tree.tiff", width = 8, height = 10, units = 'in', res = 300)
h<-max(nodeHeights(tree))
plotTree(tree,plot)
obj <- geo.legend(alpha = 0.3, cex = 1.2)
#obj$leg<-h-obj$leg
plotTree(tree,ylim=c(-0.2*Ntip(tree),Ntip(tree)),lwd=1, ftype = "off")
geo.legend(leg=obj$leg,colors=obj$colors,cex=1.2)
# axisPhylo(side = 3)
nodelabels(cex=0.60, bg = "#ffffff", frame = "circle")
graphics.off()

### -------------------------- ###
### END of script
### -------------------------- ###   

