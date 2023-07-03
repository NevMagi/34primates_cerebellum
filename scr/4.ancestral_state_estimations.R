### ----------------------------------------------------- ###
# Ancestral state reconstructions
# Roberto Toro, Katja Heuer 2018
# Script adapted by N. Magielse 
# for inclusion of cerebellar and ansiform area data - 2022
### ----------------------------------------------------- ###     

#------------------------------------
# SET-UP
#------------------------------------
library(ape)
#library(picante)
#library(caper)
library(Rphylopars)
library(phytools)
library(tidyverse)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(plyr)
library(dispRity)

# Set working directory
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/Final_code/")

source("./input/catn.R")

# Set output directory
outputdir <- "4.ancestral_state_estimations"
dir.create(outputdir)

#------------------------------------
# Load Data
#------------------------------------

## Load 10k tree, force it to be ultrametric by extending branches
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
tree13 <- force.ultrametric(read.nexus("./input/AnsiformAreaTree/consensusTree_10kTrees_Primates_Version3_13.nex"), method = "extend")
tree15 <- force.ultrametric(read.nexus("./input/AnsiformAreaTree/consensusTree_10kTrees_Primates_Version3_15.nex"), method = "extend")

## Load the phenotypic data
# We use the log-transformed data to visualize the vastly different absolute traits
# We create many variants of the same data frame here, so it should be plug-and-play later.
pheno <- read.csv2("./0.cleanInput/Phenotypes_outliersRM.csv", sep=';')  %>% dplyr::rename(species = WilsonReederName)
pheno <- pheno[c("species","CerebellarVol", "SurfaceArea" , "CerebralVol","AbsGI","FoldingLength","FoldingNumber","Lambda","Delta", "CrusVol_1_3", "CerebroCerebellar", "CerebellarCrus")]
pheno_raw <- read.csv2("./0.cleanInput/Phenotypes_raw_outliersRM.csv", sep=';')  %>% dplyr::rename(species = WilsonReederName)
pheno_raw <- pheno_raw[c("species","CerebellarVol", "SurfaceArea" , "CerebralVol","AbsGI","FoldingLength","FoldingNumber","Lambda","Delta", "CrusVol_1_3", "CerebroCerebellar", "CerebellarCrus")]

pheno.noRatio <- pheno[c("species","CerebellarVol", "SurfaceArea" , "CerebralVol","AbsGI","FoldingLength","FoldingNumber","Lambda","Delta", "CrusVol_1_3")]
pheno.AA <- pheno %>% drop_na(CrusVol_1_3)
pheno.AA.MOI <- pheno.AA[c("species","CerebellarVol","CerebralVol", "CrusVol_1_3")]
pheno.old <- pheno[c("species", "SurfaceArea" , "CerebralVol","AbsGI","FoldingLength","FoldingNumber","Lambda","Delta")]
pheno.CC <- pheno[c("species","CerebellarVol","CerebralVol")]
pheno.MOI <- pheno[c("species","CerebellarVol","CerebralVol", "CrusVol_1_3")]
pheno.outliers <- read.csv2("./0.cleanInput/Phenotypes.csv", sep=';') %>% dplyr::rename(species = WilsonReederName)
pheno.outliers <- pheno.outliers[c("species","CerebellarVol", "SurfaceArea" , "CerebralVol","AbsGI","FoldingLength","FoldingNumber","Lambda","Delta", "CrusVol_1_3", "CerebroCerebellar", "CerebellarCrus")]
pheno.old.outliers <- pheno.outliers %>% dplyr::select(species, SurfaceArea , CerebralVol, AbsGI, FoldingLength, FoldingNumber, Lambda, Delta)
phenoMedian <- read.csv2("./0.cleanInput/PhenotypesMedian_outliersRM.csv", sep=';')
phenoMedian.AA <- phenoMedian %>% drop_na(MedianCrus)


#------------------------------------
# FIT BROWNIAN MOTION MODELS
#------------------------------------

# The full table of measurements from the current study and Heuer et al., 2019.
# The Brownian Motion model is used to infer ancestral characters from extant species traits.
# We also create ancestral character estimations based on different combinations of traits, to asses its effect on node-wise estimates.
catn("fit a Brownian Motion model")
p.BM <- phylopars(trait_data=pheno, tree=tree)
sink("4.ancestral_state_estimations/4.1.BM_full.txt")
p.BM
sink()

# Since the ratios seem to be less explained by phylogeny, we exclude them to observe the effect on the BM reconstructions.
catn("fit a Brownian Motion model")
p.BM.noRatio <- phylopars(trait_data=pheno.noRatio, tree=tree)
sink("4.ancestral_state_estimations/4.2.BM_full-noRatio.txt")
p.BM.noRatio
sink()

# For cerebellar and cerebral volume only
catn("fit a Brownian Motion model")
p.BM.CC <- phylopars(trait_data=pheno.CC, tree=tree)
sink("4.ancestral_state_estimations/4.3.BM_cerebrocerebellar.txt")
p.BM.CC
sink()

## For cerebellar, cerebral, and ansiform area volumes (Measurements of Interest (MOIs)).
catn("fit a Brownian Motion model")
p.BM.MOI <- phylopars(trait_data=pheno.MOI, tree=tree)
sink("4.ancestral_state_estimations/4.4.BM_MOI.txt")
p.BM.MOI
sink()

## To show the effect that the two outliers have on the Brownian Motion models, we also fit a model for the full species data.
catn("fit a Brownian Motion model")
p.BM.outliers <- phylopars(trait_data=pheno.outliers, tree=tree)
sink("4.ancestral_state_estimations/4.5.BM_outliers.txt")
p.BM.outliers
sink()

# Lastly, we fit Brownian Motion models for only those observations that are complete (13 species, AA volumes).

# Full data
catn("fit a Brownian Motion model")
p.BM.13 <- phylopars(trait_data=pheno.AA, tree=tree13)
sink("4.ancestral_state_estimations/4.6.BM_full_13species.txt")
p.BM.13
sink()

# MOIs
catn("fit a Brownian Motion model")
p.BM.13.MOI <- phylopars(trait_data=pheno.AA.MOI, tree=tree13)
sink("4.ancestral_state_estimations/4.7.BM_MOI_13species.txt")
p.BM.13.MOI
sink()

#------------------------------------------
# SAVE ANCESTRAL STATE ESTIMATIONS TO CSV
#------------------------------------------
# ntips Number of species we have phenotypes for
ntips <- length(tree$tip.label)
# Number of nodes inside the tree (ntips-1) + including the tips
nnodes <- ntips + tree$Nnode

# ntips Number of species we have phenotypes for
ntips13 <- length(tree13$tip.label)
# Number of nodes inside the tree (ntips-1) + including the tips
nnodes13 <- ntips13 + tree13$Nnode

## Put into csv the info for all the estimated phenotypes at each internal node (=all ancestors)
catn("save data in csv")
csv <- cbind(p.BM$anc_recon[(ntips+1):nnodes,])
# lower and upper confidence interval
lci <- csv - sqrt(p.BM$anc_var[(ntips+1):nnodes,])*1.96
uci <- csv + sqrt(p.BM$anc_var[(ntips+1):nnodes,])*1.96
colnames(lci) <- paste(colnames(lci), "-95%CI", sep = "")
colnames(uci) <- paste(colnames(uci), "+95%CI", sep = "")
csv <- cbind(csv,lci)
csv <- cbind(csv,uci)
write.csv(csv, file="./4.ancestral_state_estimations/4.8.ACE_fullmodel.csv")

## Convert the MOI phenotypes into raw values, for supplemental table1.
raw.ace <- read_delim("./4.ancestral_state_estimations/4.8.ACE_fullmodel.csv", delim = ",")
raw.ace$CerebellarVol <- 10^raw.ace$CerebellarVol
raw.ace$`CerebellarVol-95%CI` <- 10^raw.ace$`CerebellarVol-95%CI`
raw.ace$`CerebellarVol+95%CI` <- 10^raw.ace$`CerebellarVol+95%CI`

raw.ace$CerebralVol <- 10^raw.ace$CerebralVol
raw.ace$`CerebralVol-95%CI` <- 10^raw.ace$`CerebralVol-95%CI`
raw.ace$`CerebralVol+95%CI` <- 10^raw.ace$`CerebralVol+95%CI`

raw.ace$CrusVol_1_3 <- 10^raw.ace$CrusVol_1_3
raw.ace$`CrusVol_1_3-95%CI` <- 10^raw.ace$`CrusVol_1_3-95%CI`
raw.ace$`CrusVol_1_3+95%CI` <- 10^raw.ace$`CrusVol_1_3+95%CI`

write_csv2(raw.ace, "./4.ancestral_state_estimations/4.8a.ACE_fullmodel_rawMOIs.csv")

## For ancestral state estimations of the MOIs, we also use the 13-species model.
catn("save data in csv")
csv <- cbind(p.BM.13.MOI$anc_recon[(ntips13+1):nnodes13,])
# lower and upper confidence interval
lci <- csv - sqrt(p.BM.13.MOI$anc_var[(ntips13+1):nnodes13,])*1.96
uci <- csv + sqrt(p.BM.13.MOI$anc_var[(ntips13+1):nnodes13,])*1.96
colnames(lci) <- paste(colnames(lci), "-95%CI", sep = "")
colnames(uci) <- paste(colnames(uci), "+95%CI", sep = "")
csv <- cbind(csv,lci)
csv <- cbind(csv,uci)
write.csv(csv, file="./4.ancestral_state_estimations/4.9.ACE_MOIs_13species.csv")


## We illustrate the effect of including different combinations of measurements for BM-model construction.

# Cerebellar and Cerebral Volumes
catn("save data in csv")
csv <- cbind(p.BM.CC$anc_recon[(ntips+1):nnodes,])
# lower and upper confidence interval
lci <- csv - sqrt(p.BM.CC$anc_var[(ntips+1):nnodes,])*1.96
uci <- csv + sqrt(p.BM.CC$anc_var[(ntips+1):nnodes,])*1.96
colnames(lci) <- paste(colnames(lci), "-95%CI", sep = "")
colnames(uci) <- paste(colnames(uci), "+95%CI", sep = "")
csv <- cbind(csv,lci)
csv <- cbind(csv,uci)
write.csv(csv, file="./4.ancestral_state_estimations/4.10.ACE_cerebrocerebellar.csv")

# Cerebellar, Cerebral, and Ansiform Area Volumes
catn("save data in csv")
csv <- cbind(p.BM.MOI$anc_recon[(ntips+1):nnodes,])
# lower and upper confidence interval
lci <- csv - sqrt(p.BM.MOI$anc_var[(ntips+1):nnodes,])*1.96
uci <- csv + sqrt(p.BM.MOI$anc_var[(ntips+1):nnodes,])*1.96
colnames(lci) <- paste(colnames(lci), "-95%CI", sep = "")
colnames(uci) <- paste(colnames(uci), "+95%CI", sep = "")
csv <- cbind(csv,lci)
csv <- cbind(csv,uci)
write.csv(csv, file="./4.ancestral_state_estimations/4.11.ACE_MOI.csv")

# Full model with outliers
catn("save data in csv")
csv <- cbind(p.BM.outliers$anc_recon[(ntips+1):nnodes,])
# lower and upper confidence interval
lci <- csv - sqrt(p.BM.outliers$anc_var[(ntips+1):nnodes,])*1.96
uci <- csv + sqrt(p.BM.outliers$anc_var[(ntips+1):nnodes,])*1.96
colnames(lci) <- paste(colnames(lci), "-95%CI", sep = "")
colnames(uci) <- paste(colnames(uci), "+95%CI", sep = "")
csv <- cbind(csv,lci)
csv <- cbind(csv,uci)
write.csv(csv, file="./4.ancestral_state_estimations/4.11.ACE_outliers.csv")

#------------------------------------------
# PLOT ANCESTRAL STATE ESTIMATIONS
#------------------------------------------

## Define a two-color scaled function for plotting ancestral state estimations along the phylogenetic tree. ##
myplotanc <- function(tree, tips, tipnames, states, mytitle) {
  names(tips) <- tipnames
  obj <- contMap(tree,tips,method="user",anc.states=states, plot="F")
  obj <- setMap(obj, colors=c("white", "darkblue"))
  plot.contMap(obj)
  axisPhylo()
  title(mytitle)
}

myplotanc.ci.Taxo <- function(tree, tips, tipnames, states, legend.title, rangeCIs) {
  names(tips) <- tipnames
  obj <- contMap(tree,tips,method="user",anc.states=states, plot="F", lims = rangeCIs)
  obj <- setMap(obj, colors=c("blue", "white", "red"))
  plot.contMap(obj, legend=FALSE,ylim=c(1-0.09*(Ntip(obj$tree)-1),Ntip(obj$tree)),
               mar=c(5.1,0.4,0.4,0.4))
  add.color.bar(36.502,obj$cols,title= legend.title, y=3,
                lims=obj$lims,digits=2,prompt=FALSE,x=0,
                y=1-0.08*(Ntip(obj$tree)-1),lwd=4,fsize=1,subtitle="")
  axisPhylo()
  title(xlab="Time to present (in million years)", adj = 0.17)
  errorbar.contMap(obj, lwd=2)
}

myplotanc.ci.noTaxo <- function(tree, tips, tipnames, states, legend.title, rangeCIs) {
  names(tips) <- tipnames
  obj <- contMap(tree,tips,method="user",anc.states=states, plot="F", ftype = "off", lims = rangeCIs)
  obj <- setMap(obj, colors=c("blue", "white", "red"))
  plot.contMap(obj, legend=FALSE,ylim=c(1-0.09*(Ntip(obj$tree)-1),Ntip(obj$tree)),
               mar=c(5.1,0.4,0.4,0.4), ftype = "off")
  add.color.bar(36.502,obj$cols,title= legend.title, y=3, fsize=1.0,
                lims=obj$lims,digits=2,prompt=FALSE,x=0,
                y=1-0.08*(Ntip(obj$tree)-1), fsize=0.8,subtitle="")
  axisPhylo()
  title(xlab="Time to present (in million years)", adj = 0.17)
  errorbar.contMap(obj, lwd=2)
}

catn("plot ancestral values")
tipnames <- rownames(p.BM$anc_recon)[1:ntips]
tipnames13 <- rownames(p.BM.13$anc_recon)[1:ntips13]

# Evolution cerebellar volume across the tree
pdf("./4.ancestral_state_estimations/4.12.BM-anc-cerebellarVolume.pdf")
myplotanc(tree, p.BM$anc_recon[1:ntips,"CerebellarVol"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CerebellarVol"], "Ancestral cerebellar volume reconstruction")
graphics.off()
tiff("./4.ancestral_state_estimations/4.12a.BM-anc-cerebellarVolume-ci.tiff", width = 4, height = 8, units = 'in', res = 300)
myplotanc.ci.noTaxo(tree, p.BM$anc_recon[1:ntips,"CerebellarVol"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CerebellarVol"], legend.title = "Cerebellar volume \n (log10 mm3)", rangeCIs = c(2.27,5.00))
graphics.off()

#---------------------------------------------------------------------------------
# INTERLUDE: CORRELATION CI (UNCERTAINTY) WITH INCREASING TIME TO PRESENT
# Showing that ancestral character estimations become uncertain further back in time
# This is not a groundbreaking observation, but good to know nonetheless.
#---------------------------------------------------------------------------------

## CIs seem to be correlated with evolutionary time. Let's check.
node.age <- tree.age(tree, order = "past", fossil = TRUE, digits = 3)
node.ci <- p.BM$anc_var %>% as.data.frame() %>% rownames_to_column(var = "elements")
node.ci.vs.age.df <- left_join(node.age, node.ci, by="elements") %>% slice_tail(n = 33) %>%
  melt(id.var="ages") %>% slice_tail(n = 363) 
node.ci.vs.age.df

# Add regression coefficients and slope.
regression=function(df){
  #setting the regression function. 
  reg_fun<-lm(formula=df$value~df$ages) #regression function
  #getting the slope, intercept, R square and adjusted R squared of 
  #the regression function (with 3 decimals).
  slope<-round(coef(reg_fun)[2],3)  
  intercept<-round(coef(reg_fun)[1],3) 
  R2<-round(as.numeric(summary(reg_fun)[8]),3)
  R2.Adj<-round(as.numeric(summary(reg_fun)[9]),3)
  c(slope,intercept,R2,R2.Adj)
}
coeff.ages <- ddply(node.ci.vs.age.df, "variable", regression)
colnames(coeff.ages) <-c ("variable","slope","intercept","R2","R2.Adj")

# New facet label names for the neuroanatomical phenotypes.
var.labs <- c("Cerebellar volume", "Cerebral surface area", "Cerebral volume", "AbsGI", "Folding length", "Folding number", "Lambda", "Delta", "Ansiform area volume", "Cerebellar/ cerebral volume", "Ansiform area/ cerebellar volume")
names(var.labs) <- c("CerebellarVol", "SurfaceArea", "CerebralVol", "AbsGI", "FoldingLength", "FoldingNumber", "Lambda", "Delta", "CrusVol_1_3", "CerebroCerebellar", "CerebellarCrus")

tiff("./4.ancestral_state_estimations/4.13.BM_ANC_CI_correlation_time.tiff", width = 12, height = 8, units = 'in', res = 300)
node.ci.vs.age <- left_join(node.age, node.ci, by="elements") %>% slice_tail(n = 33) %>% 
  melt(id.var="ages") %>% slice_tail(n = 363) %>% ggplot(aes(ages,value)) + 
  geom_point() + 
  stat_smooth() +
  geom_smooth(method = "lm") +
  facet_wrap(~variable, scales = "free", labeller = labeller(variable = var.labs)) + 
  theme_classic() +
  theme(axis.text.y = element_blank()) + 
  geom_label(data=coeff.ages, inherit.aes=FALSE, aes(x = 42, y = 4,
                                                     label=paste("slope=",slope,","," ","R^2.Adj=",R2.Adj)),label.size = NA) +
  labs(x = "Time to present (in million years)", y = "Covariance estimates") +
  theme(axis.title.x = element_text(size = 18, face = "plain"),
        axis.title.y = element_text(size = 18, face = "plain")) +
  theme(axis.title.x = element_text(margin = margin(t = 30))) +
  theme(axis.title.y = element_text(margin = margin(r = 30))) +
  #ggtitle("Variance of Neuroanatomical Measurements at Ancestral Nodes") +
  theme(plot.title = element_text(size = 24, face = "plain", hjust = 0.5, margin = margin(b = 30))) +
  theme(
    strip.text.x = element_text(
      size = 12, face = "plain"
    ),
    strip.text.y = element_text(
      size = 12, face = "plain"
    )
  )
node.ci.vs.age
graphics.off()


#------------------------------------------
# CONTINUATION: PLOT ANCESTRAL STATE ESTIMATIONS
#------------------------------------------

# Ancestral Ansiform Area Volume
pdf("./4.ancestral_state_estimations/4.14.BM-anc-ansiformareaVolume.pdf")
myplotanc(tree13, p.BM.13$anc_recon[1:ntips13,"CrusVol_1_3"], tipnames13, p.BM.13$anc_recon[(ntips13+1):nnodes13,"CrusVol_1_3"], "Ansiform Area Volume Reconstruction")
graphics.off()
tiff("./4.ancestral_state_estimations/4.14a.BM-anc-ansiformareaVolume-ci.tiff", width = 4, height = 8, units = 'in', res = 300)
myplotanc.ci.noTaxo(tree13, p.BM.13$anc_recon[1:ntips13,"CrusVol_1_3"], tipnames13, p.BM.13$anc_recon[(ntips13+1):nnodes13,"CrusVol_1_3"], "Ansiform area volume \n (log10 mm3)", rangeCIs = c(1.76, 4.94))
graphics.off()

contMap(tree,y,lims=lims,method="anc.ML",legend=1.5)
# Ancestral cerebral volume
pdf("./4.ancestral_state_estimations/4.15.BM-anc-cerebralVolume.pdf")
myplotanc(tree, p.BM$anc_recon[1:ntips,"CerebralVol"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CerebralVol"], "Cerebral Volume Reconstruction")
graphics.off()
tiff("./4.ancestral_state_estimations/4.15a.BM-anc-cerebralVolume-ci.tiff", width = 4, height = 8, units = 'in', res = 300)
myplotanc.ci.noTaxo(tree, p.BM$anc_recon[1:ntips,"CerebralVol"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CerebralVol"], "Cerebral volume \n (log10 mm3)", rangeCIs = c(3.03,5.78))
graphics.off()

# Ancestral Cerebellar/Cerebral Ratio (in percentages)
pdf("./4.ancestral_state_estimations/4.16.BM-anc-cerebelloCerebralRatio.pdf")
myplotanc(tree, p.BM$anc_recon[1:ntips,"CerebroCerebellar"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CerebroCerebellar"], "Reconstruction Ratio Cerebellar-to-Cerebral Volume")
graphics.off()
tiff("./4.ancestral_state_estimations/4.16a.BM-anc-cerebelloCerebralRatio-ci.tiff", width = 8, height = 8, units = 'in', res = 300)
myplotanc.ci.Taxo(tree, p.BM$anc_recon[1:ntips,"CerebroCerebellar"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CerebroCerebellar"], "Cerebellar/cerebral \n volume (%)", rangeCIs = c(8.06, 23.03))
graphics.off()

# Ancestral Ansiform Area/Cerebellar Ratio (in percentages)
pdf("./4.ancestral_state_estimations/4.17.BM-anc-aaCerebellarRatio.pdf")
myplotanc(tree13, p.BM.13$anc_recon[1:ntips13,"CerebellarCrus"], tipnames13, p.BM.13$anc_recon[(ntips13+1):nnodes13,"CerebellarCrus"], "Reconstruction Ratio Ansiform Area-to-Cerebellar Volume")
graphics.off()
tiff("./4.ancestral_state_estimations/4.17a.BM-anc-aaCerebellarRatio-ci.tiff", width = 8, height = 8, units = 'in', res = 300)
myplotanc.ci.Taxo(tree13, p.BM.13$anc_recon[1:ntips13,"CerebellarCrus"], tipnames13, p.BM.13$anc_recon[(ntips13+1):nnodes13,"CerebellarCrus"], "Ansiform area/ \n cerebellar volume (%)", rangeCIs = c(-2.1, 33.4))
graphics.off()


#------------------------------------------
# PLOT PHYLOGRAMS
#------------------------------------------

# Define a phenogram function (not used currently)
catn("Plot phenograms (BM model)")
plotpheno<-function(ph,tr,trait,mytitle, my.xlab, my.ylab) {
  t <- tree.age(tree, order = "past", fossil = TRUE, digits = 3)
  h<-p.BM$anc_recon[,trait]
  plot(c(min(t$ages),max(t$ages)*1.1),c(min(h),max(h)),type="n")
  length(tr$edge)
  for(i in 1:nrow(tr$edge)) {
    ed<-tr$edge[i,]
    segments(t$ages[ed[1]],h[ed[1]],t$ages[ed[2]],h[ed[2]])
  }
  for(i in 1:length(tree$tip.label)) {
    text(t$ages[i],h[i],tree$tip.label[i], adj=c(0,0), cex=0.1)
  }
  title(mytitle)
}

## Reverse time (scale). Also ultimately not used.
#xa<-c(phenoMedian$MedianCerebellum,fastAnc(tree,phenoMedian$MedianCerebellum))
#H<-nodeHeights(tree)
#H<-max(H)-H
#Y<-matrix(xa[tree$edge],nrow(tree$edge),2)
#plot.new()
#plot.window(xlim=range(H)[2:1],ylim=range(Y))
#axis(1)
#axis(2)
#for(i in 1:nrow(tree$edge)) lines(H[i,],Y[i,],lwd=2)
#title(xlab="Time in Mya",ylab="Cerebellar Volume")

## Extract characters of interest
cerebellumPhylogram <- setNames(phenoMedian$MedianCerebellum,
                                phenoMedian$WilsonReederName)
cerebrumPhylogram <- setNames(phenoMedian$MedianCerebrum,
                              phenoMedian$WilsonReederName)
cerebrocerebellarPhylogram <- setNames(phenoMedian$MedianCerebroCerebellar,
                                       phenoMedian$WilsonReederName)
ansiformPhylogram <- setNames(phenoMedian.AA$MedianCrus,
                              phenoMedian.AA$WilsonReederName)
cerebellaransiformPhylogram <- setNames(phenoMedian.AA$MedianCerebellarCrus,
                                        phenoMedian.AA$WilsonReederName)

# Give the hominoidea a unique color
plotTree(tree,pts=F,node.numbers=T) #check from what node on the hominoidea starts 
# node 47 (internal node)
tree<-paintSubTree(tree,node=47,state="2")
plot(tree)
cols<-c("black","gold"); names(cols)<-1:2
plotSimmap(tree,cols,pts=F,lwd=3,node.numbers=T)

## Plot a traitgram for the cerebellum
tiff("./4.ancestral_state_estimations/4.18.BM-anc-cerebellarVolume-phenogram.tiff", width = 10, height = 10, units = 'in', res = 300)
phenogram(tree,cerebellumPhylogram, colors = cols, ftype="i",
          spread.cost=c(1,0),fsize=0.7,xlab="Time since the root (in Mya)",
          ylab="log(Cerebellar Volume)")
graphics.off()

## Plot a traitgram for the cerebrum
tiff("./4.ancestral_state_estimations/4.19.BM-anc-cerebralVolume-phenogram.tiff", width = 10, height = 10, units = 'in', res = 300)
phenogram(tree,cerebrumPhylogram, colors = cols, ftype="i",
          spread.cost=c(1,0),fsize=0.7,xlab="Time since the root (in Mya)",
          ylab="log(Cerebral Volume)")
graphics.off()

## Plot a traitgram for the ratio of Cerebellar-to-Cerebral Volume
tiff("./4.ancestral_state_estimations/4.20.BM-anc-cerebrocerebellarRatio-phenogram.tiff", width = 10, height = 10, units = 'in', res = 300)
phenogram(tree,cerebrocerebellarPhylogram, colors = cols, ftype="i",
          spread.cost=c(1,0),fsize=0.7,xlab="Time since the root (in Mya)",
          ylab="Ratio Cerebellar/Cerebral Volume")
graphics.off()

# For the 13 species tree: give the hominoidea a unique color, as well.
plotTree(tree13,pts=F,node.numbers=T) #check from what node on the hominoidea starts 
# node 19 (internal node)
tree13<-paintSubTree(tree13,node=19,state="2")
plot(tree13)
cols<-c("black","gold"); names(cols)<-1:2
plotSimmap(tree13,cols,pts=F,lwd=3,node.numbers=T)

## Plot a traitgram for the Ansiform Area
tiff("./4.ancestral_state_estimations/4.21.BM-anc-ansiformVolume-phenogram.tiff", width = 10, height = 10, units = 'in', res = 300)
phenogram(tree13,ansiformPhylogram, colors = cols, ftype="i",
          spread.cost=c(1,0),fsize=0.7,xlab="Time since the root (in Mya)",
          ylab="log(Ansiform Area Volume)")
graphics.off()

## `Plot a traitgram for the ratio of Ansiform area to Cerebellar Volume
tiff("./4.ancestral_state_estimations/4.22.BM-anc-cerebellaransiformRatio-phenogram.tiff", width = 10, height = 10, units = 'in', res = 300)
phenogram(tree13,cerebellaransiformPhylogram, colors = cols, ftype="i",
          spread.cost=c(1,0),fsize=0.7,xlab="Time since the root (in Mya)",
          ylab="Ratio Ansiform Area/ Cerebellar Volume")
graphics.off()

### -------------------------- ###
### END of script
### -------------------------- ###   
