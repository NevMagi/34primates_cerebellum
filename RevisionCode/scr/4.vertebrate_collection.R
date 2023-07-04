############################################
### Analysis carried out in response
### to Reviewer #1. Specifically, these
### analyses seek to explore whether the main
### result holds when we only consider species
### in the Paris Vertebrate Collection.
### This analysis seeks to provide additional
### perspective on the contribution of difference
### in preparation and preparation between
### brains.This is assumed to be more similar
### in the Paris Collection than versus the
### rest of the data.
###
### NB: Macaca mulatta and Macaca fascicularis
### have been removed from the data, because
### they have a variable provenance.
############################################

library(tidyverse)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ape)
library(phytools)
library(Rphylopars)
library(dispRity)
library(nlme)

### SET UP ###
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/RevisionCode")

outputdir <- "4.vertebrate_collection"
dir.create(outputdir)

## Data input
paris <- read_csv2("./input/PhenotypesMedian_split_PARIS.csv") %>% dplyr::rename(species = WilsonReederName) %>% dplyr::rename(Mass = `Body mass species mean` ) %>% select(MedianCerebellum, MedianCerebrum, Mass, Paris) %>% as.data.frame()
paris$Paris <- as.factor(paris$Paris)
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
row.names(paris) <- tree[["tip.label"]] # Depends on alphabet sorting of your data.
paris <- na.omit(paris) %>% filter(!Paris=='Partial')
tree<-treedata(tree,paris, sort=T,warnings=T)$phy # Match tree to the data


#One could manually assign groupings and filter based on that.
PARIS <-paris %>% filter(Paris== "Yes")
PARIS <- PARIS[-c(18),] # Manually remove the rhesus macaque row
tree_paris <- treedata(tree,PARIS, sort=T,warnings=T)$phy # Match tree to the paris data
## We ommitted data from the rhesus macaque here manually, because its value was an outlier in intraspecific analysis.
## We suspect that for this specimen, ex vivo preparation compounded the already small brain (likely due to it being a juvenile specimen).


#One could manually assign groupings and filter based on that.
NOPARIS <-paris %>% filter(Paris== "No")
tree_noparis <- treedata(tree,NOPARIS, sort=T,warnings=T)$phy # Match tree to the paris data 



#####################################
### CALCULATE PICs
#####################################
Cparis <- PARIS$MedianCerebellum
Vparis <- PARIS$MedianCerebrum

Cnoparis <- NOPARIS$MedianCerebellum
Vnoparis <- NOPARIS$MedianCerebrum

# Compute contrasts
pic.C.paris <- pic(Cparis,tree_paris)
pic.V.paris <- pic(Vparis, tree_paris)

pic.C.noparis <- pic(Cnoparis,tree_noparis)
pic.V.noparis <- pic(Vnoparis, tree_noparis)


#------------------------------------
# PIC-regressions Paris
#------------------------------------
# Compute a linear regression of C on V
result <- lm(pic.C.paris ~ pic.V.paris)
sink("4.vertebrate_collection/4.1.lm-C-V.paris.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("4.vertebrate_collection/4.2.lm-C-V-mean-nointercept-plot.paris.pdf")
plot(pic.C.paris ~ pic.V.paris)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.
## They are the same; check! ##

## Check if PIC/PGLS models are the same.
bm.paris <-corBrownian(1,tree_paris)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=PARIS,correlation=bm.paris)
summary(pglsModel)
PICmodel <- lm(pic.C.paris ~pic.V.paris + 0)
sink("4.vertebrate_collection/4.3.PGLS-PIC-C-V.paris.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


#------------------------------------
# PIC-regressions non-Paris
#------------------------------------
# Compute a linear regression of C on V
result <- lm(pic.C.noparis ~ pic.V.noparis)
sink("4.vertebrate_collection/4.4.lm-C-V.noparis.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("4.vertebrate_collection/4.5.lm-C-V-mean-nointercept-plot.noparis.pdf")
plot(pic.C.noparis ~ pic.V.noparis)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.
## They are the same; check! ##

## Check if PIC/PGLS models are the same.
bm.noparis <-corBrownian(1,tree_noparis)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=NOPARIS,correlation=bm.noparis)
summary(pglsModel)
PICmodel <- lm(pic.C.noparis ~pic.V.noparis + 0)
sink("4.vertebrate_collection/4.6.PGLS-PIC-C-V.noparis.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


###########################################################
### Use pANCOVA to test deviations in intercept or slope
### and with that deviation from the general allometric 
### relationship for cerebellar-to-cerebral scaling depending
### on provenance (Paris Collection).
###########################################################


library(evomap)
library(geiger)

# This is an alternative way to subset the trees :). Very convenient
Paris <-getTips(tree,findMRCA(tree,tips = tree_paris[["tip.label"]]))
Noparis <-getTips(tree,findMRCA(tree,tips = tree_noparis[["tip.label"]]))


#------------------------------------
# pANCOVA - Paris
#------------------------------------
#For differences in slope: 
grpS<-rep("A",length(rownames(paris))) 
grpS[Noparis]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(paris)

#For differences in intercept: 
grpI<-rep("A",length(rownames(paris))) 
grpI[Noparis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(paris)


## Data subsetting
## Cerebellum vs. Cerebrum
data<-primateData; tree<-primateTree
Y<-"Brain"; X<-"Body"
data<-data[,c(which(colnames(data)==Y),which(colnames(data)==X)),drop=F]
data<-na.omit(data); data<-log(data)
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)   #match the data to the tree
colnames(data)<-c("Dependent","Independent")  
## Very hacky but bypasses weird error.
data <- data[1:33,]
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
tree <-treedata(tree,paris, sort=T,warnings=T)$phy
rownames(data) <- tree[["tip.label"]]

## Plot separate PGLS regressions for Paris and non-Paris data
## Cerebellum vs. Cerebrum

## Again, extremely hacky, but works.
data$Dependent <- paris$MedianCerebellum
data$Independent <- paris$MedianCerebrum

pdf("4.vertebrate_collection/4.7.PGLS.Paris.CC.pdf")
plot(data$Dependent~data$Independent,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=Paris,col="green",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=Noparis,col="grey",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.CC <-model.matrix(as.formula(Dependent~Independent),data)
Model.CC.s <-model.matrix(as.formula(Dependent~grpS:Independent),data) #slopes differ
Model.CC.i <-model.matrix(as.formula(Dependent~grpI + Independent),data) #intercepts differ
Model.CC.si <-model.matrix(as.formula(Dependent~grpI + grpS:Independent),data) #slopes & intercepts differ

sink("4.vertebrate_collection/4.8.ModelTesting.SlopesIntercepts.CC.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.si)
sink()


###-----------------------------------------------###
# Plot figures supplement - MNHN data
####----------------------------------------------###

### -------------------------- ###
### PGLS Plotting with Caper
### -------------------------- ###   
library(caper)
library(evomap)

# Model 1: Cerebellar volume regressed on cerebral volume in MNHN
tiff("./4.vertebrate_collection/4.9.BM-PGLS-cerebellumCerebrum-MNHN-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod1 <- gls(MedianCerebellum ~ MedianCerebrum, data=PARIS, correlation=bm.paris)
plot(mod1)
plot(PARIS$MedianCerebellum ~ PARIS$MedianCerebrum, xlab = "Cerebral volume (log10-transformed)", ylab = "Cerebellar volume (log10-transformed)")
abline(mod1)
abline(a=-0.7240586 , b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod1.ci <-gls.ci(PARIS$MedianCerebellum, PARIS$MedianCerebrum,vcv(tree_paris))
lines(pGLS.mod1.ci$CI.plot$X,pGLS.mod1.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod1.ci$CI.plot$X,pGLS.mod1.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod1.pi <-gls.pi(PARIS$MedianCerebellum, PARIS$MedianCerebrum,vcv(tree_paris), 1)
lines(pGLS.mod1.pi$PI.plot$X,pGLS.mod1.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod1.pi$PI.plot$X,pGLS.mod1.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()

# Model 2: Cerebellar volume regressed on cerebral volume in non-MNHN data.
tiff("./4.vertebrate_collection/4.10.BM-PGLS-cerebellumCerebrum-noMNHN-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod2 <- gls(MedianCerebellum ~ MedianCerebrum, data=NOPARIS, correlation=bm.noparis)
plot(mod2)
plot(NOPARIS$MedianCerebellum ~ NOPARIS$MedianCerebrum, xlab = "Cerebral volume (log10-transformed)", ylab = "Cerebellar volume (log10-transformed)")
abline(mod2)
abline(a=-0.2519140, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod2.ci <-gls.ci(NOPARIS$MedianCerebellum, NOPARIS$MedianCerebrum,vcv(tree_noparis))
lines(pGLS.mod2.ci$CI.plot$X,pGLS.mod2.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod2.ci$CI.plot$X,pGLS.mod2.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod2.pi <-gls.pi(NOPARIS$MedianCerebellum, NOPARIS$MedianCerebrum,vcv(tree_noparis), 1)
lines(pGLS.mod2.pi$PI.plot$X,pGLS.mod2.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod2.pi$PI.plot$X,pGLS.mod2.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()


