############################################
### Analyses carried out in response
### to Reviewer #3. Specifically, these
### analyses seek to explore ape versus non-
### ape scaling and how these groupings
### interact with the main PGLS by
### dummy coding these groups.
############################################

## Load packages
library(ape)
#library(geiger)
library(nlme)
library(phytools)
library(corrplot)
library(tidyverse)
library(magrittr)
library(evomap)
#library(car) #maybe not necessary (trying some things)
library(rr2) # For PGLS R2, also see R2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs from Anthony Ives.
#library(phylolm)
#library(relaimpo)
library(ggpmisc)

## Set working directory
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/RevisionCode") ## Provide path to folder with scripts and input data.

## Set output directory
outputdir <- "2.ape_pANCOVA"
dir.create(outputdir)


## Load the data and tree.
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
pheno <- read.csv2("./input/Phenotypes_outliersRM.csv", sep=';')  %>% dplyr::rename(species = WilsonReederName) %>% dplyr::select(species, Hominoidea, CerebellarVol, SurfaceArea, CerebralVol, AbsGI, FoldingLength, FoldingNumber, Lambda, Delta, CrusVol_1_3, ROC)
phenoMedian <- read.csv2("./input/PhenotypesMedian_outliersRM.csv", sep=';')  %>% dplyr::rename(species = WilsonReederName) %>% dplyr::select(species, Hominoidea, MedianCerebellum, MedianCerebrum, MedianCrus, MedianROC)

## Match the tree and data, and split the data into ape and  nonape subsets
rownames(phenoMedian) <- tree[["tip.label"]]
phenoMedianAA <- phenoMedian %>% drop_na()
tree<-treedata(tree,phenoMedian, sort=T,warnings=T)$phy # Match tree to the data. Should not change anything.
tree13<-treedata(tree,phenoMedianAA, sort=T,warnings=T)$phy # Match tree to the data. Should reduce the number of species down to 13.


#############################################################################
### We first reaxmine the cerebellar-cerebral PGLS.
#############################################################################

#One could manually assign groupings and filter based on that.
## Nonapes
NONAPES <-phenoMedian %>% filter(Hominoidea== "Nonape")
tree_nonapes <- treedata(tree,NONAPES, sort=T,warnings=T)$phy # Match tree to the nonape data

## Apes
APES <-phenoMedian %>% filter(Hominoidea== "Ape")
tree_apes <- treedata(tree,APES, sort=T,warnings=T)$phy # Match tree to the ape data

# This is an alternative way to subset the trees :). Very convenient
apes <-getTips(tree,findMRCA(tree,c("Pongo_pygmaeus","Hylobates_lar")))
nonapes<-setdiff(getTips(tree,findMRCA(tree,c("Galago_demidoff","Cercopithecus_cephus_cephus"))),apes)


###########################################
### PGLS and PIC regressions for apes
### and non-apes separately.
###########################################

#####################################
### FIRST: nonapes
#####################################

Cnonapes <- NONAPES$MedianCerebellum
Vnonapes <- NONAPES$MedianCerebrum

# Compute contrasts
# Full
pic.C.nonapes <- pic(Cnonapes,tree_nonapes)
pic.V.nonapes <- pic(Vnonapes, tree_nonapes)

## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.
bm<-corBrownian(1,tree_nonapes)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=NONAPES,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.C.nonapes ~pic.V.nonapes + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("2.ape_pANCOVA/2.1.PGLS-PIC-C-V.nonape.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


#####################################
### THEN: apes
#####################################

Capes <- APES$MedianCerebellum
Vapes<- APES$MedianCerebrum

# Compute contrasts
# Full
pic.C.apes <- pic(Capes,tree_apes)
pic.V.apes <- pic(Vapes, tree_apes)

## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.
bm<-corBrownian(1,tree_apes)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=APES,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.C.apes ~pic.V.apes + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("2.ape_pANCOVA/2.2.PGLS-PIC-C-V.ape.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


#############################################################################
### We then re-examine the ansiform area regressions.
#############################################################################

#One could manually assign groupings and filter based on that.
## Nonapes
NONAPES13 <-phenoMedianAA %>% filter(Hominoidea== "Nonape")
tree_nonapes13 <- treedata(tree13,NONAPES13, sort=T,warnings=T)$phy # Match tree to the nonape data

## Apes
APES13 <-phenoMedianAA %>% filter(Hominoidea== "Ape")
tree_apes13 <- treedata(tree13,APES13, sort=T,warnings=T)$phy # Match tree to the ape data

# This is an alternative way to subset the trees :). Very convenient
apes13 <-getTips(tree13,findMRCA(tree13,c("Pongo_pygmaeus","Gorilla_gorilla_gorilla")))
nonapes13 <-setdiff(getTips(tree13,findMRCA(tree13,c("Daubentonia_madagascariensis","Chlorocebus_sabaeus"))),apes13)

#####################################
### FIRST: nonapes
#####################################

## Ansiform areas are regressed on rest of cerebellar (ROC) volumes and cerebral volumes.
AA13nonapes <- NONAPES13$MedianCrus
ROC13nonapes <- NONAPES13$MedianROC
V13nonapes <- NONAPES13$MedianCerebrum

# Compute contrasts
pic.AA13.nonapes <- pic(AA13nonapes, tree_nonapes13)
pic.ROC13.nonapes <- pic(ROC13nonapes, tree_nonapes13)
pic.V13.nonapes <- pic(V13nonapes, tree_nonapes13)


## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.
bm<-corBrownian(1,tree_nonapes13)
pglsModel<-gls(MedianCrus~MedianROC,data=NONAPES13,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.AA13.nonapes ~pic.ROC13.nonapes + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("2.ape_pANCOVA/2.3.PGLS-PIC-AA-ROC.nonape.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.
bm<-corBrownian(1,tree_nonapes13)
pglsModel<-gls(MedianCrus~MedianCerebrum,data=NONAPES13,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.AA13.nonapes ~pic.V13.nonapes + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("2.ape_pANCOVA/2.4.PGLS-PIC-AA-V.nonape.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


#####################################
### THEN: apes
#####################################

## Ansiform areas are regressed on rest of cerebellar (ROC) volumes and cerebral volumes.
AA13apes <- APES13$MedianCrus
ROC13apes <- APES13$MedianROC
V13apes <- APES13$MedianCerebrum

# Compute contrasts
pic.AA13.apes <- pic(AA13apes, tree_apes13)
pic.ROC13.apes <- pic(ROC13apes, tree_apes13)
pic.V13.apes <- pic(V13apes, tree_apes13)


## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.
bm<-corBrownian(1,tree_apes13)
pglsModel<-gls(MedianCrus~MedianROC,data=APES13,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.AA13.apes ~pic.ROC13.apes + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("2.ape_pANCOVA/2.5.PGLS-PIC-AA-ROC.ape.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.
bm<-corBrownian(1,tree_apes13)
pglsModel<-gls(MedianCrus~MedianCerebrum,data=APES13,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.AA13.apes ~pic.V13.apes + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("2.ape_pANCOVA/2.6.PGLS-PIC-AA-V.ape.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


###########################################################
### Use pANCOVA to test deviations in intercept or slope
### and with that deviation from the general allometric 
### relationship for cerebellar-to-cerebral and cerebellar
### & cerebral-to-body mass scaling.
###########################################################


library(evomap)
library(geiger)

#For differences in slope: 
grpS<-rep("A",length(rownames(phenoMedian))) 
grpS[nonapes]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(phenoMedian)

#For differences in intercept: 
grpI<-rep("A",length(rownames(phenoMedian))) 
grpI[nonapes]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(phenoMedian)

## Data subsetting
## Cerebellum vs. Cerebrum
data<-primateData; tree<-primateTree
Y<-"Brain"; X<-"Body"
data<-data[,c(which(colnames(data)==Y),which(colnames(data)==X)),drop=F]
data<-na.omit(data)
# tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)   #match the data to the tree
colnames(data)<-c("Dependent","Independent")  
## Very hacky but bypasses weird error.
data <- data[1:34,]
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend") #Reassign the tree variable. Important.
tree <-treedata(tree,phenoMedian, sort=T,warnings=T)$phy
rownames(data) <- tree[["tip.label"]]

## Plot separate PGLS regressions for apes and nonapes
## Cerebellum vs. Cerebrum

## Again, extremely hacky, but works.
data$Dependent <- phenoMedian$MedianCerebellum
data$Independent <- phenoMedian$MedianCerebrum

pdf("2.ape_pANCOVA/2.7.PGLS.Ape-membership.CC.pdf")
plot(data$Dependent~data$Independent,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=apes,col="lightskyblue",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=nonapes,col="ivory3",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.CC <-model.matrix(as.formula(Dependent~Independent),data)
Model.CC.s <-model.matrix(as.formula(Dependent~grpS:Independent),data) #slopes differ
Model.CC.i <-model.matrix(as.formula(Dependent~grpI + Independent),data) #intercepts differ
Model.CC.si <-model.matrix(as.formula(Dependent~grpI + grpS:Independent),data) #slopes & intercepts differ

sink("2.ape_pANCOVA/2.8.ModelTesting.SlopesIntercepts.CC.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.si)
sink()



######################################
### We next re-examine the ansiform-
### cerebellar and cerebral PGLS.
#####################################

#One could manually reassign groupings and filter based on that.
## Nonapes
#NONAPES13 <-phenoMedianAA %>% filter(Hominoidea== "Nonape")
#tree_nonapes13 <- treedata(tree13,NONAPES13, sort=T,warnings=T)$phy # Match tree to the nonape data

## Apes
#APES13 <-phenoMedianAA %>% filter(Hominoidea== "Ape")
#tree_apes13 <- treedata(tree13,APES13, sort=T,warnings=T)$phy # Match tree to the ape data

# This is an alternative way to subset the trees :). Very convenient
#apes13 <-getTips(tree13,findMRCA(tree13,c("Pongo_pygmaeus","Gorilla_gorilla_gorilla")))
#nonapes13 <-setdiff(getTips(tree13,findMRCA(tree13,c("Daubentonia_madagascariensis","Chlorocebus_sabaeus"))),apes13)


## Ansiform area vs. rest of cerebellum (ROC) 
## Data subsetting
## Ansiform area vs. rest rest of cerebellum & cerebrum
data<-primateData; tree = primateTree
Y<-"Brain"; X<-"Body"
data<-data[,c(which(colnames(data)==Y),which(colnames(data)==X)),drop=F]
data<-na.omit(data)
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)   #match the data to the tree
colnames(data)<-c("Dependent","Independent")  
## Very hacky but bypasses weird error.
data <- data[1:13,]
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend") #Reassign the tree variable. Important.
tree13<-treedata(tree,phenoMedianAA, sort=T,warnings=T)$phy # Match tree to the data. Should reduce the number of species down to 13.
rownames(data) <- tree13[["tip.label"]]

## Reassign the groups, based on smaller data (only 13 species with complete data)
#For differences in slope: 
grpS<-rep("A",length(rownames(phenoMedianAA))) 
grpS[nonapes13]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(phenoMedianAA)

#For differences in intercept: 
grpI<-rep("A",length(rownames(phenoMedianAA))) 
grpI[nonapes13]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(phenoMedianAA)

## Plot separate PGLS regressions for apes and nonapes

#####################################
### FIRST: ansiform area vs. ROC
#####################################

## Again, extremely hacky, but works.
data$Dependent <- phenoMedianAA$MedianCrus
data$Independent <- phenoMedianAA$MedianROC

pdf("2.ape_pANCOVA/2.9.PGLS.Ape-membership.AAROC.pdf")
plot(data$Dependent~data$Independent,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent","Independent",data,tree13,model="BM",group=apes13,col="lightskyblue",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent","Independent",data,tree13,model="BM",group=nonapes13,col="ivory3",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.AAROC<-model.matrix(as.formula(Dependent~Independent),data)
Model.AAROC.s <-model.matrix(as.formula(Dependent~grpS:Independent),data) #slopes differ
Model.AAROC.i <-model.matrix(as.formula(Dependent~grpI + Independent),data) #intercepts differ
Model.AAROC.si <-model.matrix(as.formula(Dependent~grpI + grpS:Independent),data) #slopes & intercepts differ

sink("2.ape_pANCOVA/2.10.ModelTesting.SlopesIntercepts.AAROC.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent~Independent,vcv(tree13),Model.AAROC,Model.AAROC.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent~Independent,vcv(tree13),Model.AAROC,Model.AAROC.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent~Independent,vcv(tree13),Model.AAROC,Model.AAROC.si)
sink()

#######################################
### THEN: ansiform area vs. cerebrum
#######################################

## Again, extremely hacky, but works.
data$Dependent <- phenoMedianAA$MedianCrus
data$Independent <- phenoMedianAA$MedianCerebrum

pdf("2.ape_pANCOVA/2.11.PGLS.Ape-membership.pdf")
plot(data$Dependent~data$Independent,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent","Independent",data,tree13,model="BM",group=apes13,col="lightskyblue",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent","Independent",data,tree13,model="BM",group=nonapes13,col="ivory3",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.AAV<-model.matrix(as.formula(Dependent~Independent),data)
Model.AAV.s <-model.matrix(as.formula(Dependent~grpS:Independent),data) #slopes differ
Model.AAV.i <-model.matrix(as.formula(Dependent~grpI + Independent),data) #intercepts differ
Model.AAV.si <-model.matrix(as.formula(Dependent~grpI + grpS:Independent),data) #slopes & intercepts differ

sink("2.ape_pANCOVA/2.12.ModelTesting.SlopesIntercepts.AAV.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent~Independent,vcv(tree13),Model.AAV,Model.AAV.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent~Independent,vcv(tree13),Model.AAV,Model.AAV.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent~Independent,vcv(tree13),Model.AAV,Model.AAV.si)
sink()


#################################################
### We will also regress cerebellar scaling,
### accounting for cerebral scaling (by nor-
### malizing cerebellar volumes against
### cerebral volumes) in haplorhine species
### only to be more comparable to the study
### of R. Barton and C. Venditti (2014;
### https://doi.org/10.1016/j.cub.2014.08.056).
### Our results in 2.7 and 2.8 illustrate that
### for the whole primate sample (including strep-
### sirrhines), the difference in intercepts in
### ape versus nonape cerebello-cerebral scaling 
### is not replicated.
###############################################

library(evomap)
library(ape)

## Input
## Change csv location to somewhere consistent :) ###
Haplorhine <- read.csv2("~/Downloads/PhenotypesMedian_split_strepsirr.csv") %>% dplyr::rename(species = WilsonReederName) %>% dplyr::rename(Mass = Body.mass.species.mean) %>% select(MedianCerebellum, MedianCerebrum, Mass, Strepsirrhine, Hominoidea) %>% as.data.frame()
Haplorhine$Strepsirrhine <- as.factor(Haplorhine$Strepsirrhine)
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
row.names(Haplorhine) <- tree[["tip.label"]]
# Haplorhine <- na.omit(Haplorhine)
# tree<-treedata(tree,Haplorhine, sort=T,warnings=T)$phy # Match tree to the data
Haplorhine <-Haplorhine %>% filter(Strepsirrhine== "No") %>% select(MedianCerebellum, MedianCerebrum, Hominoidea) # One could also examine the relationship with Mass (not of interest now). In that case, add selection on the 'Mass' column, and make sure to filter the NA.
tree25 <-treedata(tree,Haplorhine, sort=T,warnings=T)$phy # Match tree to the data

## Subsetting
Haplorhine.nonape <- Haplorhine %>% filter(Hominoidea== "Nonape")
Cnonapes <- Haplorhine.nonape$MedianCerebellum
Vnonapes <-Haplorhine.nonape$MedianCerebrum

Haplorhine.ape <- Haplorhine %>% filter(Hominoidea== "Ape")
Capes <- Haplorhine.ape$MedianCerebellum
Vapes <- Haplorhine.ape$MedianCerebrum

#One could manually assign groupings and filter based on that.
tree_nonapes25 <- treedata(tree25,Haplorhine.nonape, sort=T,warnings=T)$phy # Match tree to the nonape data

## Apes
tree_apes25 <- treedata(tree25,Haplorhine.ape, sort=T,warnings=T)$phy # Match tree to the ape data

# This is an alternative way to subset the trees :). Very convenient
apes25 <-getTips(tree25,findMRCA(tree25,c("Pongo_pygmaeus","Hylobates_lar")))
nonapes25 <-setdiff(getTips(tree25,findMRCA(tree25,c("Saimiri_sciureus","Cercopithecus_cephus_cephus"))),apes25)

# Compute contrasts
pic.C.nonapes <- pic(Cnonapes, tree_nonapes25)
pic.V.nonapes <- pic(Vnonapes, tree_nonapes25)

pic.C.nonapes <- pic(Capes, tree_apes25)
pic.V.nonapes <- pic(Vapes, tree_apes25)

#####################################
### FIRST: nonapes
#####################################

## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.
bm<-corBrownian(1,tree25)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=Hapl_nonapes,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.C.nonapes ~pic.V.nonapes + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("2.ape_pANCOVA/2.13.PGLS-PIC-C-V.nonape.haplorhine-only.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

#####################################
### THEN: apes
#####################################

## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.
bm<-corBrownian(1,tree25)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=Hapl_apes,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.C.apes ~pic.V.apes + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("2.ape_pANCOVA/2.14.PGLS-PIC-C-V.ape.haplorhine-only.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


###########################################################
### Use pANCOVA to test deviations in intercept or slope
### and with that deviation from the general allometric 
### relationship for cerebellar-to-cerebral scaling.
###########################################################

## Data subsetting
## Cerebellum vs. Cerebrum
data<-primateData; tree<-primateTree
Y<-"Brain"; X<-"Body"
data<-data[,c(which(colnames(data)==Y),which(colnames(data)==X)),drop=F]
data<-na.omit(data)
# tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)   #match the data to the tree
colnames(data)<-c("Dependent","Independent")  
## Very hacky but bypasses weird error.
data <- data[1:25,]
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend") #Reassign the tree variable. Important.
tree <-treedata(tree,Haplorhine, sort=T,warnings=T)$phy
rownames(data) <- tree[["tip.label"]]

## Again, extremely hacky, but works.
data$Dependent <- Haplorhine$MedianCerebellum
data$Independent <- Haplorhine$MedianCerebrum

pdf("2.ape_pANCOVA/2.15.PGLS.Ape-membership.CC.haplorhine-only.pdf")
plot(data$Dependent~data$Independent,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent","Independent",data,tree25,model="BM",group=apes25,col="lightskyblue",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent","Independent",data,tree25,model="BM",group=nonapes25,col="ivory3",lwd=5,cex=1,pch=19)
graphics.off()

## Reassign the groups, based on smaller data (only 25 species when filtering out strepsirrhines)
#For differences in slope: 
grpS<-rep("A",length(rownames(Haplorhine))) 
grpS[nonapes25]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(Haplorhine)

#For differences in intercept: 
grpI<-rep("A",length(rownames(Haplorhine))) 
grpI[nonapes25]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(Haplorhine)

## Model construction
Model.CC <-model.matrix(as.formula(Dependent~Independent),data)
Model.CC.s <-model.matrix(as.formula(Dependent~grpS:Independent),data) #slopes differ
Model.CC.i <-model.matrix(as.formula(Dependent~grpI + Independent),data) #intercepts differ
Model.CC.si <-model.matrix(as.formula(Dependent~grpI + grpS:Independent),data) #slopes & intercepts differ

sink("2.ape_pANCOVA/2.16.ModelTesting.SlopesIntercepts.CC.haplorhine-only.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.si)
sink()


################################
### END OF SCRIPT
################################

